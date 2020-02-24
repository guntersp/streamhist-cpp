/**
 * A streaming approximate histogram based on the algorithms found in Ben-Haim &
 * Tom-Tov's Streaming Parallel Decision Tree Algorithm [1]_. Histogram bins do
 * not have a preset size. As values stream into the histogram, bins are
 * dynamically added and merged.
 *
 * The accurate method of calculating quantiles (like percentiles) requires data
 * to be sorted. Streaming histograms make it possible to approximate quantiles
 * without sorting (or even individually storing) values.
 *
 * A maximum bin size is passed as an argument to the init methods. A larger bin
 * size yields more accurate approximations at the cost of increased memory
 * utilization and performance. There is no "optimal" bin count, but somewhere
 * between 20 and 80 bins should be sufficient.
 *
 * The historgram class implented here is based on VividCortex's "Streaming
 * approximate histograms in Go" [2]_, which is released under the MIT Open
 * Source license. Additional algorithmic adjustments and methods were adapted from
 * BigLM's "Streaming Histograms for Clojure/Java" [3]_, which is released under
 * the Apache License, Version 2.0.
 *
 * References
 * ----------
 * .. [1] http://jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf
 * .. [2] https://github.com/VividCortex/gohistogram
 * .. [3] https://github.com/bigmlcom/histogram
 * .. [4] https://vividcortex.com/blog/2013/07/08/streaming-approximate-histograms/
 *
 *
 * Copyright © 2015 Carson Farmer <carsonfarmer@gmail.com>
 * Copyright © 2013 VividCortex
 * All rights reserved. MIT Licensed.
 * Copyright © 2013, 2014, 2015 BigML
 * Licensed under the Apache License, Version 2.0
 */

#ifndef STREAMHIST_H
#define STREAMHIST_H


#include "streamhist/streamhist_bin.h"
#include "streamhist/streamhist_container.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <utility>

namespace streamhist {

/// A StreamHist implementation.
template <typename ValueType>
struct StreamHist {
    using Type = ValueType;
    using Bin  = streamhist::Bin<ValueType>;

    utils::SortedListWithKey<Bin> bins;

    size_t maxbins       = 64;
    size_t total         = 0;
    bool   weighted      = false;
    double _min          = (std::numeric_limits<double>::max)();
    double _max          = -(std::numeric_limits<double>::max)();
    size_t freeze        = static_cast<size_t>(-1);
    size_t missing_count = 0;


    inline constexpr StreamHist(size_t maxbins = 64, bool weighted = false, size_t freeze = static_cast<size_t>(-1)) noexcept
        : maxbins(maxbins)
        , weighted(weighted)
        , freeze(freeze) {}


    inline constexpr bool empty() const noexcept { return total <= 0; }


    /// Add a point to the histogram.
    template <typename List>
    inline constexpr StreamHist& update(const List& n, size_t count = 1) noexcept {
        /**
         * Shortcut for updating a histogram with an iterable
         * This works for anything that supports iteration, including
         * file-like objects and readers
         * This also means that nested lists (and similar structures) will
         * be 'unpacked' and added to the histogram 'automatically'
         */
        for (auto& p : n) {
            update(p, count);  // Count is assumed to apply for all
        }

        return trim();
    }

    /// Add a point to the histogram.
    template <typename Iterator>
    inline constexpr StreamHist& update(Iterator begin, Iterator end, size_t count = 1) noexcept {
        /**
         * Shortcut for updating a histogram with an iterable
         * This works for anything that supports iteration, including
         * file-like objects and readers
         * This also means that nested lists (and similar structures) will
         * be 'unpacked' and added to the histogram 'automatically'
         */
        while (begin != end) {
            update(*begin, count);  // Count is assumed to apply for all
            ++begin;
        }

        return trim();
    }

    /// Add a point to the histogram.
    inline constexpr StreamHist& update(const ValueType& n, size_t count = 1) noexcept {
        if (std::isnan(n)) {
            // We simply keep a count of the number of missing values
            missing_count += count;
            return *this;
        }

        insert(n, count);
        return trim();
    }

    /// Add a point to the histogram.
    inline constexpr StreamHist& updateValues(const ValueType& n) noexcept { return update(n); }

    /// Add a point to the histogram.
    template <typename Container>
    inline constexpr StreamHist& updateValues(const Container& n) noexcept {
        return update(n);
    }

    /// Add a point to the histogram.
    template <typename First, typename... Values>
    inline constexpr StreamHist& updateValues(const First& n, const Values&... vs) noexcept {
        update(n);

        if constexpr (sizeof...(vs) > 0) {
            updateValues(vs...);
        }
        return trim();
    }


    /**
     * Inserts a point to the histogram.
     *
     * This method implements Steps 1-4 from Algorithm 1 (Update) in ref [1].
     *
     * Notes
     * -----
     * It is better to use `update` when inserting data into the histogram,
     * as `insert` does not automatically update the total point count, or
     * call `trim` after the insertion. For large batches of inserts, insert
     * may be more efficient, but you are responsible for updating counts
     * and trimming the bins 'manually'.
     *
     * Examples
     * --------
     * >>> # Using insert
     * >>> h = StreamHist().insert(1).insert(2).insert(3)
     * >>> h.update_total(3)
     * >>> h.trim()
     *
     * >>> # Using update
     * >>> h = StreamHist().update([1, 2, 3])
     */
    inline constexpr StreamHist& insert(const ValueType& n, size_t count = 1) noexcept {
        update_total(count);

        if (_min > n) {
            _min = n;
        }
        if (_max < n) {
            _max = n;
        }

        Bin b { n, count };

        using streamhist::utils::integrate;

        auto it(bins.find(b));
        if (it != bins.end()) {
            auto& b(*it);
            integrate(b.value, b.count, n, count);
            b.count += count;

        } else {
            if (freeze != static_cast<size_t>(-1) && total >= freeze) {
                size_t index = bins.bisect(b);

                double prev_dist = (std::numeric_limits<double>::max)();
                if (index > 0) {
                    prev_dist = n - bins[index - 1].value;
                }

                double next_dist = (std::numeric_limits<double>::max)();
                if (index && index < bins.size()) {
                    next_dist = bins[index].value - n;
                }

                if (prev_dist < next_dist) {
                    auto& b(bins[index - 1]);
                    integrate(b.value, b.count, n, count);
                    b.count += count;

                } else {
                    auto& b(bins[index]);
                    integrate(b.value, b.count, n, count);
                    b.count += count;
                }

            } else {
                bins.add(std::move(b));
            }
        }

        return *this;
    }


    /// Return the value of the cumulative distribution function at x.
    inline constexpr double cdf(double x) const noexcept { return sum(x) / static_cast<double>(total); }

    /// Return the value of the probability density function at x.
    inline constexpr double pdf(double x) const noexcept { return density(x) / static_cast<double>(total); }

    /// Return the upper (max( and lower (min) bounds of the distribution.
    inline constexpr std::pair<double, double> bounds() const noexcept {
        if (len() > 0) {
            return { _min, _max };
        }
        return {};
    }

    /// Return the number of bins in this histogram.
    inline size_t count() const noexcept { return total; }

    /**
     * Return a median for the points inserted into the histogram.
     *
     * This will be the true median whenever the histogram has less than
     * the maximum number of bins, otherwise it will be an approximation.
     */
    inline constexpr double median() const noexcept {
        if (total == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }


        /**
         * bug fix e.g. for:
         * - update([0, 0, 1, 1])
         * -> total = 4
         * -> len(bins) = 2
         * -> index out of range
         */
        // if (bins.size() >= maxbins) {
        if (bins.size() != total) {
            // Return the approximate median
            return quantiles(0.5).front();
        }

        // Return the 'exact' median when possible
        size_t mid = total / 2;
        if (total == mid * 2) {
            return (bins[mid - 1] + bins[mid]).value;
        }
        return bins[mid].value;
    }

    /**
     * Return the sample mean of the distribution.
     */
    inline constexpr double mean() const noexcept {
        if (total == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double s = 0.;  // Sum
        for (auto& b : bins) {
            s += b.value * static_cast<double>(b.count);
        }
        return s / static_cast<double>(total);
    }

    /**
     * Return the variance of the distribution.
     */
    inline constexpr double var() const noexcept {
        if (total < 2) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double s = 0.;      // Sum
        double m = mean();  // Mean
        for (auto& b : bins) {
            double v = b.value - m;
            s += (b.count * v * v);
        }
        return s / static_cast<double>(total);
    }

    /**
     * Return the minimum value in the histogram.
     */
    inline constexpr double(min)() const noexcept { return _min; }

    /**
     * Return the maximum value in the histogram.
     */
    inline constexpr double(max)() const noexcept { return _max; }

    /**
     * Merge adjacent bins to decrease bin count to the maximum value.
     *
     * This method implements Steps 5-6 from Algorithm 1 (Update) in ref [1].
     */
    inline constexpr StreamHist& trim() noexcept {
        while (bins.size() > maxbins) {
            auto it = utils::argmin_bin_diff(bins.begin(), bins.end(), weighted);
            auto bin(*it);
            it = bins.pop(it);
            bin += *it;
            bins.pop(it);
            bins.add(std::move(bin));
        }

        return *this;
    }

    /**
     * Return a string reprentation of the histogram.
     */
    template <typename Stream>
    inline constexpr Stream& str(Stream& s) const noexcept {
        if (bins.size() > 0) {
            s << "Mean\tCount\n----\t-----\n";
            for (auto& b : bins) {
                s << b.value << "\t" << b.count << "\n";
            }
            s << "----\t-----\n";
            s << "Missing values: " << missing_count << "\n";
            s << "Total count: " << total;

        } else {
            s << "Empty histogram";
        }
        return s;
    }

    /**
     * Return a string reprentation of the histogram.
     */
    inline std::string str() const noexcept {
        std::stringstream s;
        str(s);
        return s.str();
    }

    inline constexpr bool operator==(const StreamHist& o) const noexcept {
        if (missing_count != o.missing_count || maxbins != o.maxbins || weighted != o.weighted || freeze != o.freeze) {
            return false;
        }

        return bins == o.bins;
    }


    /**
     * Return the number of bins in this histogram.
     */
    inline constexpr size_t len() const noexcept { return bins.size(); }

    /**
     * Update the internally-stored total number of points.
     */
    inline constexpr void update_total(size_t size = 1) noexcept { total += size; }

    /**
     * Merge two StreamHist objects into one.
     */
    inline constexpr StreamHist operator+(const StreamHist& o) const noexcept {
        StreamHist r(*this);
        r += o;
        return r;
    }

    /**
     * Merge another StreamHist object into this one.
     */
    inline constexpr StreamHist& operator+=(const StreamHist& o) noexcept { return merge(o); }

    /**
     * Merge another StreamHist object into this one.
     *
     * This method implements Algorithm 2 (Merge) in ref [1].
     */
    inline constexpr StreamHist& merge(const StreamHist& o, size_t size = static_cast<size_t>(-1)) noexcept {

        /*if (other == 0) {  // Probably using sum here...
            return *this;  // This is a little hacky...
        }*/

        for (auto& b : o.bins) {
            bins.add(b);
        }
        total += o.total;
        if (size != static_cast<size_t>(-1)) {
            maxbins = size;
        }
        trim();
        _min = std::min(_min, o._min);
        _max = std::max(_max, o._max);
        missing_count += o.missing_count;

        return *this;
    }

private:
    template <typename Map>
    inline void describeInt(Map& r, double x) const noexcept {
        r[std::to_string(int(x * 100.)) + '%'] = quantiles(x).front();
    }

    template <typename Map, typename... Quantile>
    inline void describeInt(Map& r, double x, Quantile... quantiles) const noexcept {
        describeInt(r, x);
        if constexpr (sizeof...(quantiles) > 0) {
            describeInt(r, quantiles...);
        }
    }

public:
    /**
     * Generate various summary statistics.
     */
    inline std::map<std::string, double> describe(double quantiles) const noexcept {
        std::map<std::string, double> r;
        r["count"] = count();
        r["mean"]  = mean();
        r["var"]   = var();
        r["min"]   = min();

        describeInt(r, quantiles);

        r["max"] = max();
        return r;
    }

    /**
     * Generate various summary statistics.
     */
    template <typename... Quantile>
    inline std::map<std::string, double> describe(Quantile... quantiles) const noexcept {
        std::map<std::string, double> r;
        r["count"] = count();
        r["mean"]  = mean();
        r["var"]   = var();
        r["min"]   = min();

        if constexpr (sizeof...(quantiles) > 0) {
            describeInt(r, quantiles...);
        }
        r["max"] = max();
        return r;
    }

    /**
     * Generate various summary statistics.
     */
    template <typename List>
    inline std::map<std::string, double> describe(const List& quantiles) const noexcept {
        std::map<std::string, double> r;
        r["count"] = count();
        r["mean"]  = mean();
        r["var"]   = var();
        r["min"]   = min();

        for (auto x : quantiles) {
            describeInt(r, x);
        }
        r["max"] = max();
        return r;
    }

    struct Breaks {
        std::vector<double> counts;
        std::vector<double> bounds;
    };

    /**
     * Return output like that of numpy.histogram.
     */
    inline Breaks compute_breaks(size_t n = 50) const noexcept {
        if (n <= 0) {
            return {};
        }

        Breaks b;
        b.bounds = utils::linspace(_min, _max, n);
        b.counts.reserve(n);

        double last = 0.;

        size_t num(b.bounds.size());
        for (size_t i = 1; i < num; i++) {
            auto e(b.bounds[i]);

            auto s(sum(e));
            b.counts.push_back(s - last);
            last = s;
        }

        return b;
    }

    /**
     * Print a string reprentation of the histogram.
     */
    template <typename Stream>
    inline constexpr Stream& print_breaks(Stream& s, size_t n = 50) const noexcept {
        Breaks breaks(compute_breaks(n));

        size_t num(std::min(breaks.bounds.size(), breaks.counts.size()));

        for (size_t j = 0; j < num; j++) {
            auto c(breaks.counts[j]);
            auto b(breaks.bounds[j]);

            s << b << '\t';

            size_t w(static_cast<size_t>(c / static_cast<double>(total) * 200.));
            for (size_t i = 0; i < w; i++) {
                s << '.';
            }

            s << '\n';
        }

        return s;
    }

    /**
     * Print a string reprentation of the histogram.
     */
    inline void print_breaks(size_t n = 50) const noexcept { print_breaks(std::cout, n); }


    /**
     * Return the estimated number of points in the interval [−∞, b].
     */
    inline constexpr double sum(double x) const noexcept {
        if (x < _min) {
            return 0.;  // Sum is zero!
        } else if (x >= _max) {
            return static_cast<double>(total);
        } else {
            size_t i        = bins.bisect_right(x) - 1;  // rightmost bin lte x
            auto   bin_i    = i >= 0 ? bins[i] : Bin(_min, 0);
            auto   bin_i1   = i + 1 < bins.size() ? bins[i + 1] : Bin(_max, 0);
            double prev_sum = 0.;
            for (size_t j = 0; j < std::min(i, bins.size()); j++) {
                prev_sum += bins[j].count;
            }
            prev_sum += bin_i.count / 2.;
            return Bin::compute_sum(x, bin_i, bin_i1, prev_sum);
        }
    }

    inline constexpr double density(double p) const noexcept {
        if (p < _min || p > _max) {
            return 0.;
        } else if (p == _min && p == _max) {
            return std::numeric_limits<double>::infinity();
        } else if (bins.find(Bin(p, 0)) != bins.end()) {
            double high = std::nextafter(p, std::numeric_limits<double>::infinity());
            double low  = std::nextafter(p, -std::numeric_limits<double>::infinity());
            return (density(low) + density(high)) / 2.;
        } else {
            auto bin_i  = lower(p);
            auto bin_i1 = higher(p);

            if (!bin_i && !bin_i1) {
                return Bin::compute_density(p, Bin(_min, 0), Bin(_max, 0));
            } else if (!bin_i) {
                return Bin::compute_density(p, Bin(_min, 0), *bin_i1);
            } else if (!bin_i1) {
                return Bin::compute_density(p, *bin_i, Bin(_max, 0));
            }

            return Bin::compute_density(p, *bin_i, *bin_i1);
        }
    }

private:
    template <typename List1, typename List2>
    inline void quantilesInt(List1& result, const List2& sums, double q) const noexcept {
        if (q <= 0.) {
            result.push_back(_min);
            return;
        }

        if (q >= 1.) {
            result.push_back(_max);
            return;
        }

        auto target_sum = q * static_cast<double>(total - 1) + 1;
        auto it         = std::upper_bound(sums.begin(), sums.end(), target_sum);
        auto i(it - sums.begin());

        double l0(_min);
        size_t l1(0);
        double s(target_sum);
        if (it != sums.begin()) {
            auto& b(bins[i - 1]);
            l0 = b.value;
            l1 = b.count;
            s -= sums[i - 1];
        } else {
            s -= 1;
        }

        double r0(_max);
        size_t r1(0);
        if (it != sums.end()) {
            auto& b(bins[i]);
            r0 = b.value;
            r1 = b.count;
        }

        if (l1 <= 1 && r1 <= 1) {
            // We have exact info at this quantile.  Match linear interpolation
            // strategy of numpy.quantile().
            if (r1 > 0) {
                result.push_back(l0 + (r0 - l0) * s / static_cast<double>(r1));
            } else {
                result.push_back(l0);
            }

        } else {
            if (r1 == 1) {
                // For exact bin on RHS, compensate for trapezoid interpolation using
                // only half of count.
                r1 = 2;
            }
            double bp_ratio {};
            if (l1 == r1) {
                bp_ratio = s / static_cast<double>(l1);
            } else {
                auto d(static_cast<double>(l1) - static_cast<double>(r1));
                bp_ratio = (static_cast<double>(l1) - std::sqrt(static_cast<double>(l1 * l1) - 2 * s * d)) / d;
            }
            result.push_back(bp_ratio * (r0 - l0) + l0);
        }
    }

    template <typename List1, typename List2, typename... Quantile>
    inline void quantilesInt(List1& result, const List2& sums, double x, Quantile... quantiles) const noexcept {
        quantilesInt(result, sums, x);
        if constexpr (sizeof...(quantiles) > 0) {
            quantilesInt(result, sums, quantiles...);
        }
    }

public:
    /**
     * Return the estimated data value for the given quantile(s).
     *
     * The requested quantile(s) must be between 0 and 1. Note that even if a
     * single quantile is input, a list is always returned.
     */
    inline auto quantiles(double quantiles) const noexcept {
        auto sums(utils::prepare_sums(bins.begin(), bins.end()));

        std::vector<double> result;
        quantilesInt(result, sums, quantiles);
        return result;
    }

    /**
     * Return the estimated data value for the given quantile(s).
     *
     * The requested quantile(s) must be between 0 and 1. Note that even if a
     * single quantile is input, a list is always returned.
     */
    template <typename... Quantile>
    inline auto quantiles(Quantile... quantiles) const noexcept {
        auto sums(utils::prepare_sums(bins.begin(), bins.end()));

        std::vector<double> result;
        if constexpr (sizeof...(quantiles) > 0) {
            quantilesInt(result, sums, quantiles...);
        }
        return result;
    }

    /**
     * Return the estimated data value for the given quantile(s).
     *
     * The requested quantile(s) must be between 0 and 1. Note that even if a
     * single quantile is input, a list is always returned.
     */
    template <typename List>
    inline auto quantiles(const List& quantiles) const noexcept {
        auto sums(utils::prepare_sums(bins.begin(), bins.end()));

        std::vector<double> result;
        for (auto x : quantiles) {
            quantilesInt(result, sums, x);
        }
        return result;
    }

    inline constexpr const Bin* floor(double p) const noexcept {
        Bin    hbin(p, 0);
        size_t index = bins.bisect_left(hbin);
        if (bins.find(hbin) == bins.end()) {
            index -= 1;
        }
        if (index != static_cast<size_t>(-1)) {
            return &bins[index];
        }
        return nullptr;
    }

    inline constexpr const Bin* ceiling(double p) const noexcept {
        Bin    hbin(p, 0);
        size_t index = bins.bisect_right(hbin);
        if (bins.find(hbin) != bins.end()) {
            index -= 1;
        }
        if (index < bins.size()) {
            return &bins[index];
        }
        return nullptr;
    }

    inline constexpr const Bin* lower(double p) const noexcept {
        size_t index = bins.bisect_left(Bin(p, 0)) - 1;
        if (index != static_cast<size_t>(-1)) {
            return &bins[index];
        }
        return nullptr;
    }

    inline constexpr const Bin* higher(double p) const noexcept {
        size_t index = bins.bisect_right(Bin(p, 0));
        if (index < bins.size()) {
            return &bins[index];
        }
        return nullptr;
    }
};


}  // namespace streamhist

#endif
