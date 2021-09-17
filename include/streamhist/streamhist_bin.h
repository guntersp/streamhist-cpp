#ifndef STREAMHIST_BIN_H
#define STREAMHIST_BIN_H

#include "streamhist/streamhist_utils.h"

#include <sstream>

namespace streamhist {

/**
 * Histogram bin object.
 *
 * This class implements a simple (value, count) histogram bin pair with
 * several added features such as the ability to merge two bins, comparison
 * methods, and the ability to export and import from dictionaries. The Bin
 * class should be used in conjunction with the StreamHist.
 */
template <typename ValueType>
struct Bin {
    using Type = ValueType;

    ValueType value {};
    size_t    count = 1;


    /**
     * Create a Bin with a given mean and count.
     *
     * Parameters
     * ----------
     * value : float
     *     The mean of the bin.
     * count : int (default=1)
     *     The number of points in this bin. It is assumed that there are
     *      `count` points surrounding `value`, of which `count/2` points are
     *      to the left and `count/2` points are to the right.
     */
    inline constexpr Bin(ValueType _value, size_t _count = 1) noexcept
        : value(_value)
        , count(_count) {}


    /**
     * Create a bin instance from a dictionary.
     *
     * Parameters
     * ----------
     * d : dict
     *     The dictionary must at a minimum a `mean` or `value` key. In
     *     addition, it may contain a `count` key which contains the number
     *     of points in the bin.
     */
    template <typename Map>
    inline static constexpr Bin from_dict(const Map& d) noexcept {
        auto it(d.find("mean"));
        if (it == d.end()) {
            it = d.find("value");
        }

        if (it == d.end()) {
            return {};
        }

        auto it2(d.find("count"));
        return { *it, it2 == d.end() ? 1 : *it2 };
    }

    /**
     * Alternative method for getting the bin's mean and count.
     *
     * Parameters
     * ----------
     * index : int
     *     The index must be either 0 or 1, where 0 gets the mean (value),
     *     and 1 gets the count.
     */
    template <size_t index>
    inline constexpr auto getitem() const noexcept {
        if constexpr (index == 0) {
            return value;
        } else if constexpr (index == 1) {
            return count;
        } else {
            static_assert(index <= 1);
        }
    }

    /**
     * Simple representation of a histogram bin.
     *
     * Returns
     * -------
     * Bin(value=`value`, count=`count`) where value and count are the bin's
     * stored mean and count.
     */
    template <typename Stream>
    inline constexpr Stream& repr(Stream& s) const noexcept {
        s << "Bin(value=" << value << ", count=" << count << ")";
        return s;
    }

    /**
     * Simple representation of a histogram bin.
     *
     * Returns
     * -------
     * Bin(value=`value`, count=`count`) where value and count are the bin's
     * stored mean and count.
     */
    inline std::string repr() const noexcept {
        std::stringstream s;
        repr(s);
        return s.str();
    }

    /**
     * String representation of a histogram bin.
     */
    template <typename Stream>
    inline constexpr Stream& str(Stream& s) const noexcept {
        s << "{'mean': " << value << ", 'count': " << count << "}";
        return s;
    }

    /**
     * String representation of a histogram bin.
     */
    inline std::string str() const noexcept {
        std::stringstream s;
        str(s);
        return s.str();
    }

    /**
     * Tests for equality of two bins.
     *
     * Parameters
     * ----------
     * obj : Bin
     *     The bin to which this bin's mean is compared.
     */
    inline constexpr bool operator==(const Bin& o) const noexcept { return count == o.count && std::abs(value - o.value) < 1e-5; }

    /**
     * Tests if this bin has a lower mean than another bin.
     *
     * Parameters
     * ----------
     * obj : Bin
     *     The bin to which this bin's mean is compared.
     */
    inline constexpr bool operator<(const Bin& o) const noexcept { return value < o.value; }

    /**
     * Tests if this bin has a higher mean than another bin.
     *
     * Parameters
     * ----------
     * obj : Bin
     *     The bin to which this bin's mean is compared.
     */
    inline constexpr bool operator>(const Bin& o) const noexcept { return value > o.value; }

    /**
     * Merge this bin with another bin and return the result.
     *
     * This method implements Step 7 from Algorithm 1 (Update) in ref [1].
     *
     * Parameters
     * ----------
     * obj : Bin
     *     The bin that will be merged with this bin.
     */
    inline constexpr Bin operator+(const Bin& o) const noexcept {
        Bin b({}, count + o.count);  // Summed heights
        if (b.count != 0) {
            using streamhist::utils::combine;
            // Weighted average
            b.value = combine(value, count, o.value, o.count);
        }
        return b;
    }

    /**
     * Merge another bin into this one.
     *
     * Parameters
     * ----------
     * obj : Bin
     *     The bin that will be merged into this bin.
     */
    inline constexpr Bin& operator+=(const Bin& o) noexcept {
        *this = operator+(o);
        return *this;
    }


    /**
     * Finding the density starting from the sum.
     *
     * s = p + (1/2 + r - r^2/2)*i + r^2/2*i1
     * r = (x - m) / (m1 - m)
     * s_dx = i - (i1 - i) * (x - m) / (m1 - m)
     */
    inline static constexpr double compute_density(double p, const Bin& bin_i, const Bin& bin_i1) noexcept {
        double b_diff   = p - bin_i.value;
        double p_diff   = bin_i1.value - bin_i.value;
        double bp_ratio = b_diff / p_diff;

        double inner = (static_cast<double>(bin_i1.count) - static_cast<double>(bin_i.count)) * bp_ratio;
        return (static_cast<double>(bin_i.count) + inner) * (1.0 / (bin_i1.value - bin_i.value));
    }

    inline static constexpr double compute_quantile(double x, const Bin& bin_i, const Bin& bin_i1, double prev_sum) noexcept {
        double d = x - prev_sum;
        double a = static_cast<double>(bin_i1.count) - static_cast<double>(bin_i.count);
        if (bin_i1.count == bin_i.count) {
            double offset = d / (static_cast<double>(bin_i.count + bin_i1.count) / 2.0);
            return bin_i.value + (offset * (bin_i1.value - bin_i.value));
        } else {
            double b = 2.0 * static_cast<double>(bin_i.count);
            double c = -2.0 * d;
            double z = {};
            find_z(a, b, c, z);
            return (bin_i.value + (bin_i1.value - bin_i.value) * z);
        }
    }

    inline static constexpr double compute_sum(double x, const Bin& bin_i, const Bin& bin_i1, double prev_sum) noexcept {
        double b_diff   = x - bin_i.value;
        double p_diff   = bin_i1.value - bin_i.value;
        double bp_ratio = b_diff / p_diff;

        double i1Term = 0.5 * bp_ratio * bp_ratio;
        double iTerm  = bp_ratio - i1Term;

        double first = prev_sum + static_cast<double>(bin_i.count) * iTerm;
        double ss    = first + static_cast<double>(bin_i1.count) * i1Term;
        return ss;
    }

    inline static constexpr bool find_z(double a, double b, double c, double& root) noexcept {
        double r2 {};
        if (!utils::roots(a, b, c, root, r2)) {
            return false;
        }

        if (root >= 0. && root <= 1.) {
            return true;
        }

        if (r2 >= 0. && r2 <= 1.) {
            root = r2;
            return true;
        }

        return false;
    }
};

}  // namespace streamhist

#endif
