#ifndef STREAMHIST_UTILS_H
#define STREAMHIST_UTILS_H

/// Some useful utility functions and classes.

#include <algorithm>
#include <cmath>
#include <vector>

namespace streamhist {
namespace utils {

template <typename Bin>
inline static constexpr double diff(const Bin& a, const Bin& b, bool weighted) noexcept {
    double diff(b.value - a.value);
    if (weighted) {
        diff *= std::log(exp(1) + std::min(a.count, b.count));
    }
    return diff;
}

template <typename List>
inline static std::vector<double> bin_diff(const List& list, bool weighted = false) noexcept {
    size_t l = list.size();
    if (l < 2) {
        return {};
    }

    std::vector<double> r;
    r.reserve(l - 1);
    for (size_t i = 1; i < l; i++) {
        r.push_back(diff(list[i - 1], list[i], weighted));
    }
    return r;
}

template <typename List>
inline static size_t argmin(const List& list) noexcept {
    // Turns out Python's min and max functions are super fast!
    // http://lemire.me/blog/archives/2008/12/17/fast-argmax-in-python/

    size_t l = list.size();
    double m = (std::numeric_limits<double>::max)();

    size_t r = static_cast<size_t>(-1);
    for (size_t i = 0; i < l; i++) {
        if (list[i] < m) {
            r = i;
            m = list[i];
        }
    }
    return r;
}

template <typename Iterator>
inline static auto argmin_bin_diff(Iterator begin, Iterator end, bool weighted = false) noexcept {
    // same as argmin(bin_diff(bins, weighted));

    double m = (std::numeric_limits<double>::max)();

    auto it1(begin);
    auto it2(it1);
    ++it2;

    if (it1 == end) {
        return begin;
    }

    auto mit(begin);

    while (it2 != end) {
        auto v(diff(*it1, *it2, weighted));

        if (v < m) {
            mit = it1;
            m   = v;
        }

        ++it1;
        ++it2;
    }
    return mit;
}

template <typename Iterator>
inline static auto prepare_sums(Iterator begin, Iterator end) noexcept {

    std::vector<double> res;
    if (begin == end) {
        return res;
    }

    res.reserve(end - begin);

    auto it(begin);


    size_t sum = 0;
    while (it != end) {
        auto& b(*it);
        sum += b.count;

        res.push_back(static_cast<double>(sum) - (b.count > 1 ? static_cast<double>(b.count) * 0.5 : 0.));

        ++it;
    }
    return res;
}

/// Custom version of numpy's linspace to avoid numpy depenency.
inline static std::vector<double> linspace(double start, double stop, size_t num) noexcept {
    if (num == 1) {
        return { stop };
    }
    double h = (stop - start) / static_cast<double>(num);

    std::vector<double> r;
    r.reserve(num + 1);
    for (size_t i = 0; i < num + 1; i++) {
        r.push_back(start + h * i);
    }
    return r;
}

/// Super simple quadratic solver.
inline static constexpr bool roots(double a, double b, double c, double& r1, double& r2) noexcept {
    double d = b * b - (4.0 * a * c);
    if (d < 0) {
        return false;
    } else if (d == 0.) {
        r1 = (-b + sqrt(d)) / (2.0 * a);
        r2 = r1;
        return true;
    } else {
        r1 = (-b + sqrt(d)) / (2.0 * a);
        r2 = (-b - sqrt(d)) / (2.0 * a);
        return true;
    }
}

/// Weighted average
inline static constexpr double combine(double v1, size_t c1, double v2, size_t c2) noexcept {
    auto value = (v1 * static_cast<double>(c1) + v2 * static_cast<double>(c2));
    value /= static_cast<double>(c1 + c2);
    return value;
}

/// Weighted average
template <typename ValueType>
inline static constexpr ValueType combine(const ValueType& v1, size_t c1, const ValueType& v2, size_t c2) noexcept {
    return v1.combine(c1, v2, c2);
}

}  // namespace utils
}  // namespace streamhist

#endif
