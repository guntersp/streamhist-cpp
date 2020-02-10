#ifndef STREAMHIST_BIN_CONTAINER_H
#define STREAMHIST_BIN_CONTAINER_H

#include "streamhist/streamhist_bin.h"

#include <algorithm>
#include <vector>

namespace streamhist {
namespace utils {

template <typename _Item>
struct SortedListWithKey {
    using Item = _Item;


    inline void clear() noexcept { list.clear(); }

    inline size_t size() const noexcept { return list.size(); }

    inline void reserve(size_t n) noexcept { list.reserve(n); }

    inline Item&       operator[](size_t index) noexcept { return list[index]; }
    inline const Item& operator[](size_t index) const noexcept { return list[index]; }

    inline void add(const Item& i) noexcept {
        if (list.empty()) {
            list.push_back(i);
        } else {
            list.insert(std::upper_bound(list.begin(), list.end(), i), i);
        }
    }

    inline void add(Item&& i) noexcept {
        if (list.empty()) {
            list.push_back(std::move(i));
        } else {
            list.insert(std::upper_bound(list.begin(), list.end(), i), std::move(i));
        }
    }

    inline auto begin() noexcept { return list.begin(); }
    inline auto begin() const noexcept { return list.begin(); }

    inline auto end() noexcept { return list.end(); }
    inline auto end() const noexcept { return list.end(); }

    inline auto& front() noexcept { return list.front(); }
    inline auto& front() const noexcept { return list.front(); }

    inline auto& back() noexcept { return list.back(); }
    inline auto& back() const noexcept { return list.back(); }

    inline auto find(const Item& i) noexcept {
        auto it(std::lower_bound(list.begin(), list.end(), i));

        if (it != list.end() && !(i < *it)) {
            return it;
        }

        return list.end();
    }

    inline auto find(const Item& i) const noexcept {
        auto it(std::lower_bound(list.begin(), list.end(), i));
        return it != list.end() && !(i < *it) ? it : list.end();
    }

    inline auto pop(typename std::vector<Item>::iterator it) noexcept { return list.erase(it); }

    inline size_t bisect_left(const Item& i) const noexcept { return std::lower_bound(list.begin(), list.end(), i) - list.begin(); }

    inline size_t bisect_right(const Item& i) const noexcept { return std::upper_bound(list.begin(), list.end(), i) - list.begin(); }

    inline size_t bisect(const Item& i) const noexcept { return bisect_right(i); }

    inline constexpr bool operator==(const SortedListWithKey& o) const noexcept { return list == o.list; }

    template <size_t size>
    inline constexpr bool operator==(const Item (&o)[size]) const noexcept {
        if (list.size() != size) {
            return false;
        }

        for (size_t i = 0; i < size; i++) {
            if (!(list[i] == o[i])) {
                return false;
            }
        }
        return true;
    }

private:
    std::vector<Item> list;
};


}  // namespace utils
}  // namespace streamhist

#endif
