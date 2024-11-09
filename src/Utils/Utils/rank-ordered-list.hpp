#ifndef RankOrderedList_hpp
#define RankOrderedList_hpp

#include <algorithm>
#include <vector>

#include "index-distance-struct.hpp"

/**
 * @brief ranked lists of indices based on distance to query. uses IndexDistanceStruct
 *
 */
struct RankOrderedList {
  public:
    RankOrderedList(){};
    ~RankOrderedList(){};

    std::vector<IndexDistanceStruct> list = std::vector<IndexDistanceStruct>{};

    // add to back of stack
    inline void add(IndexDistanceStruct& ids) { list.push_back(ids); }
    inline void add(unsigned int const index, float const distance) { list.push_back(IndexDistanceStruct(index, distance)); }

    // adds elements of list to end
    // does not check for duplicate entries --> user's responsility
    inline void insert(RankOrderedList const& ref) { list.insert(list.end(), ref.list.begin(), ref.list.end()); }

    inline void sort() { std::sort(list.begin(), list.end()); }

    inline void resize(unsigned int const size) { list.resize(size); }

    inline void clear() { list.clear(); }
};

#endif  // RankOrderedList_hpp