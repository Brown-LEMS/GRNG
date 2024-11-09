#ifndef IndexDistanceStruct_hpp
#define IndexDistanceStruct_hpp

#include <cmath>

/**
 * @brief struct to hold distance from query to any other index. used for ranked lists
 * 
 */
struct IndexDistanceStruct {
  public:
    // Constructor/Destructor
    IndexDistanceStruct(){};
    IndexDistanceStruct(unsigned int index, float distance) : index(index), distance(distance){};
    ~IndexDistanceStruct(){};

    unsigned int index = 0;
    float distance = -HUGE_VAL;

    inline bool operator<(IndexDistanceStruct const& ref) {
        return (distance < ref.distance);
    }
    inline bool operator>(IndexDistanceStruct const& ref) {
        return (distance > ref.distance);
    }
};

#endif  // IndexDistanceStruct_hpp