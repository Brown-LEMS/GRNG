#ifndef IndexPairDistanceStruct_hpp
#define IndexPairDistanceStruct_hpp
/*
Copyright 2022, Brown University, Providence, RI.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose other than its incorporation into a
commercial product or service is hereby granted without fee, provided
that the above copyright notice appear in all copies and that both
that copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
// Cole Foster
// 01-26-2022

/**
 * @brief Simple Struct to hold distance between two indices. Used to hold max link distance.
 *
 */
struct IndexPairDistanceStruct {
  public:
    // Constructor/Destructor
    IndexPairDistanceStruct(){};
    IndexPairDistanceStruct(unsigned int const index1, unsigned int const index2, float const distance) : invalid(false), index1(index1), index2(index2), distance(distance){};
    ~IndexPairDistanceStruct(){};

    bool invalid = true;
    unsigned int index1 = 0;
    unsigned int index2 = 0;
    float distance = 0.0f;

    inline bool operator==(IndexPairDistanceStruct const& ref) const {
        bool same = (((index1 == ref.index1) && (index2 == ref.index2)) ||
                     ((index1 == ref.index2) && (index2 == ref.index1)));
        return same;
    }

    inline bool operator!=(IndexPairDistanceStruct const& ref) const {
        bool same = ((*this) == ref);
        return !same;
    }

    /**
     * @brief determine if instance pair matches this newly given pair of indices
     *
     * @param index1a
     * @param index2a
     * @return true
     * @return false
     */
    inline bool const isLinkMatch(unsigned int const index1a, unsigned int const index2a) const {
        bool isMatch = false;
        if (index1a == index1 && index2a == index2) {
            isMatch = true;
        } else if (index1a == index2 && index2a == index1) {
            isMatch = true;
        }
        return isMatch;
    }

    /**
     * @brief determine if instance pair matches newly given
     *
     * @param ipds
     * @return true
     * @return false
     */
    inline bool const isLinkMatch(IndexPairDistanceStruct const& ipds) const {
        return isLinkMatch(ipds.index1, ipds.index2);
    }

    /**
     * @brief Check if exact match, indices in order. important for coarse neighbors
     * 
     * @param index1a 
     * @param index2a 
     * @return true 
     * @return false 
     */
    inline bool const isExactMatch(unsigned int const index1a, unsigned int const index2a) const {
        bool isMatch = false;
        if (index1a == index1 && index2a == index2) {
            isMatch = true;
        }
        return isMatch;
    }
    
    inline bool const isExactMatch(IndexPairDistanceStruct const& ipds) const {
        return isExactMatch(ipds.index1, ipds.index2);
    }

    /**
     * @brief instance invalid if not initialized with index pair and distance
     *
     * @return true
     * @return false
     */
    inline bool const isInvalid() const { return (invalid == true); }
};

#endif  // IndexPairDistanceStruct_hpp