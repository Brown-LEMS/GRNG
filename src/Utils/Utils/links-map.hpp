#ifndef LinksMap_hpp
#define LinksMap_hpp

#include <tsl/sparse_map.h>
#include <tsl/sparse_set.h>

/**
 * @brief struct to hold links. used in stage 7 for storing links to be removed, and to not check links twice
 * 
 */
struct LinksMap
{
    public:
        LinksMap() {};
        ~LinksMap() {};

        tsl::sparse_map<unsigned int,tsl::sparse_set<unsigned int>> map = 
                tsl::sparse_map<unsigned int,tsl::sparse_set<unsigned int>> {};
        
        // check if link(index1,index2) is already stored. if so,returns true
        inline bool checkAndAddLink(unsigned int const index1, unsigned int const index2) 
        {
            bool flag_alreadyStored = false;

            // store in larger index
            if (index1 < index2) {
                if (map.find(index1) != map.end()) { 
                    if (map[index1].find(index2) != map[index1].end()) {
                        flag_alreadyStored = true;
                    } else {
                        map[index1].insert(index2);
                    }
                } else {
                    map[index1] = tsl::sparse_set<unsigned int> { index2 };
                }
            } else {        
                if (map.find(index2) != map.end()) { 
                    if (map[index2].find(index1) != map[index2].end()) {
                        flag_alreadyStored = true;
                    } else {
                        map[index2].insert(index1);
                    }
                } else {
                    map[index2] = tsl::sparse_set<unsigned int> { index1 };
                }
            }
            
            return flag_alreadyStored;
        }

        /**
         * @brief for coarse neighbor invalidation. order matters. repeat doesn't happen.
         * 
         * @param index 
         * @param indexLink 
         */
        inline void addOrderedIndexLink(unsigned int const index, unsigned int const indexLink) {
            map[index].insert(indexLink);
        }

        inline void addIndexLinks(unsigned int const index, tsl::sparse_set<unsigned int> const& links) {
            map[index].insert(links.begin(),links.end());
        }

        /**
         * @brief add links in LinksMap ref to this instance
         * 
         * @param ref 
         */
        inline void merge(LinksMap const& ref) 
        {
            tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>::const_iterator it1;
            for (it1 = ref.map.begin(); it1 != ref.map.end(); it1++) {
                unsigned int const index1 = (*it1).first;
                tsl::sparse_set<unsigned int> const& links = (*it1).second;

                map[index1].insert(links.begin(),links.end());
            }
        }

        inline bool isIdentical(LinksMap& ref)
        {
            bool flag_identical = true;
            if (map.size() != ref.map.size()) { // same number of maps?
                flag_identical = false;

            } else { // same exact links within?
                tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>::const_iterator it1;
                for (it1 = map.begin(); it1 != map.end(); it1++) {
                    unsigned int const index_this = (*it1).first; 

                    // check same index exists in both
                    if (ref.map.find(index_this) == ref.map.end()) {
                        flag_identical = false;
                        break;
                    }

                    // check same links for index
                    tsl::sparse_set<unsigned int> const& links_this = (*it1).second;
                    tsl::sparse_set<unsigned int> const& links_ref = ref.map[index_this];

                    if (links_this != links_ref) {
                        flag_identical = false;
                        break;
                    }
                }
            }

            return flag_identical;
        }

        inline void clear() {
            map.clear();
        }
};


#endif // LinksMap_hpp