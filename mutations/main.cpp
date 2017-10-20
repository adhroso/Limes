//
//  main.cpp
//  Limes with single mutation
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#include "Serial.h"

long offset_test_forward_(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2, std::size_t start_pos, Mutations &mutations) {
    std::pair<std::string::const_iterator, std::string::const_iterator> pair;
    
    Mutations m;
    
    pair = std::mismatch(first1, last1, first2, [&m, &start_pos](char a, char b) {
        if (a=='N' || b=='N') {
            return false;
        } else if (a!=b) {
            m.push_back(start_pos);
        }
        start_pos++;
        return m.size()<2;
    });
    
    m.swap(mutations);
    //return std::distance(first1, pair.first);
    
    return start_pos;
}

long offset_test_reverse_(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, std::size_t start_pos, Mutations& mutations) {
    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
    
    Mutations m;
    
    pair = std::mismatch(first1, last1, first2, [&m,&start_pos](char a, char b) {
        if (a=='N' || b=='N') {
            return false;
        } else if (a!=b) {
            m.push_back(start_pos);
        }
        start_pos--;
        return m.size()<2;
    });
    
    //because we're expanding to the left, second accuring mutation should be in the first position.
    if (m.size() > 1)
        std::swap(m[0], m[1]);
    
    m.swap(mutations);
    
    //return  std::distance(first1, pair.first);
    return start_pos;
}

int compute_range_and_start_pos_test(int left_first_mut, int left_second_mut, int right_first_mut, int right_second_mut, int &start) {
    int max_range = 0;
    
    int range1 = right_first_mut - left_first_mut - 1;
    int range2 = right_second_mut - left_second_mut - 1;
    
    if (range1 >  range2) {
        max_range = range1;
        start = left_first_mut+1;
    } else {
        max_range = range2;
        start = left_second_mut+1;
    }
    
    
    return max_range;
}

void test2() {
    std::string s1 ("AAAAAAAAAGAAAAGAAAAATAAAATAATAAAAAA");
    std::string s2 ("AAACAAAAAAACAAAAAATAAAAAAAAA");
    
    std::size_t idx1 = 24;
    std::size_t idx2 = 18;
    
    //    std::cout << s1[idx1] << std::endl;
    //    std::cout << s2[idx2] << std::endl;
    //    return;
    std::string::const_iterator d1_begin = s1.begin();
    std::string::const_iterator d2_begin = s2.begin();
    std::string::const_iterator d1_end = s1.end();
    std::string::const_iterator d2_end = s2.end();
    
    const Sequence::size_type size1 = s1.length(), size2 = s2.length();
    
    std::string::const_reverse_iterator d1_rbegin = s1.rbegin();
    std::string::const_reverse_iterator d2_rbegin = s2.rbegin();
    std::string::const_reverse_iterator d1_rend = s1.rend();
    std::string::const_reverse_iterator d2_rend = s2.rend();
    
    std::advance(d1_rbegin, size1-idx1);
    std::advance(d2_rbegin, size2-idx2);
    
    std::advance(d1_begin, idx1);
    std::advance(d2_begin, idx2);
    //--------------------------------------
    
    Mutations left_mutations, right_mutations;
    
    //offset_test_forward(d1_begin, d1_end, d2_begin, idx1, m);
    std::size_t right_most_offset = (size1-idx1 < size2-idx2) ? offset_test_forward_(d1_begin, d1_end, d2_begin, idx1, right_mutations) : offset_test_forward_(d2_begin, d2_end, d1_begin, idx1, right_mutations);
    std::size_t left_most_offset = (idx1 < idx2) ? offset_test_reverse_(d1_rbegin, d1_rend, d2_rbegin, idx1, left_mutations) : offset_test_reverse_(d2_rbegin, d2_rend, d1_rbegin, idx1, left_mutations);
    
    //short left_first_mut = 0, left_second_mut = 0;
    std::size_t max_range = 0;
    std::size_t query_start = 0, target_start = 0;
    
    if (left_mutations.size()+right_mutations.size() < 2) {                         //case 1,2,3
        std::cout << "cases 1,2,3" << std::endl;
        query_start = left_most_offset;
        target_start = idx2-(idx1-query_start);
        max_range = static_cast<int>(right_most_offset-left_most_offset-1);
        
    } else if (right_mutations.size() + left_mutations.size() > 3) {                //case 9
        const std::size_t range1 = right_mutations[0] - left_mutations[0];
        const std::size_t range2 = right_mutations[1] - left_mutations[1];
        
        if (range1 >  range2) {
            std::cout << "case 9, range1" << std::endl;
            max_range = range1;
            query_start = left_mutations[0];
            target_start = idx2-(idx1-left_mutations[0]);
            
        } else {
            std::cout << "case 9, range2" << std::endl;
            max_range = range2;
            query_start = left_mutations[1]+1;
            target_start = idx2-(idx1-left_mutations[1]);
        }
    } else if (right_mutations.size() == 1 && left_mutations.size() == 1) {         //case 4
        const int range1 = static_cast<int>(right_most_offset - left_mutations[0]);
        const int range2 = static_cast<int>(right_mutations[0] - left_most_offset);
        
        if (range1 > range2) {
            std::cout << "case 4, region 1" << std::endl;
            max_range = range1;
            query_start = left_mutations[0];
            target_start = target_start = idx2-(idx1-left_mutations[0]);
            
        } else {
            std::cout << "case 4, region 2" << std::endl;
            max_range = range2;
            query_start = left_most_offset;
            target_start = target_start = idx2-(idx1-left_most_offset);
        }
        
    } else if (right_mutations.size() + left_mutations.size() > 2) {                //case 7,8
        
        if (right_mutations.size() > 1) {
            
            const std::size_t range1 = static_cast<int>(right_mutations[0] - left_most_offset);
            const std::size_t range2 = right_mutations[1] - left_mutations[0];
            
            if (range1 > range2) {
                std::cout << "case 7 range1" << std::endl;
                max_range = range1;
                query_start = left_most_offset;
                target_start = target_start = idx2-(idx1-left_most_offset);
                
            } else {
                std::cout << "case 7 range2" << std::endl;
                max_range = range2;
                query_start = left_mutations[0];
                target_start = idx2-(idx1-left_mutations[0]);
            }
            
        } else {
            const std::size_t range1 = right_mutations[0] - left_mutations[0];
            const std::size_t range2 = static_cast<int>(right_most_offset - left_mutations[1]);
        
            if (range1 > range2) {
                std::cout << "case 8 range1" << std::endl;
                max_range = range1;
                query_start = left_mutations[0];
                target_start = target_start = idx2-(idx1-left_mutations[0]);
                
            } else {
                std::cout << "case 8 range2" << std::endl;
                max_range = range2;
                query_start = left_mutations[1];
                target_start = idx2-(idx1-left_mutations[1]);
            }
        }
    } else {                                                                        //case 5,6
        
        if (right_mutations.empty()) {                                              //case 5
            std::cout << "case 5" << std::endl;
            max_range = static_cast<int>(right_most_offset-left_mutations[0]);
            query_start = left_mutations[0];
            target_start = target_start = idx2-(idx1-left_mutations[0]);
            
        } else {                                                                    //case 6
            std::cout << "case 6" << std::endl;
            max_range = static_cast<int>(right_mutations[1] - left_most_offset);
            query_start = left_most_offset;
            target_start = target_start = idx2-(idx1-left_most_offset);
        }
    }
    
    std::cout << "query range: " << max_range << ", start pos: " << query_start << std::endl;
    std::cout << "target range: " << max_range << ",start pos: " << target_start << std::endl;
    
    std::cout << "query:" << s1.substr(query_start, max_range) << std::endl;
    std::cout << "target:" << s2.substr(target_start, max_range) << std::endl;
    
}

int main(int argc, const char * argv[]) {
//    test2();
//    return 0;

    if (argc < 6)  {
        std::cerr << "usage: ./Limes <genome 1 dir> <genome 2 dir> <file_ext> [fa]> <limes file> <progress file>" << std::endl;
        return EXIT_FAILURE;
    } else {
        const std::string dir1 = argv[1];
        const std::string dir2 = argv[2];
        const std::string fileType = argv[3];
        const std::string pathToLimes = argv[4];
        const std::string pathToProgress = argv[5];
        Files g1 = retrieve_directory_content(dir1, fileType);
        Files g2 = retrieve_directory_content(dir2, fileType);
        if (g1.empty() || g2.empty()) return EXIT_FAILURE;
        
        bool g1Assembled = true;//isAssembled(g1);
        bool g2Assebmled = true;//isAssembled(g2);
        if (g1Assembled && g2Assebmled) {
            std::cout << "Both are assembled" << std::endl;
            run(g1, g2, pathToLimes, pathToProgress);
        }
        else if((g1Assembled && !g2Assebmled) || (!g1Assembled && g2Assebmled)) {
            std::cout << "one is assembled but not the other" << std::endl;
            
            if (g1Assembled && !g2Assebmled)    run(g2.front(), g1, pathToLimes, pathToProgress);
            else    run(g1.front(), g2, pathToLimes, pathToProgress);
            
        } else {
            std::vector<std::string>::size_type i, j, size_g1 = g1.size(), size_g2 = g2.size();
            std::cout << "first genome: " << dir1 << " has " << size_g1 << " number of files" << std::endl;
            std::cout << "second genome: " << dir2 << " has " << size_g2 << " number of files" << std::endl;
            std::cout << "neither is assembled" << std::endl;
            for (i = 0; i < size_g1; ++i) {
                for (j = 0; j < size_g2; ++j) {
                    run(g1.front(), g2.front(), pathToLimes, pathToProgress);
                }
            }
        }
        return EXIT_SUCCESS;
    }
}

#pragma mark -
