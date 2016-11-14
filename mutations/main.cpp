//
//  main.cpp
//  Limes with single mutation
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#include "Serial.h"


int main(int argc, const char * argv[]) {
//    bool valid = false;
//    std::string w = "GGGAGGGAGGGAGG";
//    std::cout << encodeWord((const u_char*)w.c_str(), WORDSIZE, valid) << std::endl;
//    return 0;
    
//    std::string s1 ("CTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCT");
//    std::string s2 ("TTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCTTCTCCT");
//    
//    std::string::iterator s1_begin = s1.begin();
//    std::advance(s1_begin, 0);
//    
//    std::string::iterator s2_begin = s2.begin();
//    std::advance(s2_begin, 0);
//    
//    int mut_pos = 0;
//    std::pair<std::string::iterator, std::string::iterator> p = std::mismatch(s1_begin, s1.end(), s2_begin);
//    std::cout << "default: " << std::distance(s1_begin, p.first) << ", right offset: " << offset(s1_begin, s1.end(), s2_begin, mut_pos) << ", right_mut: " << mut_pos << std::endl;
//
//    std::string::reverse_iterator s1_rbegin = s1.rbegin();
//    std::advance(s1_rbegin, s1.size());
//    
//    std::string::reverse_iterator s2_rbegin = s2.rbegin();
//    std::advance(s2_rbegin, s2.size());
//    
//    mut_pos=0;
//    std::pair<std::string::reverse_iterator, std::string::reverse_iterator> p2 = std::mismatch(s1_rbegin, s1.rend(), s2_rbegin);
//    std::cout << "default: " << std::distance(s1_rbegin, p2.first) << ", left_offset: " << offset(s1_rbegin, s1.rend(), s2_rbegin, mut_pos) << ", left_mut: " << mut_pos << std::endl;
//
//    
//    
//    
//    
//    
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
