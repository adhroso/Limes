//
//  main.cpp
//  Limes2.0
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#include "Serial.h"


int main(int argc, const char * argv[]) {
    if (argc < 8)  {
        std::cerr << "usage: ./Limes [unassembled: 0, assembled: 1] <genome 1 dir> [unassembled: 0, assembled: 1] <genome 2 dir> <file_ext> [fa]> <limes file> <progress file>" << std::endl;
        return EXIT_FAILURE;
    } else {
        const bool g1Assembled = atoi(argv[1]);
        const std::string dir1 = argv[2];
        
        const bool g2Assebmled = atoi(argv[3]);
        const std::string dir2 = argv[4];
        
        const std::string fileType = argv[5];
        const std::string pathToLimes = argv[6];
        const std::string pathToProgress = argv[7];
        
        Files g1 = retrieve_directory_content(dir1, fileType);
        Files g2 = retrieve_directory_content(dir2, fileType);
        if (g1.empty() || g2.empty()) return EXIT_FAILURE;
        
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
            std::cout << "first genome: " << dir1 << " has " << size_g1 << " files" << std::endl;
            std::cout << "second genome: " << dir2 << " has " << size_g2 << " files" << std::endl;
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
