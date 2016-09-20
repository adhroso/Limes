//
//  main.cpp
//  Limes2.0
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#include "Serial.h"


int main(int argc, const char * argv[]) {
    int rc;
    if (argc < 5)  {
        std::cerr << "usage: ./Limes <genome 1 dir> <genome 2 dir> <file_ext> [fa]> <limes file> <progress file>" << std::endl;
        rc = EXIT_FAILURE;
    } else {
        const std::string dir1 = argv[1];
        const std::string dir2 = argv[2];
        const std::string fileType = argv[3];
        const std::string pathToLimes = argv[4];
        const std::string pathToProgress = argv[5];
        Files g1 = retrieve_directory_content(dir1, fileType);
        Files g2 = retrieve_directory_content(dir2, fileType);
        if (g1.empty() || g2.empty()) return EXIT_FAILURE;
        
        bool g1Assembled = isAssembled(g1);
        bool g2Assebmled = isAssembled(g2);
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
        rc = EXIT_SUCCESS;
    }
    //extract_unique_proteins();
    //find_all_limes_serial_approach(argv[1], argv[2], argv[3], argv[4]);
    //find_all_limes_serial_approach(argv[1], argv[2], argv[3], argv[4], argv[5]);
    //find_all_limes_with_sequence(argv[1], argv[2], argv[3], argv[4]);
    //find_all_limes_between_unassembled_genomes_cluster(argv[1], argv[2], argv[3], argv[4], argv[5]);
    return rc;
}

#pragma mark -
