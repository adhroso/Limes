//
//  Util.cpp
//  Limes2.0
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#include "Util.h"

bool isAssembled(const File &file) {
    return !(file.find("dna.nonchromosomal") != file.npos);
}

bool isAssembled(Files &g) {
    Files::size_type i;
    for (i = 0; i < g.size(); ++i) {
        const bool assembled = isAssembled(g[i]);
        if (!assembled)
            return false;
    }
    return true;
}

bool file_exists(const Path &path_to_file) {
    struct stat st;
    
    int status = lstat(path_to_file.c_str(), &st);
    if (status != -1) {
        return S_ISREG(st.st_mode);
    }
    return false;
}

bool directory_file_exists(const Path &path_to_file) {
    struct stat st;
    
    int status = lstat(path_to_file.c_str(), &st);
    if (status != -1) {
        return S_ISDIR(st.st_mode);
    }
    return false;
}

