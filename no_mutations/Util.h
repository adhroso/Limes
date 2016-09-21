//
//  Util.hpp
//  Limes2.0
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#ifndef Util_hpp
#define Util_hpp

#include <iostream>
#include <fstream>
#include <stdarg.h>
#include <cstdarg>
#include <sys/stat.h>

#include "Types.h"


bool file_exists(const Path &);
bool directory_file_exists(const Path &);
    
bool isAssembled(Files &g);
bool isAssembled(const File &file);

void reverse_complement(Sequence &);


#endif /* Util_hpp */
