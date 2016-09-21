//
//  IO.hpp
//  Limes2.0
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#ifndef IO_hpp
#define IO_hpp

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "Types.h"


Files retrieve_directory_content( const Path &, const Extension &);

Sequence loadDataWithContentsOFile(const File&);
Sequence loadDataWithContentsOFile(const File &file, Header &header);
void load_next_batch(std::ifstream &, std::vector<std::pair<Header, Sequence> > &, const int batch_size=1000);

void print_string_stdout(const char *s);
void info(const char *fmt,...);

#endif /* IO_hpp */
