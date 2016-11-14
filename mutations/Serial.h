//
//  Serial.h
//  Limes
//
//  Created by Andi Dhroso on 10/30/13.
//  Copyright (c) 2013 Andi Dhroso. All rights reserved.
//

#ifndef __Limes__Serial__
#define __Limes__Serial__

#include <fstream>
#include <assert.h>
#include <stdarg.h>
#include <cstdarg>
#include <bitset>

#include "Types.h"
#include "IO.h"
#include "Util.h"
#include "Timing.hpp"

typedef std::pair<std::string::size_type, int> Lime;
typedef std::vector<Lime> Limes;

extern Limes limeObjs_target;
extern Limes limeObjs_query;
extern LookupTable vec2D;
extern scottgs::Timing timer;

void remove_invalid_lime_candidates(Candidates &, const Sequence &, const Sequence &);
void remove_invalid_lime_candidates_reverse(Candidates &, const Query &, const Target &);
void find_lime_candidates(Chromosome &query_chr, Candidates &candidates);

void find_limes_in_forward_direction(QueryChr &, const Query &, const Target &, Candidates &);
void find_limes_in_reverse_direction(QueryChr &, const Query &, const Target &, Candidates &candidates);
void find_limes(const Target &, const Query &, const QueryChr&);

long offset(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2, int &mut_pos);
long offset(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, int &mut_pos);

void run(const Files &, const Files &, const Output &, const Progress &);
void run(const Sequence &, const Files &g, const Output &, const Progress &);
void run(const File &, const File &, const Output &, const Progress &);
u_int32_t encodeWord(const u_char* word, unsigned short len, bool &valid);
#endif /* defined(__Limes__Serial__) */
