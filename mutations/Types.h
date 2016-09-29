//
//  main.h
//  Limes
//
//  Created by Andi Dhroso on 9/20/13.
//  Copyright (c) 2013 Andi Dhroso. All rights reserved.
//

#ifndef Limes_main_h
#define Limes_main_h

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <stdarg.h>
#include <cstdarg>

#pragma mark - Definitions
#define WORDSIZE 12
#define SEQLENGTH 100

#define BA_ENCODED_A 0x00
#define BA_ENCODED_C 0x01
#define BA_ENCODED_G 0x02
#define BA_ENCODED_T 0x03

typedef int Hash;
typedef int Location;
typedef int FirstMutation;
typedef int SecondMutation;
typedef std::pair<FirstMutation, SecondMutation> MutationPos;
typedef std::vector<Hash> Chromosome;
typedef std::vector<Hash> QueryChr;
typedef std::vector<Hash> TargetChr;

typedef std::string Sequence;
typedef std::string Header;
typedef std::string Query;
typedef std::string Target;

typedef std::string File;
typedef std::string Output;
typedef std::string Progress;
typedef std::string Path;
typedef std::string Extension;
typedef std::vector<File> Files;


class Element {
public:
    Element() : id(-1), idx(-1){}
    Element(const Hash ID, const std::size_t index) : id(ID), idx(index){}
    void operator()(const Hash ID, const std::size_t index){
        id = ID;
        idx = index;
    }
    int id;
    std::size_t idx;
};
typedef std::vector<Element> Pos;
typedef Pos Elements;
typedef std::vector<Pos> LookupTable;

class Candidate {
public:
    Candidate():idx_1(-1),idx_2(-1) {}
    Candidate(const std::size_t idx1, const std::size_t idx2):idx_1(idx1),idx_2(idx2) {}
    
    //idx_1 is the index from first chromosome, and idx_2 is the index from the second chromosome
    std::size_t idx_1, idx_2;
};
typedef std::vector<Candidate> Candidates;
#pragma mark - 

#endif


