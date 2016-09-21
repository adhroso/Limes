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

typedef int INDEX;
typedef std::string WORD;
typedef std::string RawData;
typedef std::vector<WORD> NumberHash;
typedef std::vector<int> Chromosome;
typedef std::map<WORD, INDEX> StringHash;



class Element {
public:
    Element() : id(-1), idx(-1){}
    Element(const int ID, const int index) : id(ID), idx(index){}
    void operator()(const int ID, const int index){
        id = ID;
        idx = index;
    }
    int id, idx;
};
typedef std::vector<Element> Pos;
typedef std::vector<Pos> Vector2D;

class Candidate {
public:
    Candidate():idx_1(-1),idx_2(-1) {}
    Candidate(const int idx1, const int idx2):idx_1(idx1),idx_2(idx2) {}
    int idx_1, idx_2;   //idx_1 is the index from first chromosome, and idx_2 is the index from the second chromosome
};
typedef std::vector<Candidate> Candidates;

#pragma mark - Pre-process
void initializeStringHash();
void initializeNumberHash();
void initializeHash();
//void initializeVec2D(const Chromosome &);

#pragma mark - IO
RawData loadDataWithContentsOFile(const std::string&);
RawData loadDataWithContentsOFile(const std::string &file, std::string &sequenceHeader);
void print_string_stdout(const char *s);
void info(const char *fmt,...);

#endif
