//
//  Serial.cpp
//  Limes
//
//  Created by Andi Dhroso on 10/30/13.
//  Copyright (c) 2013 Andi Dhroso. All rights reserved.
//

#include "Serial.h"

Lime limeObjs_target;
Lime limeObjs_query;

Vector2D vec2D;
scottgs::Timing timer;

#pragma mark - 
/**
 transforms a 'word' +w+ of length +len+ from the 4-letter nucleotide
 alphabet (ACGT) to a dense encoding
 @note skips any letter not in {A,C,G,T}
 @param word - the 'word' to encode
 @param len - the length of the word, if +len+<=8 the leftmost bits
 are guaranteed to be 0x0000 and the value can safely be stored in a 16-bit
 variable. Be aware of endianess during conversion and use ntohl().
 @return returns a minimal, unique representation
 */
u_int32_t encodeWord(const u_char* word, unsigned short len, bool &valid){
	assert(len<=sizeof(unsigned int)*4);
    valid = true;
	unsigned int buffer=0;
	unsigned short i;
	for(i=0;i<len;i++){
		switch(word[i]){
			case 'A':
			case 'a':
				buffer=(buffer<<2)|BA_ENCODED_A;
				break;
			case 'C':
			case 'c':
				buffer=(buffer<<2)|BA_ENCODED_C;
				break;
			case 'G':
			case 'g':
				buffer=(buffer<<2)|BA_ENCODED_G;
				break;
			case 'T':
			case 't':
				buffer=(buffer<<2)|BA_ENCODED_T;
				break;
			default: //skips any other letter!
				valid = false;
                break;
		}
	}
	return buffer;
}

void reverse_complement(RawData &raw_data) {
    std::reverse(raw_data.begin(), raw_data.end());
    std::for_each(raw_data.begin(), raw_data.end(), [](char &c){
        switch (c) {
            case 'A':
            case 'a':
                c='T';
                break;
            case 'C':
            case 'c':
                c='G';
                break;
            case 'G':
            case 'g':
                c='C';
                break;
            case 'T':
            case 't':
                c='A';
                break;
            default:
                break;
        }
    });
}

std::vector<std::string> list_directory_content( const std::string &dir_path, const std::string &ext) {
    std::vector<std::string> fileNames;
    fileNames.reserve(50000);
    
    struct dirent* dp = NULL;
    struct stat st;
    
    DIR* dirp = opendir(dir_path.c_str());
    while (dirp) {
        if ((dp = readdir(dirp)) != NULL) {
            std::string s (dir_path+dp->d_name);
            int status = lstat(s.c_str(), &st);
            if (status != -1) {
                if (S_ISREG(st.st_mode) && (s.substr(s.length()-ext.length()) == ext)) {
                    fileNames.push_back(s);
                    //info("%s\n", s.c_str());
                }
            }
        } else {
            break;
        }
    }
    
    if(dirp)
        closedir(dirp);
    else {
        std::cerr << "error, directory path not valid: \n" << dir_path << std::endl;
    }
    
    return fileNames;
}

Chromosome generateHashWithContentsOfData_approach_1(const RawData &data) {
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    size_t length = data.size();
    Chromosome chromosome;
    chromosome.reserve(length);
    
    bool valid;
    size_t i;
    for (i = 0; i < length; i+=chunk) {
        const std::string w = data.substr(i,WORDSIZE);
        const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
        chromosome.push_back(valid ? hash : -1);
    }
    return chromosome;
}

void generateHashWithContentsOfData_approach_1(const RawData &data, Chromosome &c) {
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    size_t length = data.size();
    
    Chromosome chromosome;
    chromosome.reserve(length);
    
    bool valid;
    size_t i;
    for (i = 0; i < length; i+=chunk) {
        const std::string w = data.substr(i,WORDSIZE);
        const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
        chromosome.push_back(valid ? hash : -1);
    }
    chromosome.swap(c);
}

Chromosome generateHashWithContentsOfData_approach_1_B(const RawData &data) {
    size_t length = data.size()-WORDSIZE+1;
    bool valid;
    
    Chromosome chromosome2;
    chromosome2.resize(length);
    for (size_t i = 0; i < length; ++i) {
        const std::string w = data.substr(i,WORDSIZE);
        const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
        chromosome2[i] = valid ? hash : -1;
    }
    return chromosome2;
/*
    __block Chromosome chromosome((int)length,-1);
    size_t strideCount = 8;
    size_t stride = floor(length / strideCount);
    
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_apply(strideCount, queue, ^(size_t idx){
        size_t i = idx*stride;
        size_t stop = i+stride;
        bool valid;
        do {
            const std::string w = data.substr(i,WORDSIZE);
            const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
            chromosome[i] = valid ? hash : -1;
        } while (++i < stop);
    });
    
    for (size_t i = length - (length % stride); i < length; ++i) {
        const std::string w = data.substr(i,WORDSIZE);
        const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
        chromosome[i] = valid ? hash : -1;
    }
    //chromosome.erase(std::remove(chromosome.begin(), chromosome.end(), -1),chromosome.end());
    //    std::cout << "parallel chromosome size = " << chromosome.size() << std::endl;
    //    std::cout << "parallel chromosome size = " << chromosome2.size() << std::endl;
    //    std::cout << "chromosomes are " << (verify_hash_generation(chromosome2, chromosome)? "equal" : "not equal") << std::endl;
    
    return chromosome;
 */
}

void initializeVec2D_approach_1(const Chromosome &chr, std::string s="") {  //do we need second parameter (s)
    vec2D.clear();
    int size = static_cast<int>(pow(4, WORDSIZE));
    //initialize vec2D - holds location(s) of each possible word
    vec2D.reserve(size);
    for (int i = 0; i < size; ++i) {
        Pos v;
        v.reserve(4000);
        vec2D.push_back(v);
    }
    
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    size = static_cast<int>(chr.size()-chunk);
    for (int i = 0; i < size; ++i) {
        assert(i+chunk < chr.size());
        const int begin = chr[i];
        const int end = chr[i+chunk];
        if (begin > -1 && end > -1) {            // no negative values should be in Element.id
            Pos & pos = vec2D[begin];
            pos.push_back(Element(end, i));    //index (i) refers to the beginning
        }
    }
}

long offset(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2 ) {
    std::pair<std::string::const_iterator, std::string::const_iterator> pair;
    pair = std::mismatch(first1, last1, first2, [](char a, char b){
        switch (a) {
            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
                return a == b;
            default:
                return false;
        }
    });
    return  std::distance(first1, pair.first);
}

long offset(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2 ) {
    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
    pair = std::mismatch(first1, last1, first2, [](char a, char b){
        switch (a) {
            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
                return a == b;
            default:
                return false;
        }
    });
    return  std::distance(first1, pair.first);
}

#pragma mark -
//d1 = query
//d2 = target
void remove_invalid_lime_candidates(Candidates &candidates, const RawData &d1, const RawData &d2) {
    if (candidates.empty()) return;
    
    Candidates::size_type length = candidates.size();
    const RawData::size_type size1 = d1.length(), size2 = d2.length();
    
    std::string::const_iterator d1_begin, d1_end=d1.end(), d2_begin, d2_end = d2.end();
    std::string::const_reverse_iterator d1_rbegin, d1_rend=d1.rend(), d2_rbegin, d2_rend = d2.rend();
    
    for (Candidates::size_type i = 0; i < length; ++i) {
        Candidate &c = candidates[i];
        const int idx_1 = c.idx_1;
        const int idx_2 = c.idx_2;
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (left to right)        //
        //////////////////////////////////////////////////////////
        d1_begin = d1.begin();
        d2_begin = d2.begin();
        std::advance(d1_begin, idx_1);
        std::advance(d2_begin, idx_2);
//        if (print && idx_2 == 1) {
//            
//            out << "Forward direction - forward" << std::endl;
//            out << "i = " << i <<  std::endl;
//            out << "idx 1 = " << idx_1 << std::endl;
//            out << "idx 2 = " << idx_2 << std::endl;
//            
//            out << "query iterator distance from beggining: " << std::distance(d1.begin(),d1_begin) << std::endl;
//            out << "target iterator distance from beggining: " << std::distance(d2.begin(),d2_begin) << std::endl;
//            
//            out << "Query size " << d1.size() << std::endl;
//            out << "Target size " << d2.size() << std::endl;
//        }

        const int queryDistanceToEnd = static_cast<int>(size1-idx_1);
        const int targetDistanceToEnd = static_cast<int>(size2-idx_2);
        //const long right_offset = offset(d1_begin, d1_end, d2_begin);
        const long right_offset = (queryDistanceToEnd < targetDistanceToEnd) ? offset(d1_begin, d1_end, d2_begin) : offset(d2_begin, d2_end, d1_begin);
        if (right_offset < SEQLENGTH/2 - WORDSIZE/2)
            continue;
        //////////////////////////////////////////////////////////
        // scan in the forward direction (right to left)        //
        //////////////////////////////////////////////////////////
        d1_rbegin = d1.rbegin();
        d2_rbegin = d2.rbegin();
        std::advance(d1_rbegin, size1-idx_1);
        std::advance(d2_rbegin, size2-idx_2);
//        if (print && idx_2 == 1) {
//            out << "Reverse direction - reverse" << std::endl;
//            out << "i = " << i << std::endl;
//            out << "idx 1 = " << idx_1 << std::endl;
//            out << "idx 2 = " << idx_2 << std::endl;
//            
//            out << "query iterator distance from beggining: " << std::distance(d1.rbegin(),d1_rbegin) << std::endl;
//            out << "target iterator distance from beggining: " << std::distance(d2.rbegin(),d2_rbegin) << std::endl;
//            
//            out << "Query size " << d1.size() << std::endl;
//            out << "Target size " << d2.size() << std::endl;
//            
//        }

        //const long left_offset = offset(d1_rbegin, d1_rend, d2_rbegin);
        const long left_offset = (idx_1 < idx_2) ? offset(d1_rbegin, d1_rend, d2_rbegin) : offset(d2_rbegin, d2_rend, d1_rbegin);
        
        const long diff = right_offset+left_offset;
        if (diff  >= SEQLENGTH) {
            limeObjs_target.push_back(std::make_pair(idx_2-left_offset,diff));
            limeObjs_query.push_back(std::make_pair(idx_1-left_offset,diff));
        }
        //else
        //    c.idx_1=-1;
        
    }
    //used only for debugging
    //candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [](const Candidate& c) { return c.idx_1 < 0;}), candidates.end());
}

//void remove_invalid_lime_candidates(Candidates &candidates, const RawData &d1, const RawData &d2) {
//    if (candidates.empty()) return;
//    
//    Candidates::size_type length = candidates.size();
//    
//    const RawData::size_type size1 = d1.length(), size2 = d2.length();
//    std::string::const_iterator d1_begin, d1_end=d1.end(), d2_begin, d2_end=d2.end(), start_1, start_2;
//    std::string::const_reverse_iterator d1_rbegin, d1_rend=d1.rend(), d2_rbegin, d2_rend=d2.rend(), start_1r, start_2r;
//    
//    for (Candidates::size_type i = 0; i < length; ++i) {
//        Candidate &c = candidates[i];
//        const int idx_1 = c.idx_1;
//        const int idx_2 = c.idx_2;
//        
//        //////////////////////////////////////////////////////////
//        // scan in the forward direction (left to right)        //
//        //////////////////////////////////////////////////////////
//        d1_begin = d1.begin();
//        std::advance(d1_begin, idx_1);
//        
//        d2_begin = d2.begin();
//        std::advance(d2_begin, idx_2);
//        
//        start_1 = d1_begin;
//        start_2 = d2_begin;
//        while (d1_begin != d1_end && d2_begin != d2_end) {
//            bool valid;
//            switch (*d1_begin) {
//                case 'A':   case 'C':   case 'G':   case 'T':
//                case 'a':   case 'c':   case 'g':   case 't':
//                    valid = true;
//                    break;
//                default:
//                    valid = false;
//                    break;
//            }
//            
//            if (valid && *d1_begin == *d2_begin) {
//                ++d1_begin;
//                ++d2_begin;
//            } else {
//                break;
//            }
//        }
//        const long right_offset = std::distance(start_1, d1_begin);
//        
//        
//        //////////////////////////////////////////////////////////
//        // scan in the forward direction (right to left)        //
//        //////////////////////////////////////////////////////////
//        d1_rbegin = d1.rbegin();
//        std::advance(d1_rbegin, size1-idx_1);
//        
//        d2_rbegin = d2.rbegin();
//        std::advance(d2_rbegin, size2-idx_2);
//        
//        start_1r = d1_rbegin;
//        start_2r = d2_rbegin;
//        while (d1_rbegin != d1_rend && d2_rbegin != d2_rend) {
//            bool valid;
//            switch (*d1_rbegin) {
//                case 'A':   case 'C':   case 'G':   case 'T':
//                case 'a':   case 'c':   case 'g':   case 't':
//                    valid = true;
//                    break;
//                default:
//                    valid = false;
//                    break;
//            }
//            
//            if (valid && *d1_rbegin == *d2_rbegin) {
//                ++d1_rbegin;
//                ++d2_rbegin;
//            } else {
//                break;
//            }
//        }
//        const long left_offset = std::distance(start_1r, d1_rbegin);
//        const int diff = static_cast<int>(right_offset+left_offset);
//        if (diff  >= SEQLENGTH) {
//            limeObjs_target.push_back(std::make_pair(idx_2-left_offset,diff));
//            limeObjs_query.push_back(std::make_pair(idx_1-left_offset,diff));
//        } else
//            c.idx_1=-1;
//        
//    }
//    //used only for debugging
//    candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [](const Candidate& c) { return c.idx_1 < 0;}), candidates.end());
//}

void remove_invalid_lime_candidates_reverse(Candidates &candidates, const RawData &query_data, const RawData &target_data) {
    if (candidates.empty()) return;
    
    Candidates::size_type length = candidates.size();
    const RawData::size_type size1 = query_data.length(), size2 = target_data.length();
    
    std::string::const_iterator d1_begin, d1_end=query_data.end(), d2_begin, d2_end = target_data.end();
    std::string::const_reverse_iterator d1_rbegin, d1_rend=query_data.rend(), d2_rbegin, d2_rend = target_data.rend();
    
    //pick up any left over iterations
    for (Candidates::size_type i = 0; i < length; ++i) {
        Candidate &c = candidates[i];
        const int idx_1 = c.idx_1;
        const int idx_2 = c.idx_2;
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (left to right)        //
        //////////////////////////////////////////////////////////
        d1_begin = query_data.begin();
        d2_begin = target_data.begin();
        std::advance(d1_begin, idx_1);
        std::advance(d2_begin, idx_2);
//        if ( print && idx_2 == 1) {
//            out << "Forward direction - reverse" << std::endl;
//            out << "i = " << i <<  std::endl;
//            out << "idx 1 = " << idx_1 << std::endl;
//            out << "idx 2 = " << idx_2 << std::endl;
//            
//            out << "query iterator distance from beggining: " << std::distance(query_data.begin(),d1_begin) << std::endl;
//            out << "target iterator distance from beggining: " << std::distance(target_data.begin(),d2_begin) << std::endl;
//            
//            out << "Query size " << query_data.size() << std::endl;
//            out << "Target size " << target_data.size() << std::endl;
//        }
        const int queryDistanceToEnd = static_cast<int>(size1-idx_1);
        const int targetDistanceToEnd = static_cast<int>(size2-idx_2);
        //const long right_offset = offset(d1_begin, d1_end, d2_begin);
        const long right_offset = (queryDistanceToEnd < targetDistanceToEnd) ? offset(d1_begin, d1_end, d2_begin) : offset(d2_begin, d2_end, d1_begin);
        if (right_offset < SEQLENGTH/2 - WORDSIZE/2)
            continue;
        
        //////////////////////////////////////////////////////////
        // scan in the reverse direction (right to left)        //
        //////////////////////////////////////////////////////////
        d1_rbegin = query_data.rbegin();
        d2_rbegin = target_data.rbegin();
        std::advance(d1_rbegin, size1-idx_1);
        std::advance(d2_rbegin, size2-idx_2);

//        if (print && idx_2 == 1){
//            out << "Reverse direction - reverse" << std::endl;
//            out << "i = " << i << std::endl;
//            out << "idx 1 = " << idx_1 << std::endl;
//            out << "idx 2 = " << idx_2 << std::endl;
//            
//            out << "query iterator distance from beggining: " << std::distance(query_data.rbegin(),d1_rbegin) << std::endl;
//            out << "target iterator distance from beggining: " << std::distance(target_data.rbegin(),d2_rbegin) << std::endl;
//            
//            out << "Query size " << query_data.size() << std::endl;
//            out << "Target size " << target_data.size() << std::endl;
//
//        }
        //const long left_offset = offset(d1_rbegin, d1_rend, d2_rbegin);
        const long left_offset = (idx_1 < idx_2) ? offset(d1_rbegin, d1_rend, d2_rbegin) : offset(d2_rbegin, d2_rend, d1_rbegin);
        
        const long diff = right_offset+left_offset;
        if (diff  >= SEQLENGTH) {
            const RawData::size_type query_data_size = query_data.size();
            const RawData::size_type query_start = query_data_size - (idx_1+right_offset);
            limeObjs_target.push_back(std::make_pair(idx_2-left_offset, diff));
            limeObjs_query.push_back(std::make_pair(query_start, diff));
         }
        //else {
        //    c.idx_1=-1;
        //}
    }
    //used only for debugging
    //candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [](const Candidate& c) { return c.idx_1 < 0;}), candidates.end());
}


//void remove_invalid_lime_candidates_reverse(Candidates &candidates, const RawData &query_data, const RawData &target_data) {
//    if (candidates.empty()) return;
//    
//    Candidates::size_type length = candidates.size();
//    
//    const RawData::size_type size1 = query_data.length(), size2 = target_data.length();
//    std::string::const_iterator d1_begin, d1_end=query_data.end(), d2_begin, d2_end=target_data.end(), start_1, start_2;
//    std::string::const_reverse_iterator d1_rbegin, d1_rend=query_data.rend(), d2_rbegin, d2_rend=target_data.rend(), start_1r, start_2r;
//
//    //pick up any left over iterations
//    for (Candidates::size_type i = 0; i < length; ++i) {
//        Candidate &c = candidates[i];
//        const int idx_1 = c.idx_1;
//        const int idx_2 = c.idx_2;
//        
//        //////////////////////////////////////////////////////////
//        // scan in the forward direction (left to right)        //
//        //////////////////////////////////////////////////////////
//        d1_begin = query_data.begin();
//        std::advance(d1_begin, idx_1);
//        
//        d2_begin = target_data.begin();
//        std::advance(d2_begin, idx_2);
//        
//        start_1 = d1_begin;
//        start_2=d2_begin;
//        while (d1_begin != d1_end && d2_begin != d2_end) {
//            bool valid;
//            switch (*d1_begin) {
//                case 'A':   case 'a':   case 'C':   case 'c':
//                case 'G':   case 'g':   case 'T':   case 't':
//                    valid = true;
//                    break;
//                default:
//                    valid = false;
//                    break;
//            }
//            
//            if (valid && *d1_begin == *d2_begin) {
//                ++d1_begin;
//                ++d2_begin;
//            } else {
//                break;
//            }
//        }
//        const long right_offset = std::distance(start_1, d1_begin);
//        
//        
//        //////////////////////////////////////////////////////////
//        // scan in the reverse direction (right to left)        //
//        //////////////////////////////////////////////////////////
//        d1_rbegin = query_data.rbegin();
//        std::advance(d1_rbegin, size1-idx_1);
//        
//        d2_rbegin = target_data.rbegin();
//        std::advance(d2_rbegin, size2-idx_2);
//        
//        start_1r = d1_rbegin;
//        start_2r = d2_rbegin;
//        while (d1_rbegin != d1_rend && d2_rbegin != d2_rend) {
//            bool valid;
//            switch (*d1_rbegin) {
//                case 'A':   case 'a':   case 'C':   case 'c':
//                case 'G':   case 'g':   case 'T':   case 't':
//                    valid = true;
//                    break;
//                default:
//                    valid = false;
//                    break;
//            }
//            
//            if (valid && *d1_rbegin == *d2_rbegin) {
//                ++d1_rbegin;
//                ++d2_rbegin;
//            } else {
//                break;
//            }
//        }
//        const long left_offset = std::distance(start_1r, d1_rbegin);
//        const int diff = static_cast<int>(right_offset+left_offset);
//        if (diff  >= SEQLENGTH) {
//            const RawData::size_type query_data_size = query_data.size();
//            const RawData::size_type query_start = query_data_size - (idx_1+right_offset);
//            limeObjs_target.push_back(std::make_pair(idx_2-left_offset, diff));
//            limeObjs_query.push_back(std::make_pair(query_start, diff));
//        } else {
//            c.idx_1=-1;
//        }
//    }
//    //used only for debugging
//    candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [](const Candidate& c) { return c.idx_1 < 0;}), candidates.end());
//}

void generate_lime_candidates(const int i, const Chromosome &query_chr, Candidates &candidates ) {
    //Limes limes;
    const int chr1_begin = query_chr[i];
    if (chr1_begin < 0)
        return;
    
    const Pos & vec2DPos = vec2D[chr1_begin];
    if (vec2DPos.empty()) return;
    
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    const int length = static_cast<int>(vec2DPos.size());
    
    const int chr1_end = query_chr[i+1];
    for (int j = 0; j < length; ++j) {
        const Element & element = vec2DPos[j];
        const int chr2_end = element.id;
        if (chr1_end == chr2_end) { //chr1_end might be negative, but element.id (chr2_end) should never be negative
            candidates.push_back(Candidate(i*chunk, element.idx));
        }
    }
}

void find_candidate_limes(Chromosome &query_chr, Candidates &candidates) {
    const int iterations = static_cast<int>(query_chr.size()-1);
//    const int percent = ceil(iterations*.10);
    
//    std::cout << "\nProgress..." << std::endl;
    for (int i = 0; i < iterations; ++i) {
//        if ((i % percent) == 0)
//            info("Amount complete = %d\n", i);
        generate_lime_candidates(i, query_chr, candidates);
    }
}

void find_limes_in_forward_direction(Chromosome &query_chr, const RawData &query_data, const RawData &target_data, Candidates &candidates) {
    //////////////////////////
    // Find candidate limes //
    //////////////////////////
    //timer.split();
    find_candidate_limes(query_chr, candidates);
    //std::cout << "\nNumber of candidate limes = " << candidates.size() << std::endl;
    //std::cout << "First stage processing time = " << timer.getSplitElapsedTime() << std::endl;
    
    //////////////////////////////////
    // Remove invalid candidates    //
    //////////////////////////////////
    //timer.split();
    remove_invalid_lime_candidates(candidates, query_data, target_data);
    // std::cout << "\nNumber of candidate limes = " << candidates.size() << std::endl;
    //std::cout << "Second stage processing time = " << timer.getSplitElapsedTime() << std::endl;
    //std::cout << "**************************************************************" << std::endl;
}

void find_limes_in_reverse_direction(Chromosome &query_chr, const RawData &query_data, const RawData &target_data, Candidates &candidates) {
    //////////////////////////
    // first stage filter   //
    //////////////////////////
//    timer.split();
    find_candidate_limes(query_chr, candidates);

//    std::cout << "\nNumber of candidate limes = " << candidates.size() << std::endl;
//    std::cout << "First stage processing time = " << timer.getSplitElapsedTime() << std::endl;
    
    //////////////////////////
    // second stage filter  //
    //////////////////////////
//    timer.split();
    remove_invalid_lime_candidates_reverse(candidates, query_data, target_data);
//    std::cout << "\nNumber of candidate limes = " << candidates.size() << std::endl;
//    std::cout << "Second stage processing time = " << timer.getSplitElapsedTime() << std::endl;
//    std::cout << "**************************************************************" << std::endl;
}

//target    = does not change
//query     = we scan it in the forward direction, then in the reverse direction
//a         = generated by target
void find_all_limes_with_pair(const RawData &target_data, const RawData &query_data, const Chromosome& target_chr) {
//    timer.split();
    RawData query_data_cloned (query_data);
    Chromosome query_chr = generateHashWithContentsOfData_approach_1(query_data_cloned);    //skips every chunk size (SEQLENGTH/2 - WORDSIZE/2)
//    std::cout << "Translation time for rat chromosome = " << timer.getSplitElapsedTime() << std::endl;
    
    Candidates candidates;
    candidates.reserve(6000000);
    
    //////////////////////////
    // forward scan         //
    //////////////////////////
    find_limes_in_forward_direction(query_chr, query_data_cloned, target_data, candidates);
    candidates.clear();

    //////////////////////////
    // reverse-complement   //
    //////////////////////////
//    timer.split();
    reverse_complement(query_data_cloned);
//    std::cout << "Chromosome reverse-complement processing time = " << timer.getSplitElapsedTime() << std::endl;
//    timer.split();
    
    //////////////////////////
    // reverse scan         //
    //////////////////////////
    generateHashWithContentsOfData_approach_1(query_data_cloned, query_chr);
//    std::cout << "Chromosome translation processing time = " << timer.getSplitElapsedTime() << std::endl;
    
    find_limes_in_reverse_direction(query_chr, query_data_cloned, target_data, candidates);
}

void find_all_limes_with_pair(const RawData &target_data,  RawData &query_data, const Chromosome& target_chr) {
    //    timer.split();
    Chromosome query_chr = generateHashWithContentsOfData_approach_1(query_data);    //skips every chunk size (SEQLENGTH/2 - WORDSIZE/2)
    //    std::cout << "Translation time for rat chromosome = " << timer.getSplitElapsedTime() << std::endl;
    
    Candidates candidates;
    candidates.reserve(6000000);
    
    //////////////////////////
    // forward scan         //
    //////////////////////////
    find_limes_in_forward_direction(query_chr, query_data, target_data, candidates);
    candidates.clear();
    
    //////////////////////////
    // reverse-complement   //
    //////////////////////////
    //    timer.split();
    reverse_complement(query_data);
    //    std::cout << "Chromosome reverse-complement processing time = " << timer.getSplitElapsedTime() << std::endl;
    //    timer.split();
    
    //////////////////////////
    // reverse scan         //
    //////////////////////////
    generateHashWithContentsOfData_approach_1(query_data, query_chr);
    //    std::cout << "Chromosome translation processing time = " << timer.getSplitElapsedTime() << std::endl;
    
    find_limes_in_reverse_direction(query_chr, query_data, target_data, candidates);
}

bool optimize(const std::vector<std::string> &g1, const std::vector<std::string> &g2) {
    std::vector<std::string>::size_type length = g1.size() < g2.size() ? g1.size() : g2.size();
    int g1Votes = 0, g2Votes = 0;
    for (std::vector<std::string>::size_type i = 0; i < 10 && i < length; ++i) {
        std::vector<std::string>::size_type g1Size = loadDataWithContentsOFile(g1[i]).size();
        std::vector<std::string>::size_type g2Size = loadDataWithContentsOFile(g2[i]).size();
        g1Votes += g1Size < g2Size ? -1 : 1;
    }
    return g1Votes < g2Votes;
}

void generate_limes(const std::vector<std::string> &g1, const std::vector<std::string> &g2, const std::string &pathToLimes, const std::string &pathToProgress) {
    timer.start();
    const bool shouldSwap = false;//optimize(g1, g2);
    std::vector<std::string>genome1 = g1, genome2 = g2;
    if (shouldSwap)
        genome1.swap(genome2);
    
    const int totalComparisons = static_cast<int>(g1.size()*g2.size());
    const std::vector<std::string>::size_type g1Size = genome1.size();
    const std::vector<std::string>::size_type g2Size = genome2.size();
    
    limeObjs_target.reserve(100000);
    limeObjs_query.reserve(1000000);
    Chromosome target_chr, query_chr;
    RawData target_data, query_data;
    
    std::string targetSeqHeader, querySeqHeader;
    
    std::ofstream out (pathToLimes.c_str());
    std::ofstream progress (pathToProgress.c_str());
    //double io_time = 0;
    
    int counter = 1;
    //outter loop = the larger raw data size
    for (std::vector<std::string>::size_type i = 0; i < g1Size; ++i) {
        const std::string targetSeqFile(genome1[i]);
        //timer.split();
        target_data = loadDataWithContentsOFile(targetSeqFile, targetSeqHeader);
        //io_time += timer.getSplitElapsedTime();

        target_chr = generateHashWithContentsOfData_approach_1_B(target_data);          //skips by 1 letter
        initializeVec2D_approach_1(target_chr, target_data);                            //we want to minimize vec2D generation
        
        for (std::vector<std::string>::size_type j = 0; j < g2Size; ++j) {
            const std::string querySeqFile(genome2[j]);
            //timer.split();
            query_data = loadDataWithContentsOFile(querySeqFile, querySeqHeader);
            //io_time += timer.getSplitElapsedTime();
            
            timer.split();
            progress << "Processing..." << counter << "/" << totalComparisons << "\n";
            progress << targetSeqFile << "\n" << querySeqFile << std::endl;
            
            find_all_limes_with_pair(target_data,query_data,target_chr);
            
            progress << "Processing time = " << timer.getSplitElapsedTime() << "\n" << std::endl;
            
            //write limes to file
            if (!limeObjs_target.empty()) {
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
                out << "#1" << querySeqHeader << std::endl;
                out << "#start_1" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_query.begin(), q_it, [&out](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                out << "#2" << targetSeqHeader << std::endl;
                out << "#start_2" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_target.begin(), t_it, [&out](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                limeObjs_target.clear();
                limeObjs_query.clear();
            }
            counter++;
        }
    }
    //progress << "Total io time = " << io_time << std::endl;
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    out.close();
    progress.close();
    timer.stop();
}

void generate_limes(const std::string &file, const std::vector<std::string> &g, const std::string &limes, const std::string &progress_file) {
    if ( !file_exists(file)  ) {
		std::cerr << file << std::endl;
		std::cerr << "error: path to file not found...exiting." << std::endl;
		return ;
	}
    
    timer.start();
    limeObjs_target.reserve(100000);
    limeObjs_query.reserve(100000);
    
    std::vector<std::pair<std::string, std::string> > target_sequences; //first=header, second=sequence
    std::vector<std::pair<std::string, std::string> > query_sequences;  //first=header, second=sequence
    
    //first parameter
    std::ifstream in_query (file.c_str());
    load_next_batch(in_query, query_sequences, 1000000);
    std::vector<std::pair<std::string, std::string> >::size_type number_of_query_seq = query_sequences.size();
    in_query.close();
    
    //second parameter
    std::vector<std::string>::size_type number_of_files = g.size();
    
    RawData query_sequence, target_sequence;
    std::string target_header;
    
    std::ofstream out (limes.c_str());
    std::ofstream progress (progress_file.c_str());
    for (std::vector<std::string>::size_type i = 0; i < number_of_files; ++i) {
        const std::string target_file(g[i]);
        target_sequence = loadDataWithContentsOFile(target_file, target_header);
        const Chromosome target_chr = generateHashWithContentsOfData_approach_1_B(target_sequence);     //skips by 1 letter
        initializeVec2D_approach_1(target_chr, target_sequence);                                        //we want to minimize vec2D generation
        
        timer.split();
        for (std::vector<std::pair<std::string, std::string> >::size_type n = 0; n < number_of_query_seq; ++n) {
            const std::string query_header (query_sequences[n].first);
            query_sequence = query_sequences[n].second;
            
            find_all_limes_with_pair(target_sequence,query_sequence,target_chr);
            
            //write limes to file
            if (!limeObjs_target.empty()) {
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
                out << "(1)" << query_header << std::endl;
                out << "start (1)" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_query.begin(), q_it, [&out,&query_sequence](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                out << "(2)" << target_header << std::endl;
                out << "start (2)" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_target.begin(), t_it, [&out, &target_sequence](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                limeObjs_target.clear();
                limeObjs_query.clear();
            }
        }
        progress << i+1 << "/" << number_of_files << "\t" << number_of_query_seq << "\nProcessing  time = " << timer.getSplitElapsedTime() << std::endl << std::endl;
    }
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    out.close();
    progress.close();
    timer.stop();
}

void generate_limes(const std::string &queryGenome, const std::string &targetG2, const std::string &limes, const std::string &progress_file) {
    if ( !file_exists(queryGenome)) {
		std::cerr << queryGenome << std::endl;
		std::cerr << "error: path to file not found...exiting." << std::endl;
		return ;
	} else if (!file_exists(targetG2)) {
		std::cerr << targetG2 << std::endl;
		std::cerr << "error: path to file not found...exiting." << std::endl;
    }
    
    timer.start();
    limeObjs_target.reserve(100000);
    limeObjs_query.reserve(100000);
    
    std::vector<std::pair<std::string, std::string> > target_sequences; //first=header, second=sequence
    std::vector<std::pair<std::string, std::string> > query_sequences;  //first=header, second=sequence
    
    //query sequences
    std::ifstream in_query (queryGenome.c_str());
    load_next_batch(in_query, query_sequences, 1000000);
    in_query.close();
    
    //target sequences
    std::ifstream in_target (targetG2.c_str());
    load_next_batch(in_target, target_sequences, 1000000);
    in_target.close();
    
    if (target_sequences.size() > query_sequences.size())
        target_sequences.swap(query_sequences);
    
    std::vector<std::pair<std::string, std::string> >::size_type number_of_query_seq = query_sequences.size();
    std::vector<std::pair<std::string, std::string> >::size_type number_of_target_seq = target_sequences.size();
    
    RawData query_sequence, target_sequence;
    std::string target_header;
    
    std::ofstream out (limes.c_str());
    std::ofstream progress (progress_file.c_str());
    
//    scottgs::Timing outputTimer;
//    outputTimer.start();
    double lookuptable_time = 0;
    for (std::vector<std::pair<std::string, std::string> >::size_type i = 0; i < number_of_target_seq; ++i) {
        target_header = target_sequences[i].first;
        const RawData target_sequence = target_sequences[i].second;
        
        timer.split();
        const Chromosome target_chr = generateHashWithContentsOfData_approach_1_B(target_sequence);     //skips by 1 letter
        initializeVec2D_approach_1(target_chr, target_sequence);                                        //we want to minimize vec2D generation
        const double tmp = timer.getSplitElapsedTime();
        lookuptable_time += tmp;
        std::cout << i+1 << "/" << number_of_target_seq << "\t" << number_of_query_seq << "\nLookup table time = " << tmp << std::endl << std::endl;
        
        timer.split();
        for (std::vector<std::pair<std::string, std::string> >::size_type n = 0; n < number_of_query_seq; ++n) {
            const std::string query_header (query_sequences[n].first);
            query_sequence = query_sequences[n].second;
            
            find_all_limes_with_pair(target_sequence,query_sequence,target_chr);
            
            //write limes to file
            if (!limeObjs_target.empty()) {
//                outputTimer.split();
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
                out << "(1)" << query_header << std::endl;
                out << "start (1)" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_query.begin(), q_it, [&out,&query_sequence](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                out << "(2)" << target_header << std::endl;
                out << "start (2)" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_target.begin(), t_it, [&out, &target_sequence](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                limeObjs_target.clear();
                limeObjs_query.clear();
            }
        }
        progress << i+1 << "/" << number_of_target_seq << "\t" << number_of_query_seq << "\nProcessing  time = " << timer.getSplitElapsedTime() << std::endl << std::endl;
    }
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    out.close();
    progress.close();
    timer.stop();
}

void find_all_limes_serial_approach(const std::string &path_to_genome_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output, const std::string &progressFile) {
    timer.start();
    
    std::vector<std::string> genome_a, genome_b;
    
    genome_a = list_directory_content(path_to_genome_a, ext);
    genome_b = list_directory_content(path_to_genome_b, ext);
    const std::vector<std::string>::size_type genome_a_size = genome_a.size();
    const std::vector<std::string>::size_type genome_b_size = genome_b.size();
    
    limeObjs_target.reserve(100000);
    limeObjs_query.reserve(1000000);
    Chromosome target_chr, query_chr;
    RawData target_data, query_data;
    
    std::string targetSeqHeader, querySeqHeader;
    std::ofstream out (output.c_str());
    std::ofstream progress (progressFile.c_str());
    //outter loop = the smaller raw data size
    //outter loop = the larger raw data size
    //for now, enforce this via the order of the input (arguments)
    for (std::vector<std::string>::size_type i = 0; i < genome_a_size; ++i) {
        const std::string file1(genome_a[i]);
        //const std::string file1(path_to_genome_a);
        
//      timer.split();
        target_data = loadDataWithContentsOFile(file1, targetSeqHeader);
//      target_data = loadDataWithContentsOFile(file1);
//        std::cout << "\nLoading time for human = " << timer.getSplitElapsedTime() << std::endl;
        
//        timer.split();
        target_chr = generateHashWithContentsOfData_approach_1_B(target_data);          //skips by 1 letter
//        std::cout << "Translation time for human chromosome = " << timer.getSplitElapsedTime() << std::endl;
        
//        timer.split();
        initializeVec2D_approach_1(target_chr, target_data);                            //we want to minimize vec2D generation
//        std::cout << "Vec2D initialization time " << timer.getSplitElapsedTime() << std::endl<<std::endl;
        
        
        for (std::vector<std::string>::size_type j = 0; j < genome_b_size; ++j) {
            const std::string file2(genome_b[j]);
            //const std::string file2(path_to_genome_b);
            timer.split();
            query_data = loadDataWithContentsOFile(file2, querySeqHeader);
            //query_data = loadDataWithContentsOFile(file2);
//            std::cout << "\nLoading time for rat = " << timer.getSplitElapsedTime() << std::endl;
            
            timer.split();
            progress << "Processing...\n" << file1.c_str() << "\n" << file2.c_str() << std::endl;
            find_all_limes_with_pair(target_data,query_data,target_chr);
            progress << "Processing time = " << timer.getSplitElapsedTime() << "\n" << std::endl;
            
            //write limes to file
            if (!limeObjs_target.empty()) {
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
                out << "(1)" << querySeqHeader << std::endl;
                out << "start (1)" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_query.begin(), q_it, [&out](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                out << "(2)" << targetSeqHeader << std::endl;
                out << "start (2)" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_target.begin(), t_it, [&out](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                limeObjs_target.clear();
                limeObjs_query.clear();
            }
        }
    }
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    out.close();
    progress.close();
    timer.stop();
}

void load_next_batch(std::ifstream &in, std::vector<std::pair<std::string, std::string> > &sequences, const int batch_size) {
    int counter = 0;
    std::vector<std::pair<std::string, std::string> > data;
    data.reserve(batch_size);

    std::string line;
    while (counter < batch_size && in.good()) {
        std::getline(in, line);
        const std::string header(line);
        std::string seq ("");

        while (in.peek() != '>' && in.good()) {
            std::getline(in, line);
            seq.append(line);
        }
        data.push_back(std::make_pair(header, seq));
        counter++;
    }
    data.swap(sequences);
}

bool file_exists(const std::string &path_to_file) {
    struct stat st;
    
    int status = lstat(path_to_file.c_str(), &st);
    if (status != -1) {
        return S_ISREG(st.st_mode);
    }
    return false;
}

bool directory_file_exists(const std::string &path_to_file) {
    struct stat st;
    
    int status = lstat(path_to_file.c_str(), &st);
    if (status != -1) {
        return S_ISDIR(st.st_mode);
    }
    return false;
}


//Larger genome should be first parameter and smaller genome should be the second parameter
void find_all_limes_between_unassembled_genomes_cluster(const std::string &path_to_file, const std::string &path_to_genome, const std::string &ext, const std::string &output, const std::string &progress_file) {
    if ( !file_exists(path_to_file)  ) {
		std::cerr << path_to_file << std::endl;
		std::cerr << "error: path to file not found...exiting." << std::endl;
		return ;
	}

	if ( !directory_file_exists(path_to_genome)) {
		std::cerr << "error: path to directory not found...exiting." << std::endl;
		return ;
	}

    timer.start();
    limeObjs_target.reserve(1000);
    limeObjs_query.reserve(1000);
    
    std::vector<std::pair<std::string, std::string> > target_sequences; //first=header, second=sequence
    std::vector<std::pair<std::string, std::string> > query_sequences;  //first=header, second=sequence
    
    //////////////////////////////
    //          input           //
    //////////////////////////////
    //first parameter
    std::ifstream in_query (path_to_file.c_str());
    load_next_batch(in_query, query_sequences, 26000);
    std::vector<std::pair<std::string, std::string> >::size_type number_of_query_seq = query_sequences.size();
    in_query.close();
    
    //second parameter
    std::vector<std::string> list_of_files;
    list_of_files = list_directory_content(path_to_genome, ext);
    std::vector<std::string>::size_type number_of_files = list_of_files.size();
    
    //////////////////////////////
    //          output          //
    //////////////////////////////
    //output
    std::ofstream progress (progress_file.c_str());
    std::ofstream out (output.c_str());
   //std::ifstream in_target;
    for (std::vector<std::string>::size_type i = 0; i < number_of_files; ++i) {
        const std::string target_file(list_of_files[i]);
        std::ifstream in_target (target_file.c_str());
        load_next_batch(in_target, target_sequences, 26000);
        in_target.close();
       
        const std::vector<std::pair<std::string, std::string> >::size_type number_of_seq = target_sequences.size();
        for (std::vector<std::pair<std::string, std::string> >::size_type j = 0; j < number_of_seq; ++j) {
            const std::string target_header (target_sequences[j].first);
            const RawData target_sequence = target_sequences[j].second;
            const Chromosome target_chr = generateHashWithContentsOfData_approach_1_B(target_sequence);          //skips by 1 letter
            initializeVec2D_approach_1(target_chr, target_sequence);                                             //we want to minimize vec2D generation
                
            timer.split();
            for (std::vector<std::pair<std::string, std::string> >::size_type n = 0; n < number_of_query_seq; ++n) {
                const std::string query_header (query_sequences[n].first);
                RawData query_sequence = query_sequences[n].second;
                
                //find all limes
                find_all_limes_with_pair(target_sequence,query_sequence,target_chr);
                
                //write limes to file
                if (!limeObjs_target.empty()) {
                    assert(limeObjs_target.size() == limeObjs_query.size());
                    std::sort(limeObjs_query.begin(), limeObjs_query.end());
                    std::sort(limeObjs_target.begin(), limeObjs_target.end());
                    const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                    const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                    
                    out << "(1)" << query_header << std::endl;
                    out << "start (1)" << "\t" << "length" << std::endl;
                    std::for_each(limeObjs_query.begin(), q_it, [&out,&query_sequence](const std::pair<std::string::size_type, int> &l){
                        out << l.first << "\t" << l.second << std::endl;
                    });
                    out << "(2)" << target_header << std::endl;
                    out << "start (2)" << "\t" << "length" << std::endl;
                    std::for_each(limeObjs_target.begin(), t_it, [&out, &target_sequence](const std::pair<std::string::size_type, int> &l){
                        out << l.first << "\t" << l.second << std::endl;
                    });
                    limeObjs_target.clear();
                    limeObjs_query.clear();
                }
            }
            progress << i+1 << "/" << number_of_files << "\t" << number_of_query_seq << "\nProcessing  time = " << timer.getSplitElapsedTime() << std::endl << std::endl;
        }
    }
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    out.close();
    progress.close();
    timer.stop();
}


//Larger genome should be first parameter and smaller genome should be the second parameter
void find_all_limes_between_unassembled_genomes_v2(const std::string &path_to_genome_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output, const std::string &progress_) {
    if ( ! (directory_file_exists(path_to_genome_a) && directory_file_exists(path_to_genome_b)) ) return ;
    
    timer.start();
    std::vector<std::string> genome_a_files, genome_b_files;
    genome_a_files = list_directory_content(path_to_genome_a, ext);
    genome_b_files = list_directory_content(path_to_genome_b, ext);
    
    limeObjs_target.reserve(100000);
    limeObjs_query.reserve(100000);
    Chromosome target_chr, query_chr;
    
    std::ofstream progress (progress_.c_str());
    std::ofstream out (output.c_str());
    
    std::vector<std::pair<std::string, std::string> > target_sequences;
    std::vector<std::pair<std::string, std::string> > query_sequences;
    
    std::ifstream in_target;
    for (std::vector<std::string>::size_type i = 0; i < genome_a_files.size(); ++i) {
        const std::string target_file(genome_a_files[i]);
        std::ifstream in_target (target_file.c_str());
        while (in_target.good()) {
            load_next_batch(in_target, target_sequences, 10);
            for (std::vector<std::pair<std::string, std::string> >::size_type j = 0; j < target_sequences.size(); ++j) {
                const std::string target_header (target_sequences[j].first);
                RawData target_sequence = target_sequences[j].second;
                
                target_chr = generateHashWithContentsOfData_approach_1_B(target_sequence);          //skips by 1 letter
                initializeVec2D_approach_1(target_chr, target_sequence);                            //we want to minimize vec2D generation
                
                timer.split();
                for (std::vector<std::string>::size_type k = 0; k < genome_b_files.size(); ++k) {
                    const std::string query_file(genome_b_files[i]);
                    std::ifstream in_query (query_file.c_str());
                    std::vector<std::pair<std::string, std::string> >::size_type n;
                    while (in_query.good()) {
                        load_next_batch(in_query, query_sequences, 360000);
                        for (n = 0; n < query_sequences.size(); ++n) {
                            const std::string query_header (query_sequences[n].first);
                            RawData query_sequence = query_sequences[n].second;
                            
                            find_all_limes_with_pair(target_sequence,query_sequence,target_chr);
                            
                            //write limes to file
                            if (!limeObjs_target.empty()) {
                                assert(limeObjs_target.size() == limeObjs_query.size());
                                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                                const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                                const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                                
                                out << "(1)" << query_header << std::endl;
                                out << "start (1)" << "\t" << "length" << std::endl;
                                std::for_each(limeObjs_query.begin(), q_it, [&out,&query_sequence](const std::pair<std::string::size_type, int> &l){
                                    out << l.first << "\t" << l.second << std::endl;
                                });
                                out << "(2)" << target_header << std::endl;
                                out << "start (2)" << "\t" << "length" << std::endl;
                                std::for_each(limeObjs_target.begin(), t_it, [&out, &target_sequence](const std::pair<std::string::size_type, int> &l){
                                    out << l.first << "\t" << l.second << std::endl;
                                });
                                limeObjs_target.clear();
                                limeObjs_query.clear();
                            }
                        }
                    }
                    in_query.close();
                    progress << "Target files processed: " << i+1 <<  "\nQuery files processed: " << k+1 << "\nProcessing  time = " << timer.getSplitElapsedTime() << std::endl << std::endl;
                }
            }
        }
        in_target.close();
    }
    
    out.close();
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    timer.stop();
}

void find_all_limes_with_sequence(const std::string &path_to_sequences_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output) {
    if ( ! (file_exists(path_to_sequences_a) && directory_file_exists(path_to_genome_b)) ) return ;
    
    timer.start();
    std::vector<std::string> genome_b;
    genome_b = list_directory_content(path_to_genome_b, ext);
    
    //query
    std::ifstream in (path_to_sequences_a.c_str());
    std::vector<std::pair<std::string, std::string> > data;
    std::vector<std::pair<std::string, std::string> > sequences_a;
    
    //target
    const std::vector<std::string>::size_type genome_b_size = genome_b.size();
    
    limeObjs_query.reserve(10000);
    limeObjs_target.reserve(10000);
    Chromosome target_chr, query_chr;
    RawData target_data, query_data;
    
    std::ofstream out (output.c_str(), std::ios::app);
    std::ofstream progress ("/Users/andi/Documents/Projects/Limes/Limes_Serial/progress.txt", std::ios::app);
    for (std::vector<std::string>::size_type i = 0; i < genome_b_size; ++i) {
        const std::string file1(genome_b[i]);
        target_data = loadDataWithContentsOFile(file1);
        target_chr = generateHashWithContentsOfData_approach_1_B(target_data);
        initializeVec2D_approach_1(target_chr, target_data);
        
        while (!in.eof()) {
            load_next_batch(in, sequences_a, 1000);
            for (std::vector<std::pair<std::string, std::string> >::size_type j = 0; j < sequences_a.size(); ++j) {
                const std::string header (sequences_a[j].first);
                query_data = sequences_a[j].second;
                
                //timer.split();
                //info("Processing...\n%s\n%s\n", file1.c_str(), header.c_str());
                find_all_limes_with_pair(target_data,query_data,target_chr);
                //info("Processing time = %f\n\n", timer.getSplitElapsedTime());
                
                //write limes to file
                if (!limeObjs_target.empty()) {
                    assert(limeObjs_target.size() == limeObjs_query.size());
                    std::sort(limeObjs_query.begin(), limeObjs_query.end());
                    std::sort(limeObjs_target.begin(), limeObjs_target.end());
                    const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                    const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                    
                    out << "(1)" << header << std::endl;
                    out << "start (1)" << "\t" << "length" << std::endl;
                    std::for_each(limeObjs_query.begin(), q_it, [&out](const std::pair<std::string::size_type, int> &l){
                        out << l.first << "\t" << l.second << std::endl;
                    });
                    out << "(2)" << file1 << std::endl;
                    out << "start (2)" << "\t" << "length" << std::endl;
                    std::for_each(limeObjs_target.begin(), t_it, [&out](const std::pair<std::string::size_type, int> &l){
                        out << l.first << "\t" << l.second << std::endl;
                    });
                    limeObjs_target.clear();
                    limeObjs_query.clear();
                }
            }
        }
        progress << i << " out of " << genome_b_size << "\nProcessing  time = " << timer.getTotalElapsedTime() << std::endl << std::endl;
        in.clear();
        in.seekg (0, std::ios::beg);
    }
    progress.close();
    out.close();
    in.close();
    
    std::cout << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    timer.stop();
}

//void find_all_limes_between_unassembled_genomes(const std::string &path_to_genome_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output, const std::string &progress) {
//    timer.start();
//    
//    std::vector<std::string> genome_a, genome_b;
//    
//    genome_a = list_directory_content(path_to_genome_a, ext);
//    genome_b = list_directory_content(path_to_genome_b, ext);
//    const std::vector<std::string>::size_type genome_a_size = genome_a.size();
//    const std::vector<std::string>::size_type genome_b_size = genome_b.size();
//    
//    std::vector<std::pair<std::string, std::string> > target_sequences;
//    std::vector<std::pair<std::string, std::string> > query_sequences;
//    
//    limeObjs_target.reserve(100000);
//    limeObjs_query.reserve(1000000);
//    Chromosome target_chr, query_chr;
//    RawData target_data, query_data;
//    
//    std::ofstream out (output.c_str(), std::ios::app);
//    
//    //outter loop = the smaller raw data size
//    //outter loop = the larger raw data size
//    //for now, enforce this via the order of the input (arguments)
//    for (int i = 0; i < genome_a_size; ++i) {
//        const std::string file1(genome_a[i]);
//        //const std::string file1(path_to_genome_a);
//        
//        //        timer.split();
//        target_data = loadDataWithContentsOFile(file1);
//        //        std::cout << "\nLoading time for human = " << timer.getSplitElapsedTime() << std::endl;
//        
//        //        timer.split();
//        target_chr = generateHashWithContentsOfData_approach_1_B(target_data);          //skips by 1 letter
//        //        std::cout << "Translation time for human chromosome = " << timer.getSplitElapsedTime() << std::endl;
//        
//        //        timer.split();
//        initializeVec2D_approach_1(target_chr, target_data);                            //we want to minimize vec2D generation
//        //        std::cout << "Vec2D initialization time " << timer.getSplitElapsedTime() << std::endl<<std::endl;
//        
//        
//        for (int j = 0; j < genome_b_size; ++j) {
//            const std::string file2(genome_b[j]);
//            //const std::string file2(path_to_genome_b);
//            timer.split();
//            query_data = loadDataWithContentsOFile(file2);
//            //            std::cout << "\nLoading time for rat = " << timer.getSplitElapsedTime() << std::endl;
//            
//            timer.split();
//            info("Processing...\n%s\n%s\n", file1.c_str(), file2.c_str());
//            find_all_limes_with_pair(target_data,query_data,target_chr);
//            info("Processing time = %f\n\n", timer.getSplitElapsedTime());
//            
//            //write limes to file
//            if (!limeObjs_target.empty()) {
//                assert(limeObjs_target.size() == limeObjs_query.size());
//                std::sort(limeObjs_query.begin(), limeObjs_query.end());
//                std::sort(limeObjs_target.begin(), limeObjs_target.end());
//                const Lime::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
//                const Lime::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
//                
//                out << "(1)" << file2 << std::endl;
//                out << "start (1)" << "\t" << "length" << std::endl;
//                std::for_each(limeObjs_query.begin(), q_it, [&out](const std::pair<std::string::size_type, int> &l){
//                    out << l.first << "\t" << l.second << std::endl;
//                });
//                out << "(2)" << file1 << std::endl;
//                out << "start (2)" << "\t" << "length" << std::endl;
//                std::for_each(limeObjs_target.begin(), t_it, [&out](const std::pair<std::string::size_type, int> &l){
//                    out << l.first << "\t" << l.second << std::endl;
//                });
//                limeObjs_target.clear();
//                limeObjs_query.clear();
//            }
//        }
//    }
//    out.close();
//    std::cout << "Total time = " << timer.getTotalElapsedTime() << std::endl;
//    timer.stop();
//}




#pragma mark -


