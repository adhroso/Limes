//
//  Serial.cpp
//  Limes
//
//  Created by Andi Dhroso on 10/30/13.
//  Copyright (c) 2013 Andi Dhroso. All rights reserved.
//

#include "Serial.h"

Limes limeObjs_target;
Limes limeObjs_query;

LookupTable vec2D;
scottgs::Timing timer;

#pragma mark - Utility functions
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
	//assert(len<=sizeof(unsigned int)*4);
    valid = true;
	unsigned int buffer=0;
	unsigned short i;
	for(i=0;i<len;i++){
		switch(word[i]){
			case 'A':   case 'a':
				buffer=(buffer<<2)|BA_ENCODED_A;
				break;
                
			case 'C':   case 'c':
				buffer=(buffer<<2)|BA_ENCODED_C;
				break;
                
			case 'G':   case 'g':
				buffer=(buffer<<2)|BA_ENCODED_G;
				break;
                
			case 'T':   case 't':
				buffer=(buffer<<2)|BA_ENCODED_T;
				break;
                
			default: //skips any other letter!
				valid = false;
                break;
		}
	}
	return buffer;
}

void reverse_complement(Sequence &raw_data) {
    std::reverse(raw_data.begin(), raw_data.end());
    std::for_each(raw_data.begin(), raw_data.end(), [](char &c){
        switch (c) {
            case 'A':   case 'a':
                c='T';
                break;
                
            case 'C':   case 'c':
                c='G';
                break;
                
            case 'G':   case 'g':
                c='C';
                break;
                
            case 'T':   case 't':
                c='A';
                break;
                
            default:
                break;
        }
    });
}

/**
 Linearized vector generation
 */
Chromosome generateVectorIndices(const Sequence &data) {
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    std::size_t length = data.size();
    Chromosome chromosome;
    chromosome.reserve(length);
    
    bool valid;
    for (size_t i = 0; i < length; i+=chunk) {
        const std::string w = data.substr(i,WORDSIZE);
        const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
        chromosome.push_back(valid ? hash : -1);
    }
    return chromosome;
}

void generateVectorIndecies(const Sequence &data, Chromosome &c) {
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    size_t length = data.size();
    
    Chromosome chromosome;
    chromosome.reserve(length);
    
    bool valid;
    for (size_t i = 0; i < length; i+=chunk) {
        const std::string w = data.substr(i,WORDSIZE);
        const u_int32_t hash = encodeWord((const u_char*)w.c_str(), WORDSIZE, valid);
        chromosome.push_back(valid ? hash : -1);
    }
    chromosome.swap(c);
}

/**
    Lookup table generation
 */
Chromosome generateLookupTableIndices(const Sequence &data) {
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
}

void initializeLookupTable(const TargetChr &chr, std::string s="") {  //do we need second parameter (s)
    vec2D.clear();
    
    //////////////////////////////
    //Size of lookup table
    //////////////////////////////
    int size = static_cast<int>(pow(4, WORDSIZE));
    
    //initialize lookuptable - holds location(s) of each possible word
    vec2D.reserve(size);
    for (int i = 0; i < size; ++i) {
        Pos v;
        v.reserve(4000);
        vec2D.push_back(v);
    }
    
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    const size_t length = chr.size()-chunk;
    for (std::size_t i = 0; i < length; ++i) {
        assert(i+chunk < chr.size());
        const Hash begin = chr[i];
        const Hash end = chr[i+chunk];
    
        if (begin > -1 && end > -1) {            // no negative values should be in Element.id
            Pos & pos = vec2D[begin];
            pos.push_back(Element(end, i));      //index (i) refers to the beginning
            
            Pos & pos2 = vec2D[end];
            pos2.push_back(Element(begin * -1, i)); //multiply by -1 to simply flag elements to detect elements with mutations in the prefix
        }
    }
}

//long offset(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2, int &mut_pos) {
//    std::pair<std::string::const_iterator, std::string::const_iterator> pair;
//    int count = 0;
//    int pos = 0;
//    int tmp = 0;
//    pair = std::mismatch(first1, last1, first2, [&count, &pos, &tmp](char a, char b){
//        switch (a) {
//            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
//                if (a!=b) {
//                    count++;
//                    
//                    if (count<2) tmp=pos;
//                }
//                
//                pos++;
//                return count<2 ;
//            default:
//                tmp=pos;
//                return false;
//        }
//    });
//    mut_pos=tmp;
//    //std::cout << "pos: " << pos << std::endl;
//    return  std::distance(first1, pair.first);
//}
long offset(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2, std::size_t start_pos, Mutations &mutations) {
    std::pair<std::string::const_iterator, std::string::const_iterator> pair;
    
    Mutations m;
    
    pair = std::mismatch(first1, last1, first2, [&m, &start_pos](char a, char b) {
        if (a=='N' || b=='N') {
            return false;
        } else if (a!=b) {
            m.push_back(start_pos);
        }
        start_pos++;
        return m.size()<2;
    });
    
    m.swap(mutations);
    return std::distance(first1, pair.first);
}

long offset(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, std::size_t start_pos, Mutations& mutations) {
    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
    
    Mutations m;
    
    pair = std::mismatch(first1, last1, first2, [&m,&start_pos](char a, char b) {
        if (a=='N' || b=='N') {
            return false;
        } else if (a!=b) {
            m.push_back(start_pos);
        }
        start_pos--;
        return m.size()<2;
    });
    
    //because we're expanding to the left, second accuring mutation should be in the first position.
    if (m.size() > 1)
        std::swap(m[0], m[1]);
    
    m.swap(mutations);
    
    return  std::distance(first1, pair.first);
}

//long offset(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, int &mut_pos) {
//    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
//    int count = 0;
//    int pos = 0;
//    int tmp = 0;
//    pair = std::mismatch(first1, last1, first2, [&count,&pos,&tmp](char a, char b) {
//        switch (a) {
//            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
//                if (a!=b) {
//                    count++;
//                    
//                    if (count<2) tmp=pos;
//                }
//                
//                pos++;
//                return count<2 ;
//            default:
//                tmp=pos;
//                return false;
//        }
//    });
//    mut_pos=tmp;
//    return  std::distance(first1, pair.first);
//}

//long offset(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, const bool allow_mut) {
//    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
//    int count = allow_mut? 0 : 1;
//    pair = std::mismatch(first1, last1, first2, [&count](char a, char b){
//        switch (a) {
//            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
//                if (a!=b) count++;
//                return count<2 ;
//            default:
//                return false;
//        }
//    });
//    return  std::distance(first1, pair.first);
//}

#pragma mark - Fine grain filter
std::size_t offset_test_forward(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2, std::size_t start_pos, Mutations &mutations) {
    std::pair<std::string::const_iterator, std::string::const_iterator> pair;
    
    Mutations m;
    
    pair = std::mismatch(first1, last1, first2, [&m, &start_pos](char a, char b) {
        if (a=='N' || b=='N') {
            return false;
        } else if (a!=b) {
            m.push_back(start_pos);
        }
        start_pos++;
        return m.size()<2;
    });
    
    m.swap(mutations);
    //return std::distance(first1, pair.first);
    
    return start_pos;
}

std::size_t offset_test_reverse(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, std::size_t start_pos, Mutations& mutations) {
    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
    
    Mutations m;
    
    pair = std::mismatch(first1, last1, first2, [&m,&start_pos](char a, char b) {
        if (a=='N' || b=='N') {
            return false;
        } else if (a!=b) {
            m.push_back(start_pos);
        }
        start_pos--;
        return m.size()<2;
    });
    
    //because we're expanding to the left, second accuring mutation should be in the first position.
    if (m.size() > 1)
        std::swap(m[0], m[1]);
    
    m.swap(mutations);
    
    //return  std::distance(first1, pair.first);
    return start_pos;
}

//std::ofstream debug("/Users/andi/Applications/Research/Projects/Limes/Xcode/Limes/data/test/test_keys/output/debug_12/debug_12.limes.txt");
/**
 Phase II - Fine grain filter
    d1 = query
    d2 = target
 */
void remove_invalid_lime_candidates(Candidates &candidates, const Query &d1, const Target &d2) {
    if (candidates.empty()) return;
    
    Candidates::size_type length = candidates.size();
    const Sequence::size_type size1 = d1.length(), size2 = d2.length();
    
    std::string::const_iterator d1_begin, d1_end=d1.end(), d2_begin, d2_end = d2.end();
    std::string::const_reverse_iterator d1_rbegin, d1_rend=d1.rend(), d2_rbegin, d2_rend = d2.rend();
    
    Mutations left_mutations, right_mutations;
    for (Candidates::size_type i = 0; i < length; ++i) {
        Candidate &c = candidates[i];
        const std::size_t idx_1 = c.idx_1;
        const std::size_t idx_2 = c.idx_2;
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (left to right)        //
        //////////////////////////////////////////////////////////
        d1_begin = d1.begin();
        d2_begin = d2.begin();
        std::advance(d1_begin, idx_1);
        std::advance(d2_begin, idx_2);
        
        const std::size_t right_most_offset = (size1-idx_1 < size2-idx_2) ? offset_test_forward(d1_begin, d1_end, d2_begin, idx_1, right_mutations) : offset_test_forward(d2_begin, d2_end, d1_begin, idx_1, right_mutations);
        if (right_most_offset < SEQLENGTH/2 - WORDSIZE/2) continue;
 
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (right to left)        //
        //////////////////////////////////////////////////////////
        d1_rbegin = d1.rbegin();
        d2_rbegin = d2.rbegin();
        std::advance(d1_rbegin, size1-idx_1);
        std::advance(d2_rbegin, size2-idx_2);
        const std::size_t left_most_offset = (idx_1 < idx_2) ? offset_test_reverse(d1_rbegin, d1_rend, d2_rbegin, idx_1, left_mutations) : offset_test_reverse(d2_rbegin, d2_rend, d1_rbegin, idx_1, left_mutations);
        
        std::string case_region = "";
        std::size_t max_range = 0;
        std::size_t query_start = 0, target_start = 0;
        if (left_mutations.size()+right_mutations.size() < 2) {                         //case 1,2,3
//            std::cout << "cases 1,2,3" << std::endl;
            case_region = "cases 1,2,3";
            
            query_start = left_most_offset;
            target_start = idx_2-(idx_1-query_start);
            max_range = right_most_offset-left_most_offset-1;
            
        } else if (right_mutations.size() + left_mutations.size() > 3) {                //case 9
            const std::size_t range1 = right_mutations[0] - left_mutations[0];
            const std::size_t range2 = right_mutations[1] - left_mutations[1];
            
            if (range1 >  range2) {
//                std::cout << "case 9, range1" << std::endl;
                case_region = "case 9, range1";
                
                max_range = range1;
                query_start = left_mutations[0];
                target_start = idx_2-(idx_1-left_mutations[0]);
                
            } else {
//                std::cout << "case 9, range2" << std::endl;
                case_region = "case 9, range2";
                
                max_range = range2;
                query_start = left_mutations[1];
                target_start = idx_2-(idx_1-left_mutations[1]);
            }
        } else if (right_mutations.size() == 1 && left_mutations.size() == 1) {         //case 4
            const int range1 = static_cast<int>(right_most_offset - left_mutations[0]);
            const int range2 = static_cast<int>(right_mutations[0] - left_most_offset);
            
            if (range1 > range2) {
//                std::cout << "case 4, region 1" << std::endl;
                case_region = "case 4, region 1";
                
                max_range = range1;
                query_start = left_mutations[0];
                target_start = target_start = idx_2-(idx_1-left_mutations[0]);
                
            } else {
//                std::cout << "case 4, region 2" << std::endl;
                case_region = "case 4, region 2";
                
                max_range = range2;
                query_start = left_most_offset;
                target_start = target_start = idx_2-(idx_1-left_most_offset);
            }
            
        } else if (right_mutations.size() + left_mutations.size() > 2) {                //case 7,8
            
            if (right_mutations.size() > 1) {
                
                const std::size_t range1 = static_cast<int>(right_mutations[0] - left_most_offset);
                const std::size_t range2 = right_mutations[1] - left_mutations[0];
                
                if (range1 > range2) {
//                    std::cout << "case 7 range1" << std::endl;
                    
                    case_region = "case 7 range1";
                    
                    max_range = range1;
                    query_start = left_most_offset;
                    target_start = target_start = idx_2-(idx_1-left_most_offset);
                    
                } else {
//                    std::cout << "case 7 range2" << std::endl;
                    case_region = "case 7 range2";
                    
                    max_range = range2;
                    query_start = left_mutations[0];
                    target_start = idx_2-(idx_1-left_mutations[0]);
                }
                
            } else {
                const std::size_t range1 = right_mutations[0] - left_mutations[0];
                const std::size_t range2 = static_cast<int>(right_most_offset - left_mutations[1]);
                
                if (range1 > range2) {
//                    std::cout << "case 8 range1" << std::endl;
                    
                    case_region = "case 8 range1";
                    
                    max_range = range1;
                    query_start = left_mutations[0];
                    target_start = target_start = idx_2-(idx_1-left_mutations[0]);
                    
                } else {
//                    std::cout << "case 8 range2" << std::endl;
                    case_region = "case 8 range2";
                    max_range = range2;
                    query_start = left_mutations[1];
                    target_start = idx_2-(idx_1-left_mutations[1]);
                }
            }
        } else {                                                                        //case 5,6
            
            if (right_mutations.empty()) {                                              //case 5
//                std::cout << "case 5" << std::endl;
                case_region = "case 5";
                
                max_range = static_cast<int>(right_most_offset-left_mutations[0]);
                query_start = left_mutations[0];
                target_start = target_start = idx_2-(idx_1-left_mutations[0]);
                
            } else {                                                                    //case 6
//                std::cout << "case 6" << std::endl;
                case_region = "case 6";
                
                max_range = static_cast<int>(right_mutations[1] - left_most_offset);
                query_start = left_most_offset;
                target_start = target_start = idx_2-(idx_1-left_most_offset);
            }
        }
        
        if (max_range  >= SEQLENGTH) {
            //std::cout << case_region << std::endl;
            //std::cout << "left_most_offset: " << left_most_offset << ", right_most_offset: " << right_most_offset << std::endl;
            //std::cout << "left_mutations[0]: " << left_mutations[0] << ", left_mutations[1]: " << left_mutations[1] << std::endl;
            //std::cout << "right_mutations[0]: " << right_mutations[0] << ", right_mutations[1]:" << right_mutations[1] << std::endl;

            //std::cout << "query range: " << max_range << ", start pos: " << query_start << ", idx 1: " << idx_1 << std::endl;
            //std::cout << "target range: " << max_range << ",start pos: " << target_start << ", idx 2: " << idx_2 << std::endl;
            
            //debug << "query start " << query_start << ", length " <<  max_range << std::endl;
            //debug << "target start " << target_start << ", length " <<  max_range << std::endl;

            //debug << "query:" << d1.substr(query_start, max_range) << std::endl;
            //debug << "target:" << d2.substr(target_start, max_range) << std::endl << std::endl;
        }
        
        if (max_range  >= SEQLENGTH) {
            limeObjs_target.push_back(std::make_pair(target_start,max_range));
            limeObjs_query.push_back(std::make_pair(query_start,max_range));
        }
    }
}

void remove_invalid_lime_candidates_reverse(Candidates &candidates, const Query &d1, const Target &d2) {
    if (candidates.empty()) return;
    
    Candidates::size_type length = candidates.size();
    const Sequence::size_type size1 = d1.length(), size2 = d2.length();
    
    std::string::const_iterator d1_begin, d1_end=d1.end(), d2_begin, d2_end = d2.end();
    std::string::const_reverse_iterator d1_rbegin, d1_rend=d1.rend(), d2_rbegin, d2_rend = d2.rend();
    
    Mutations left_mutations, right_mutations;
    for (Candidates::size_type i = 0; i < length; ++i) {
        Candidate &c = candidates[i];
        const std::size_t idx_1 = c.idx_1;
        const std::size_t idx_2 = c.idx_2;
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (left to right)        //
        //////////////////////////////////////////////////////////
        d1_begin = d1.begin();
        d2_begin = d2.begin();
        std::advance(d1_begin, idx_1);
        std::advance(d2_begin, idx_2);
        
        const std::size_t right_most_offset = (size1-idx_1 < size2-idx_2) ? offset_test_forward(d1_begin, d1_end, d2_begin, idx_1, right_mutations) : offset_test_forward(d2_begin, d2_end, d1_begin, idx_1, right_mutations);
        if (right_most_offset < SEQLENGTH/2 - WORDSIZE/2) continue;
        
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (right to left)        //
        //////////////////////////////////////////////////////////
        d1_rbegin = d1.rbegin();
        d2_rbegin = d2.rbegin();
        std::advance(d1_rbegin, size1-idx_1);
        std::advance(d2_rbegin, size2-idx_2);
        const std::size_t left_most_offset = (idx_1 < idx_2) ? offset_test_reverse(d1_rbegin, d1_rend, d2_rbegin, idx_1, left_mutations) : offset_test_reverse(d2_rbegin, d2_rend, d1_rbegin, idx_1, left_mutations);
        
        std::string case_region = "";
        std::size_t max_range = 0;
        std::size_t query_start = 0, target_start = 0;
        if (left_mutations.size()+right_mutations.size() < 2) {                         //case 1,2,3
            //            std::cout << "cases 1,2,3" << std::endl;
            case_region = "cases 1,2,3";

            //size - (idx+offset-idx)
            query_start = d1.length() - (idx_1+right_most_offset-idx_1);
            target_start = idx_2-(idx_1-left_most_offset);
            max_range = right_most_offset-left_most_offset-1;
            
        } else if (right_mutations.size() + left_mutations.size() > 3) {                //case 9
            const std::size_t range1 = right_mutations[0] - left_mutations[0];
            const std::size_t range2 = right_mutations[1] - left_mutations[1];
            
            if (range1 >  range2) {
                //                std::cout << "case 9, range1" << std::endl;
                case_region = "case 9, range1";
                
                //size - (idx+offset-idx)
                query_start = d1.length() - (idx_1+right_mutations[0]-idx_1);
                target_start = idx_2-(idx_1-left_mutations[0]);
                max_range = range1;
                
            } else {
                //                std::cout << "case 9, range2" << std::endl;
                case_region = "case 9, range2";
                
                //size - (idx+offset-idx)
                query_start = d1.length() - (idx_1+right_mutations[1]-idx_1);
                target_start = idx_2-(idx_1-left_mutations[1]);
                max_range = range2;
                
            }
        } else if (right_mutations.size() == 1 && left_mutations.size() == 1) {         //case 4
            const std::size_t range1 = right_most_offset - left_mutations[0];
            const std::size_t range2 = right_mutations[0] - left_most_offset;
            
            if (range1 > range2) {
                //                std::cout << "case 4, region 1" << std::endl;
                case_region = "case 4, region 1";
                
                //size - (idx+offset-idx)
                query_start = d1.length() - (idx_1+right_most_offset-idx_1);
                target_start = idx_2-(idx_1-left_mutations[0]);
                max_range = range1;
                
            } else {
                //                std::cout << "case 4, region 2" << std::endl;
                case_region = "case 4, region 2";
                
                //size - (idx+offset-idx)
                query_start = d1.length() - (idx_1+right_mutations[0]-idx_1);
                target_start = idx_2-(idx_1-left_most_offset);
                max_range = range2;
            }
            
        } else if (right_mutations.size() + left_mutations.size() > 2) {                //case 7,8
            
            if (right_mutations.size() > 1) {                                           //case 7
                
                const std::size_t range1 = right_mutations[0] - left_most_offset;
                const std::size_t range2 = right_mutations[1] - left_mutations[0];
                
                if (range1 > range2) {
                    //                    std::cout << "case 7 range1" << std::endl;
                    
                    case_region = "case 7 range1";
                    
                    //size - (idx+offset-idx)
                    query_start = d1.length() - (idx_1+right_mutations[0]-idx_1);
                    target_start = idx_2-(idx_1-left_most_offset);
                    max_range = range1;
                    
                } else {
                    //                    std::cout << "case 7 range2" << std::endl;
                    case_region = "case 7 range2";
                    
                    //size - (idx+offset-idx)
                    query_start = d1.length() - (idx_1+right_mutations[1]-idx_1);
                    target_start = idx_2-(idx_1-left_mutations[0]);
                    max_range = range2;
                }
                
            } else {                                                                  //case 8
                const std::size_t range1 = right_mutations[0] - left_mutations[0];
                const std::size_t range2 = right_most_offset - left_mutations[1];
                
                if (range1 > range2) {
                    //                    std::cout << "case 8 range1" << std::endl;
                    
                    case_region = "case 8 range1";
                    
                    //size - (idx+offset-idx)
                    query_start = d1.length() - (idx_1+right_mutations[0]-idx_1);
                    target_start = target_start = idx_2-(idx_1-left_mutations[0]);
                    max_range = range1;
                    
                } else {
                    //                    std::cout << "case 8 range2" << std::endl;
                    case_region = "case 8 range2";
                    
                    //size - (idx+offset-idx)
                    query_start = d1.length() - (idx_1+right_most_offset-idx_1);
                    target_start = idx_2-(idx_1-left_mutations[1]);
                    max_range = range2;
                }
            }
        } else {                                                                        //case 5,6
            
            if (right_mutations.empty()) {                                              //case 5
                //                std::cout << "case 5" << std::endl;
                case_region = "case 5";
                
                //size - (idx+offset-idx)
                query_start = d1.length() - (idx_1+right_most_offset-idx_1);
                target_start = target_start = idx_2-(idx_1-left_mutations[0]);
                max_range = right_most_offset-left_mutations[0];
                
            } else {                                                                    //case 6
                //                std::cout << "case 6" << std::endl;
                case_region = "case 6";
                
                //size - (idx+offset-idx)
                query_start = d1.length() - (idx_1+right_mutations[1]-idx_1);
                target_start = idx_2-(idx_1-left_most_offset);
                max_range = right_mutations[1] - left_most_offset;
            }
        }
        
        if (max_range  >= SEQLENGTH) {
            //std::cout << case_region << std::endl;
            //std::cout << "left_most_offset: " << left_most_offset << ", right_most_offset: " << right_most_offset << std::endl;
            
            //std::cout << "query range: " << max_range << ", start pos: " << query_start << ", idx 1: " << idx_1 << std::endl;
            //std::cout << "target range: " << max_range << ", start pos: " << target_start << ", idx 2: " << idx_2 << std::endl;
            
//            Query rc = d1;
//            reverse_complement(rc);
//            debug << "query start " << query_start << ", length " <<  max_range << std::endl;
//            debug << "target start " << target_start << ", length " <<  max_range << std::endl;
            
//            debug << "target seq: " << d2.substr(target_start, max_range) << std::endl;
//            debug << "query: " << rc.substr(query_start, max_range) << std::endl << std::endl;
            //std::cout << "query:" << d1.substr(query_start, max_range) << std::endl;
        }
        
        
        if (max_range  >= SEQLENGTH) {
            limeObjs_target.push_back(std::make_pair(target_start,max_range));
            limeObjs_query.push_back(std::make_pair(query_start,max_range));
        }
    }
}

#pragma mark - Coarse grain filter
/**
    Phase I - Coarse grain filter
 */
void find_candidates(const Hash seqA_begin, const Hash seqA_end, const std::size_t index, Candidates &candidates);
void generate_lime_candidates(const std::size_t i, const Chromosome &query_chr, Candidates &candidates ) {
    
    const Hash seqA_begin = query_chr[i];
    const Hash seqA_end = query_chr[i+1];
    if (seqA_begin < 0 || seqA_end < 0)
        return;
    
    const Pos & vec2DPos = vec2D[seqA_begin];
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    const std::size_t length = vec2DPos.size();
    
    const std::bitset<WORDSIZE*2> endA(seqA_end);
    
    /**if (i*chunk == 84743296) {
        std::cout << "query begin hash: " << seqA_begin << std::endl;
        std::cout << "query end hash: " << seqA_end << std::endl;
    }
    */
    for (int j = 0; j < length; ++j) {
        const Element & element = vec2DPos[j];
        const Hash seqB_end = element.id;
        
        if (seqB_end < 0) continue; //skip anything flagged for mutations in the prefix
        
        std::bitset<WORDSIZE*2> endB(seqB_end);
        
        
      /*  if (i*chunk == 84743296 && seqB_end == 84743340) {
            std::cout << "target end hash: " << seqB_end << ", target end hash location: " << element.idx << std::endl;
        }
        */
        
//        if (std::abs(i*chunk-14756104) < 20 && seqA_begin == 2105376 && seqA_end == 2631720 && seqB_end == seqA_end) {
//            std::cout << "pos: " << i*chunk << ", motif hash: " << seqA_begin << std::endl;
//        }
        
        // Perform exclusive OR to get the differences between the two keys
        endB ^= endA;
        
        //////////////////////////////////////////////////////////////////////////
        // One character occupies two bits, hence less than 3
        // Increment by 2 to avoid comparing bits from different characters
        //////////////////////////////////////////////////////////////////////////
        if (endB.none()) {
            candidates.push_back(Candidate(i*chunk, element.idx));
            //if (i*chunk == 84743296)
                //std::cout << "i: " << i << ", id: " << element.id << ", idx: " << element.idx << std::endl;
        }
        else if(endB.count() < 3) {
            int count = 0;
            for (int k = 0; k < WORDSIZE*2; k+=2) if (endB[k] || endB[k+1]) count++;
            
            if(count < 2) candidates.push_back(Candidate(i*chunk, element.idx));
        }
    }
    
    find_candidates(seqA_begin, seqA_end, i*chunk, candidates);
}

//0.
//S1 = {1	"3	10"	15}             - query
//S2 = {5	"3	10"	7	"4	10"	23} - target

//1. get subset using 10 as the index giving subset S2={3,7,4,23}
//2. ignore keys that are identical to 3 from S giving S2={7,4,23}
//3. do bitwise comparison between 3 and 7,4,23 and select those that are at most 2 bit different - in this case giving 4.
//4. calculate index for 3
//5. index for 4 should be idx.
void find_candidates(const Hash seqA_begin, const Hash seqA_end, const std::size_t index, Candidates &candidates) {
    //1,2
    Pos vec2DPos = vec2D[seqA_end];
    vec2DPos.erase(std::remove_if(vec2DPos.begin(), vec2DPos.end(), [](const Element &e){
        return e.id > 0;
    }), vec2DPos.end());
    
    
    const std::bitset<WORDSIZE*2> begin_A(seqA_begin);
    for (int i = 0; i < vec2DPos.size(); ++i) {
        const Element e = vec2DPos[i];
        std::bitset<WORDSIZE*2> begin_B(e.id * -1);
        
        //3.
        begin_B ^= begin_A;
        if(begin_B.count() < 3) {
            int count = 0;
            for (int k = 0; k < WORDSIZE*2; k+=2) if (begin_B[k] || begin_B[k+1]) count++;
            
            //4,5
            if(count < 2) candidates.push_back(Candidate(index, e.idx));
        }
    }
}

void find_lime_candidates(Chromosome &query_chr, Candidates &candidates) {
    
    const std::size_t iterations = query_chr.size()-1;
    for (std::size_t i = 0; i < iterations; ++i) {
        generate_lime_candidates(i, query_chr, candidates);
    }
}

void find_limes_in_forward_direction(Chromosome &query_chr, const Query &query_data, const Target &target_data, Candidates &candidates) {
    //////////////////////////
    // Find candidate limes //
    //////////////////////////
    find_lime_candidates(query_chr, candidates);
    
    //////////////////////////////////
    // Remove invalid candidates    //
    //////////////////////////////////
    remove_invalid_lime_candidates(candidates, query_data, target_data);
}

void find_limes_in_reverse_direction(Chromosome &query_chr, const Query &query_data, const Target &target_data, Candidates &candidates) {
    //////////////////////////
    // first stage filter   //
    //////////////////////////
    find_lime_candidates(query_chr, candidates);
    
    //////////////////////////
    // second stage filter  //
    //////////////////////////
    remove_invalid_lime_candidates_reverse(candidates, query_data, target_data);
}


/**
 target    = does not change
 query     = we scan it in the forward direction, then in the reverse direction
 a         = generated by target?
*/
void find_limes(const Target &target_data, const Query &query_data, const Chromosome& target_chr) {
    //debug << "#Forward" << std::endl;
    //skips every chunk size (SEQLENGTH/2 - WORDSIZE/2)
    Sequence query_data_cloned (query_data);
    Chromosome query_chr = generateVectorIndices(query_data_cloned);
    
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
    reverse_complement(query_data_cloned);
    
    //////////////////////////
    // reverse scan         //
    //////////////////////////
    //debug << "#Reverse" << std::endl;
    generateVectorIndecies(query_data_cloned, query_chr);
    find_limes_in_reverse_direction(query_chr, query_data_cloned, target_data, candidates);
    //debug.close();
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

#pragma mark - Start of algorithm
void run(const Files &g1, const Files &g2, const Output &pathToLimes, const Progress &pathToProgress) {
    std::cout << "Length search key: " << WORDSIZE << std::endl;
    
    timer.start();
    const bool shouldSwap = false;//optimize(g1, g2);
    std::vector<std::string>genome1 = g1, genome2 = g2;
    if (shouldSwap)
        genome1.swap(genome2);
    
    const int totalComparisons = static_cast<int>(g1.size()*g2.size());
    const Files::size_type g1Size = genome1.size();
    const Files::size_type g2Size = genome2.size();
    
    limeObjs_target.reserve(100000);
    limeObjs_query.reserve(1000000);
    Chromosome target_chr, query_chr;
    Sequence target_data, query_data;
    
    Header targetSeqHeader, querySeqHeader;
    
    std::ofstream out (pathToLimes.c_str());
    std::ofstream progress (pathToProgress.c_str());
    
    int counter = 1;
    
    //note: outter loop should be the larger raw data size
    for (Files::size_type i = 0; i < g1Size; ++i) {
        const File targetSeqFile(genome1[i]);
        
        //////////////////////////////////////////////
        // Skips by 1 letter
        // goal: Need to minimize vec2D generation
        //////////////////////////////////////////////
        //timer.split();
        target_data = loadDataWithContentsOFile(targetSeqFile, targetSeqHeader);
        target_chr = generateLookupTableIndices(target_data);
        initializeLookupTable(target_chr);
        //std::cout << "Load and lookup table generation time: " << timer.getSplitElapsedTime() << std::endl;
        
        for (Files::size_type j = 0; j < g2Size; ++j) {
            const File querySeqFile(genome2[j]);
            
            timer.split();
            query_data = loadDataWithContentsOFile(querySeqFile, querySeqHeader);
            progress << "Processing..." << counter << "/" << totalComparisons << "\n";
            progress << targetSeqFile << "\n" << querySeqFile << std::endl;
            
            find_limes(target_data,query_data,target_chr);
            
            progress << "Processing time = " << timer.getSplitElapsedTime() << "\n" << std::endl;
            
            //////////////////////////////////////
            // Write limes to file
            //////////////////////////////////////
            if (!limeObjs_target.empty()) {
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Limes::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Limes::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
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
    progress << "Total time = " << timer.getTotalElapsedTime() << std::endl;
    out.close();
    progress.close();
    timer.stop();
}

void run(const File &file, const Files &g, const File &limes, const Progress &progress_file) {
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
    
    Sequence query_sequence, target_sequence;
    Header target_header;
    
    std::ofstream out (limes.c_str());
    std::ofstream progress (progress_file.c_str());
    for (std::vector<std::string>::size_type i = 0; i < number_of_files; ++i) {
        const File target_file(g[i]);
        
        //////////////////////////////////////////////
        // Skips by 1 letter
        // goal: Need to minimize vec2D generation
        //////////////////////////////////////////////
        //timer.split();
        target_sequence = loadDataWithContentsOFile(target_file, target_header);
        const Chromosome target_chr = generateLookupTableIndices(target_sequence);
        initializeLookupTable(target_chr, target_sequence);
        //std::cout << "load and lookup table generation time: " << timer.getSplitElapsedTime() << std::endl;
        
        timer.split();
        for (std::vector<std::pair<std::string, std::string> >::size_type n = 0; n < number_of_query_seq; ++n) {
            const std::string query_header (query_sequences[n].first);
            query_sequence = query_sequences[n].second;
            
            find_limes(target_sequence,query_sequence,target_chr);
            
            //////////////////////////////////////
            // Write limes to file
            //////////////////////////////////////
            if (!limeObjs_target.empty()) {
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Limes::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Limes::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
                out << "#1" << query_header << std::endl;
                out << "#start_1" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_query.begin(), q_it, [&out,&query_sequence](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                out << "#2" << target_header << std::endl;
                out << "#start_2" << "\t" << "length" << std::endl;
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

void run(const File &queryGenome, const File &targetG2, const File &limes, const File &progress_file) {
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
    
    Sequence query_sequence, target_sequence;
    std::string target_header;
    
    std::ofstream out (limes.c_str());
    std::ofstream progress (progress_file.c_str());
    
    scottgs::Timing outputTimer;
    outputTimer.start();
    for (std::vector<std::pair<std::string, std::string> >::size_type i = 0; i < number_of_target_seq; ++i) {
        target_header = target_sequences[i].first;
        const Sequence target_sequence = target_sequences[i].second;
        
        
        //////////////////////////////////////////////
        // Skips by 1 letter
        // goal: Need to minimize vec2D generation
        //////////////////////////////////////////////
        const Chromosome target_chr = generateLookupTableIndices(target_sequence);
        initializeLookupTable(target_chr, target_sequence);
        
        timer.split();
        for (std::vector<std::pair<std::string, std::string> >::size_type n = 0; n < number_of_query_seq; ++n) {
            const std::string query_header (query_sequences[n].first);
            query_sequence = query_sequences[n].second;
            
            find_limes(target_sequence,query_sequence,target_chr);
            
            //////////////////////////////////////
            // Write limes to file
            //////////////////////////////////////
            if (!limeObjs_target.empty()) {
                outputTimer.split();
                assert(limeObjs_target.size() == limeObjs_query.size());
                std::sort(limeObjs_query.begin(), limeObjs_query.end());
                std::sort(limeObjs_target.begin(), limeObjs_target.end());
                const Limes::iterator  q_it = std::unique(limeObjs_query.begin(), limeObjs_query.end());
                const Limes::iterator  t_it = std::unique(limeObjs_target.begin(), limeObjs_target.end());
                
                out << "#1" << query_header << std::endl;
                out << "#start_1" << "\t" << "length" << std::endl;
                std::for_each(limeObjs_query.begin(), q_it, [&out,&query_sequence](const std::pair<std::string::size_type, int> &l){
                    out << l.first << "\t" << l.second << std::endl;
                });
                out << "#2" << target_header << std::endl;
                out << "start_2" << "\t" << "length" << std::endl;
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

#pragma mark -


