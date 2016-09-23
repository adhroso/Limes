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
	assert(len<=sizeof(unsigned int)*4);
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

void generateVectorIndecies(const Sequence &data, Chromosome &c) {
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
    int size = static_cast<int>(pow(4, WORDSIZE));
   
    //initialize lookuptable - holds location(s) of each possible word
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
            pos.push_back(Element(end, i));      //index (i) refers to the beginning
            
            Pos & pos2 = vec2D[end];
            pos2.push_back(Element(begin,i));
        }
    }
}

long offset(std::string::const_iterator first1, std::string::const_iterator last1, std::string::const_iterator first2, int &mut_pos) {
    std::pair<std::string::const_iterator, std::string::const_iterator> pair;
    int count = 0;
    int pos = 0;
    int tmp = 0;
    pair = std::mismatch(first1, last1, first2, [&count, &pos, &tmp](char a, char b){
        switch (a) {
            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
                if (a!=b) {
                    count++;
                    
                    if (count<2) tmp=pos;
                }
                
                pos++;
                return count<2 ;
            default:
                return false;
        }
    });
    mut_pos=tmp;
    return  std::distance(first1, pair.first);
}

long offset(std::string::const_reverse_iterator first1, std::string::const_reverse_iterator last1, std::string::const_reverse_iterator first2, const bool allow_mut) {
    std::pair<std::string::const_reverse_iterator, std::string::const_reverse_iterator> pair;
    int count = allow_mut? 0 : 1;
    pair = std::mismatch(first1, last1, first2, [&count](char a, char b){
        switch (a) {
            case 'A':   case 'C':   case 'G':   case 'T': case 'a':   case 'c':   case 'g':   case 't':
                if (a!=b) count++;
                return count<2 ;
            default:
                return false;
        }
    });
    return  std::distance(first1, pair.first);
}

#pragma mark - Fine grain filter
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

        const int queryDistanceToEnd = static_cast<int>(size1-idx_1);
        const int targetDistanceToEnd = static_cast<int>(size2-idx_2);
        
        int mut_pos = 0;
        const long right_offset = (queryDistanceToEnd < targetDistanceToEnd) ? offset(d1_begin, d1_end, d2_begin,mut_pos) : offset(d2_begin,d2_end,d1_begin,mut_pos);
        if (right_offset < SEQLENGTH/2 - WORDSIZE/2)
            continue;
        
        //////////////////////////////////////////////////////////
        // scan in the forward direction (right to left)        //
        //////////////////////////////////////////////////////////
        d1_rbegin = d1.rbegin();
        d2_rbegin = d2.rbegin();
        std::advance(d1_rbegin, size1-idx_1);
        std::advance(d2_rbegin, size2-idx_2);

        const bool ALLOW_MUT = mut_pos==right_offset ? true : false;
        const long left_offset = (idx_1 < idx_2) ? offset(d1_rbegin, d1_rend, d2_rbegin,ALLOW_MUT) : offset(d2_rbegin, d2_rend, d1_rbegin,ALLOW_MUT);
        
        const long diff = right_offset+left_offset;
        if (diff  >= SEQLENGTH) {
            limeObjs_target.push_back(std::make_pair(idx_2-left_offset,diff));
            limeObjs_query.push_back(std::make_pair(idx_1-left_offset,diff));
            
            if(idx_2-left_offset == 11047460 && diff == 113) {
                std::cout << "key poss:" << idx_2 << ", seq start pos: " << idx_2-left_offset << ", seq: " << d2.substr(idx_2-left_offset, diff) << std::endl;
                //180715756
            }
        }
    }
}

void remove_invalid_lime_candidates_reverse(Candidates &candidates, const Query &query_data, const Target &target_data) {
    if (candidates.empty()) return;
    
    Candidates::size_type length = candidates.size();
    const Sequence::size_type size1 = query_data.length(), size2 = target_data.length();
    
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

        const int queryDistanceToEnd = static_cast<int>(size1-idx_1);
        const int targetDistanceToEnd = static_cast<int>(size2-idx_2);
        
        int mut_pos = 0;
        const long right_offset = (queryDistanceToEnd < targetDistanceToEnd) ? offset(d1_begin, d1_end, d2_begin,mut_pos) : offset(d2_begin, d2_end, d1_begin,mut_pos);

        const bool ALLOW_MUT = mut_pos==right_offset ? true : true;
        
        if (idx_2 == 11047462 && idx_1 == 14756104) {
            std::cout << "\nRight offset: " << right_offset << std::endl;
            std::cout << "allow mutations: " << (ALLOW_MUT ? "yes" : "no") << ", mut pos: " << mut_pos << std::endl;
        }
        
        if (right_offset < SEQLENGTH/2 - WORDSIZE/2)
            continue;
        
        //////////////////////////////////////////////////////////
        // scan in the reverse direction (right to left)        //
        //////////////////////////////////////////////////////////
        d1_rbegin = query_data.rbegin();
        d2_rbegin = target_data.rbegin();
        std::advance(d1_rbegin, size1-idx_1);
        std::advance(d2_rbegin, size2-idx_2);

        const long left_offset = (idx_1 < idx_2) ? offset(d1_rbegin, d1_rend, d2_rbegin,ALLOW_MUT) : offset(d2_rbegin, d2_rend, d1_rbegin,ALLOW_MUT);
        if (idx_2 == 11047462 && idx_1 == 14756104) {
            std::cout << "left offset: " << left_offset << std::endl << std::endl ;
        }
        
        
        const long diff = right_offset+left_offset;
        if (diff  >= SEQLENGTH) {
            const Sequence::size_type query_data_size = query_data.size();
            const Sequence::size_type query_start = query_data_size - (idx_1+right_offset);
            limeObjs_target.push_back(std::make_pair(idx_2-left_offset, diff));
            limeObjs_query.push_back(std::make_pair(query_start, diff));
            
            if(idx_2-left_offset == 11047460 && diff == 113) {
                std::cout << "found in reverse complement" << std::endl;
                std::cout << "human: key pos:" << idx_2 << ", seq start pos: " << idx_2-left_offset << ", seq: " << target_data.substr(idx_2-left_offset, diff) << std::endl;
                std::cout << "mouse: key pos:" << idx_1 << ", seq start pos: " << idx_1-left_offset << ", seq: " << query_data.substr(idx_1-left_offset, diff) << std::endl << std::endl;
                
                bool valid = false;
                std::cout << "human: hash first key pos:" << idx_2 << ", pos of key in lime:" << idx_2-(idx_2-left_offset) << ", motif: " << target_data.substr(idx_2, WORDSIZE) << ", motif hash: "<< encodeWord((const u_char*)target_data.substr(idx_2, WORDSIZE).c_str(), WORDSIZE, valid) << std::endl;
                std::cout << "mouse: hash first key pos:" << idx_1 << ", pos of key in lime:" << idx_1-(idx_1-left_offset) << ", motif: " << query_data.substr(idx_1, WORDSIZE) << ", motif hash: " << encodeWord((const u_char*)query_data.substr(idx_1, WORDSIZE).c_str(), WORDSIZE, valid) << std::endl << std::endl;
                
                const int chunk = SEQLENGTH/2 - WORDSIZE/2;
                std::cout << "human: hash second key pos:" << idx_2+(SEQLENGTH/2 - WORDSIZE/2) << ", pos of key in lime:" <<  (idx_2+chunk)-(idx_2-left_offset) << ", motif: " << target_data.substr(idx_2+chunk, WORDSIZE) << ", motif hash: " << encodeWord((const u_char*)target_data.substr(idx_2+chunk, WORDSIZE).c_str(), WORDSIZE, valid) << std::endl;
                std::cout << "mouse: hash second key pos:" << idx_1+(SEQLENGTH/2 - WORDSIZE/2) << ", pos of key in lime:" << (idx_1+chunk)-(idx_1-left_offset) <<  ", motif: " << query_data.substr(idx_1+chunk, WORDSIZE) << ", motif hash: " << encodeWord((const u_char*)query_data.substr(idx_1+chunk, WORDSIZE).c_str(), WORDSIZE, valid) << std::endl << std::endl;

                
                
                std::cout << "human: first key motif: "  << target_data.substr(idx_2, WORDSIZE+2) << ", motif hash: " << encodeWord((const u_char*)target_data.substr(idx_2, WORDSIZE+2).c_str(), WORDSIZE+2, valid) << std::endl;
                std::cout << "mouse: first key motif: " <<  query_data.substr(idx_1, WORDSIZE+2) << ", motif hash: " << encodeWord((const u_char*)query_data.substr(idx_1, WORDSIZE+2).c_str(), WORDSIZE+2, valid) << std::endl << std::endl;

                std::cout << "human: hash second key motif: " << target_data.substr(idx_2+chunk, WORDSIZE+2) << ", motif hash: " << encodeWord((const u_char*)target_data.substr(idx_2+chunk, WORDSIZE+2).c_str(), WORDSIZE+2, valid) << std::endl;
                std::cout << "mouse: hash second key motif: " << query_data.substr(idx_1+chunk, WORDSIZE+2) << ", motif hash: " << encodeWord((const u_char*)query_data.substr(idx_1+chunk, WORDSIZE+2).c_str(), WORDSIZE+2, valid) << std::endl << std::endl;
                
            }
         }
    }
}

#pragma mark - Coarse grain filter
/**
    Phase I - Coarse grain filter
 */
void find_candidates(const int seqA_begin, const int seqA_end, const int index, Candidates &candidates);
void generate_lime_candidates(const int i, const Chromosome &query_chr, Candidates &candidates ) {
    
    const int seqA_begin = query_chr[i];
    if (seqA_begin < 0)
        return;
    
    const int seqA_end = query_chr[i+1];
    
    const Pos & vec2DPos = vec2D[seqA_begin];
    const int chunk = SEQLENGTH/2 - WORDSIZE/2;
    const int length = static_cast<int>(vec2DPos.size());
    
    const std::bitset<WORDSIZE*2> endA(seqA_end);
    
    for (int j = 0; j < length; ++j) {
        const Element & element = vec2DPos[j];
        const int seqB_end = element.id;
        std::bitset<WORDSIZE*2> endB(seqB_end);
        
//        if (std::abs(i*chunk-14756104) < 20 && seqA_begin == 2105376 && seqA_end == 2631720 && seqB_end == seqA_end) {
//            std::cout << "pos: " << i*chunk << ", motif hash: " << seqA_begin << std::endl;
//        }

         endB ^= endA;
        
        //////////////////////////////////////////////////////////////////////////
        // One character occupies two bits, hence less than 3
        // Increment by 2 to avoid comparing bits from different characters
        //////////////////////////////////////////////////////////////////////////
        if (endB.none())
            candidates.push_back(Candidate(i*chunk, element.idx));
        else if(endB.count() < 3) {
            int count = 0;
            for (int k = 0; k < WORDSIZE*2; k+=2) {
                if (endB[k] || endB[k+1]) count++;
            }
            
            if(count < 2)
                candidates.push_back(Candidate(i*chunk, element.idx));
        }
    }
    
    find_candidates(seqA_begin, seqA_end, i*chunk, candidates);
}

//1. get subset using 10 as the index giving subset S={3,7,4,23}
//2. ignore keys that are identical to 3 from S giving S={7,4,23}
//3. do bitwise comparison between 3 and 7,4,23 and select those that are at most 2 bit different - in this case giving 4.
//4. calculate index for 3
//5. index for 4 should be idx.
void find_candidates(const int seqA_begin, const int seqA_end, const int index, Candidates &candidates) {
    //1,2
    Pos vec2DPos = vec2D[seqA_end];
    vec2DPos.erase(std::remove_if(vec2DPos.begin(), vec2DPos.end(), [seqA_begin](const Element &e){
        return e.id == seqA_begin;
    }), vec2DPos.end());
    
    const std::bitset<WORDSIZE*2> begin_A(seqA_begin);
    for (int i = 0; i < vec2DPos.size(); ++i) {
        const Element e = vec2DPos[i];
        std::bitset<WORDSIZE*2> begin_B(e.id);
        
        //3.
        begin_B ^= begin_A;
        if(begin_B.count() < 3) {
            int count = 0;
            for (int k = 0; k < WORDSIZE*2; k+=2) {
                if (begin_B[k] || begin_B[k+1]) count++;
            }
            
            //4,5
            if(count < 2)
                candidates.push_back(Candidate(index, e.idx));
        }
    }
    
}

void find_lime_candidates(Chromosome &query_chr, Candidates &candidates) {
    
    const int iterations = static_cast<int>(query_chr.size()-1);
    for (int i = 0; i < iterations; ++i) {
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
    generateVectorIndecies(query_data_cloned, query_chr);
    find_limes_in_reverse_direction(query_chr, query_data_cloned, target_data, candidates);
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
        target_data = loadDataWithContentsOFile(targetSeqFile, targetSeqHeader);

        //////////////////////////////////////////////
        // Skips by 1 letter
        // goal: Need to minimize vec2D generation
        //////////////////////////////////////////////
        target_chr = generateLookupTableIndices(target_data);
        initializeLookupTable(target_chr);
        
        for (Files::size_type j = 0; j < g2Size; ++j) {
            const File querySeqFile(genome2[j]);
            query_data = loadDataWithContentsOFile(querySeqFile, querySeqHeader);
            
            timer.split();
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
        target_sequence = loadDataWithContentsOFile(target_file, target_header);
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


