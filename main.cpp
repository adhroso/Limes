//
//  main.cpp
//  Limes_Serial
//
//  Created by Andi Dhroso on 11/4/13.
//  Copyright (c) 2013 Andi Dhroso. All rights reserved.
//

#include "Serial.h"

#pragma mark -
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// used for Xingyan's unique gene and protein extraction
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void load_fasta_seqs(std::ifstream &in, std::vector<std::pair<std::string, std::string> > &sequences, const int batch_size=1000) {
    std::string line;
    int counter = 0;
    std::vector<std::pair<std::string, std::string> > data;
    data.reserve(batch_size);
    
    while (counter < batch_size && !in.eof()) {
        std::getline(in, line);
        const std::string header(line);
        std::string seq ("");
        while (!in.eof() && in.peek() != '>' ) {
            std::getline(in, line);
            seq.append(line);
        }
        data.push_back(std::make_pair(header, seq));
        counter++;
    }
    data.swap(sequences);
}

void load_data(std::ifstream &in, std::vector<std::string> &results, const int batch_size=1000) {
    std::string line;
    int counter = 0;
    std::vector<std::string> data;
    data.reserve(batch_size);
    
    while (counter < batch_size && !in.eof()) {
        std::getline(in, line);
        data.push_back(line);
        counter++;
    }
    data.swap(results);
}

void extract_unique_proteins() {
    std::ifstream in ("/Users/andi/Documents/Dropbox/dataset/SVNs_797_cance_genomes_all_mutations/all_proteins.fa");
    std::vector<std::pair<std::string, std::string> > sequences;
    load_fasta_seqs(in, sequences, 30000);
    std::sort(sequences.begin(), sequences.end(), [](const std::pair<std::string, std::string> &lhs, const std::pair<std::string, std::string> &rhs){
        return lhs.first < rhs.first;
    });
    std::for_each(sequences.begin(), sequences.end(), [](std::pair<std::string, std::string> &p){ p.first.erase(0,1);});
    in.close();
    
    in.open("/Users/andi/Documents/Dropbox/dataset/SVNs_797_cance_genomes_all_mutations/genes_unique.txt");
    std::vector<std::string> genes;
    load_data(in, genes, 30000);
    in.close();
    
    std::vector<std::pair<std::string, std::string> >::iterator it;
    for (int i = 0; i < genes.size(); ++i) {
        const std::string gene = genes[i];
        it = std::find_if(sequences.begin(), sequences.end(), [gene]( const std::pair<std::string, std::string> &p){
            return p.first == gene;
        });
        assert(it != sequences.end());
        std::cout << ">" << it->first << std::endl;
        std::cout << it->second << std::endl;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::string> retrieve_directory_content( const std::string &dir_path, const std::string &ext) {
    std::string dir = dir_path;
    if (dir.back() != '/')
        dir.append("/");
    
    const std::string file_extension = "." + ext;
    std::vector<std::string> files;
    files.reserve(50000);
    
    struct dirent* dp = NULL;
    struct stat st;
    
    DIR* dirp = opendir(dir.c_str());
    while (dirp) {
        if ((dp = readdir(dirp)) != NULL) {
            std::string s (dir+dp->d_name);
            int status = lstat(s.c_str(), &st);
            if (status != -1) {
                if (S_ISREG(st.st_mode) && (s.substr(s.length()-file_extension.length()) == file_extension)) {
                    files.push_back(s);
                    //info("%s\n", s.c_str());
                }
            }
        } else {
            break;
        }
    }
    
    if(dirp)
        closedir(dirp);
    else
        std::cerr << "Error, directory path is not valid: \n" << dir << std::endl;
    
    return files;
}

bool isAssembled(const std::string &file) {
    return !(file.find("dna.nonchromosomal") != file.npos);
}

bool isAssembled(std::vector<std::string> &g) {
    std::vector<std::string>::size_type i;
    for (i = 0; i < g.size(); ++i) {
        const bool assembled = isAssembled(g[i]);
        if (!assembled)
            return false;
    }
    return true;
}

int main(int argc, const char * argv[])
{
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
        std::vector<std::string> g1 = retrieve_directory_content(dir1, fileType);
        std::vector<std::string> g2 = retrieve_directory_content(dir2, fileType);
        if (g1.empty() || g2.empty()) return EXIT_FAILURE;
        
        bool g1Assembled = true;//isAssembled(g1);
        bool g2Assebmled = true;//isAssembled(g2);
        if (g1Assembled && g2Assebmled) {
            std::cout << "Both are assembled" << std::endl;
            generate_limes(g1, g2, pathToLimes, pathToProgress);
      
        } else if((g1Assembled && !g2Assebmled) || (!g1Assembled && g2Assebmled)) {
            std::cout << "one is assembled but not the other" << std::endl;
            
            if (g1Assembled && !g2Assebmled)    generate_limes(g2.front(), g1, pathToLimes, pathToProgress);
            else    generate_limes(g1.front(), g2, pathToLimes, pathToProgress);
        
        } else {
            std::vector<std::string>::size_type i, j, size_g1 = g1.size(), size_g2 = g2.size();
            std::cout << "first genome: " << dir1 << " has " << size_g1 << " number of files" << std::endl;
            std::cout << "second genome: " << dir2 << " has " << size_g2 << " number of files" << std::endl;
            std::cout << "neither is assembled" << std::endl;
            for (i = 0; i < size_g1; ++i) {
                for (j = 0; j < size_g2; ++j) {
                    generate_limes(g1.front(), g2.front(), pathToLimes, pathToProgress);
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

#pragma mark - IO
RawData loadDataWithContentsOFile(const std::string &file) {
    std::ifstream in(file.c_str());
    if (in.is_open()) {
        //        info("Loading file %s...\n", file.c_str());
        
        in.seekg (0, in.end);
        long length = in.tellg();
        //info("file size = %d\n",length);
        
        std::string content;
        content.resize(length);
        in.seekg (0, in.beg);
        in.read(&content[0], content.size());
        in.close();
        //info("string length = %d\n",content.length());
        //info("numer of charaters extracted = %d\n",in.gcount());
        
        //remove fasta header (assumes sequence has header file)
        std::string::iterator it = std::find(content.begin(), content.end(), '\n');
        content.erase(content.begin(), it+1);
        
        //remove newline characters
        content.erase(std::remove(content.begin(), content.end(), '\n'), content.end());
        
        //remove invalid characters such as 'N'
        //content.erase(std::remove(content.begin(), content.end(), 'N'), content.end());
        
        //        info("Done loading %s.\n", file.c_str());
        return content;
    } else {
        info("Check file: %s\n", file.c_str());
        exit(EXIT_FAILURE);
    }
    
    return "";
}
RawData loadDataWithContentsOFile(const std::string &file, std::string &sequenceHeader) {
    std::ifstream in(file.c_str());
    if (in.is_open()) {
        //        info("Loading file %s...\n", file.c_str());
        
        in.seekg (0, in.end);
        long length = in.tellg();
        //info("file size = %d\n",length);
        
        std::string content;
        content.resize(length);
        in.seekg (0, in.beg);
        in.read(&content[0], content.size());
        in.close();
        //info("string length = %d\n",content.length());
        //info("numer of charaters extracted = %d\n",in.gcount());
        
        //get sequence header
        size_t pos = content.find("\n");
        sequenceHeader = content.substr(0,pos);
        
        //remove fasta header (assumes sequence has header file)
        std::string::iterator it = std::find(content.begin(), content.end(), '\n');
        content.erase(content.begin(), it+1);
        
        //remove newline characters
        content.erase(std::remove(content.begin(), content.end(), '\n'), content.end());
        
        //remove invalid characters such as 'N'
        //content.erase(std::remove(content.begin(), content.end(), 'N'), content.end());
        
        //        info("Done loading %s.\n", file.c_str());
        return content;
    } else {
        info("Check file: %s\n", file.c_str());
        exit(EXIT_FAILURE);
    }
    
    return "";
}

void print_string_stdout(const char *s)
{
    fputs(s,stdout);
    fflush(stdout);
}

void (*lime_print_string) (const char *) = &print_string_stdout;
void info(const char *fmt,...)
{
    char buf[BUFSIZ];
    va_list ap;
    va_start(ap,fmt);
    vsprintf(buf,fmt,ap);
    va_end(ap);
    (*lime_print_string)(buf);
}
#pragma mark -

