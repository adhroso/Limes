//
//  IO.cpp
//  Limes2.0
//
//  Created by Adhroso on 7/22/16.
//  Copyright © 2016 Andi Dhroso. All rights reserved.
//

#include "IO.h"

Sequence loadDataWithContentsOFile(const File &file, Header &header) {
    std::ifstream in(file.c_str());
    if (in.is_open()) {
        // info("Loading file %s...\n", file.c_str());
        
        in.seekg (0, in.end);
        long length = in.tellg();
        //info("file size = %d\n",length);
        
        Sequence content;
        content.resize(length);
        in.seekg (0, in.beg);
        in.read(&content[0], content.size());
        in.close();
        // info("string length = %d\n",content.length());
        // info("numer of charaters extracted = %d\n",in.gcount());
        
        // get sequence header
        size_t pos = content.find("\n");
        header = content.substr(0,pos);
        
        // remove fasta header (assumes sequence has header file)
        Sequence::iterator it = std::find(content.begin(), content.end(), '\n');
        content.erase(content.begin(), it+1);
        
        // remove newline characters
        content.erase(std::remove(content.begin(), content.end(), '\n'), content.end());
        
        // remove invalid characters such as 'N'
        // content.erase(std::remove(content.begin(), content.end(), 'N'), content.end());
        
        // info("Done loading %s.\n", file.c_str());
        return content;
    } else {
        info("Check file: %s\n", file.c_str());
        exit(EXIT_FAILURE);
    }
    
    return "";
}

Sequence loadDataWithContentsOFile(const File &file) {
    std::ifstream in(file.c_str());
    if (in.is_open()) {
        // info("Loading file %s...\n", file.c_str());
        
        in.seekg (0, in.end);
        long length = in.tellg();
        // info("file size = %d\n",length);
        
        Sequence content;
        content.resize(length);
        in.seekg (0, in.beg);
        in.read(&content[0], content.size());
        in.close();
        // info("string length = %d\n",content.length());
        // info("numer of charaters extracted = %d\n",in.gcount());
        
        // remove fasta header (assumes sequence has header file)
        std::string::iterator it = std::find(content.begin(), content.end(), '\n');
        content.erase(content.begin(), it+1);
        
        // remove newline characters
        content.erase(std::remove(content.begin(), content.end(), '\n'), content.end());
        
        // remove invalid characters such as 'N'
        // content.erase(std::remove(content.begin(), content.end(), 'N'), content.end());
        
        // info("Done loading %s.\n", file.c_str());
        return content;
    } else {
        info("Check file: %s\n", file.c_str());
        exit(EXIT_FAILURE);
    }
    
    return "";
}

void load_next_batch(std::ifstream &in, std::vector<std::pair<Header, Sequence> > &sequences, const int batch_size) {
    int counter = 0;
    std::vector<std::pair<Header, Sequence> > data;
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

Files retrieve_directory_content( const Path &dir_path, const Extension &ext) {
    Path dir = dir_path;
    if (dir.back() != '/') dir.append("/");
    
    const Extension file_extension = "." + ext;
    Files files;
    files.reserve(50000);
    
    struct dirent* dp = NULL;
    struct stat st;
    
    DIR* dirp = opendir(dir.c_str());
    while (dirp) {
        if ((dp = readdir(dirp)) != NULL) {
            File s (dir+dp->d_name);
            int status = lstat(s.c_str(), &st);
            if (status != -1) {
                //match file extension with what is provided
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

void print_string_stdout(const char *s) {
    fputs(s,stdout);
    fflush(stdout);
}

void (*lime_print_string) (const char *) = &print_string_stdout;
void info(const char *fmt,...) {
    char buf[BUFSIZ];
    va_list ap;
    va_start(ap,fmt);
    vsprintf(buf,fmt,ap);
    va_end(ap);
    (*lime_print_string)(buf);
}

#pragma mark -
