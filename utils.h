//
//  utils.h
//  C++
//
//  Created by Dibash Chhetri on 12/1/11.
//  Copyright (c) 2011 University of Connecticut. All rights reserved.
//

#ifndef C___utils_h
#define C___utils_h

#include <iostream>
template<typename ForwardIterator>
void print(ForwardIterator begin, ForwardIterator end){
    while( begin != end){
        std::cout << *begin << " ";
        ++begin;
    }
}

//returns the nuber of matches
int countMatches(const std::string& lhs, const std::string& rhs){
    int count = 0;
    for(int i = 0; i < lhs.size(); ++i){
        if(lhs[i] == rhs[i]) ++count;
    }
    return count;
}

int random(int min, int max){
    return rand() * (max-min) + min;
}

void skipInput(std::istream& stream, size_t n, const char delim = '\t'){
    char ch;
    while(n > 0 && stream.get(ch))
    {
        if(ch == delim){
            --n;
        }
        if(stream && stream.peek() == delim){
            std::cout << "tab found\n";
        }
    }
}
#endif
