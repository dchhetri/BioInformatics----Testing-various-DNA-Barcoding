//
//  main.cpp
//  C++
//
//  Created by Dibash Chhetri on 11/11/11.
//  Copyright (c) 2011 University of Connecticut. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "utils.h"

using namespace std;

#define DEBUG_MODE false

//number of elements to remove from the loaded data for testing
const int SAMPLE_SPACE_SIZE = 50;

struct GeneElement;

typedef std::vector<GeneElement> GeneList;

//generates a sample of random sequence to be removed
std::vector<int> generateSampleSpace(const size_t SAMPLE_SIZE, size_t minRange,size_t maxRange);

//removes the sequence sample
GeneList removeSelection(GeneList& list, const std::vector<int>& sampleSpace);

//reads the genes from file
GeneList readGeneList(const std::string& filePath);

struct GeneElement{
    string phylum_name;
    string class_name;
    string order_name;
    string species_name;
    string nuclotides;
    
    friend istream& operator>>(istream& stream, GeneElement& gene){
        string line;
        getline(stream,line);
        stringstream strstream(line);
        string temp;
        int wordCount = 0;
        while(getline(strstream,line,'\t')){
            ++wordCount;
            switch(wordCount){
                case 6: gene.phylum_name = line; break;
                case 7: gene.class_name = line; break;
                case 8: gene.order_name = line; break;
                case 9: gene.species_name = line; break;
                case 18: gene.nuclotides = line; break;
            }
        }
        return stream;
    }
    string toString(const std::string& delim = "\t")const{ 
        return phylum_name + delim + class_name + delim + order_name + delim + species_name + delim + nuclotides;
    }
};


//find the closet matched gene
template<typename Heuristic>
GeneElement matchGene(const GeneList& geneList, const GeneElement& element, const Heuristic& dist);
int main (int argc, const char * argv[])
{
    using namespace std;
    srand( (unsigned) time(0) );
    const std::string filePath = "/Users/MacBoss/Desktop/Dna-database/iBOL_phase_1.25.tsv";
    
    //read data from file
    GeneList geneList = readGeneList(filePath);
    //get random sample indices
    std::vector<int> sampleSequence = generateSampleSpace(SAMPLE_SPACE_SIZE,0,geneList.size());
    //remove selected sample
    GeneList removedSelection = removeSelection(geneList,sampleSequence);
    
    int numberOfMatches = 0;
    for(int i = 0; i < removedSelection.size(); ++i){
        cout << "At i = " << i << endl;
        const GeneElement& targetGene = removedSelection[i];
        GeneElement geneMatched = matchGene(geneList,removedSelection[i],countMatches);
        bool isSameClass = targetGene.class_name == geneMatched.class_name;
#if DEBUG_MODE    
        if(!isSameClass){
        cout << setw(10) << "Matched gene = " << geneMatched.toString() << endl;
        cout << setw(10) << "Target gene = " << targetGene.toString() << endl;
        cout << setw(10) << "Matched Score = " << countMatches(geneMatched.nuclotides,targetGene.nuclotides) << " out of " << targetGene.nuclotides.size() << endl;
        cout << setw(10) << "Same class = " << boolalpha << isSameClass << endl;
        cout << "\nPress Enter....";
        cin.get();
        }
#endif
        if(isSameClass) ++numberOfMatches;
    }
    
    //calculate data
    size_t numberOfMismatches = removedSelection.size() - numberOfMatches;
    cout << setw(10) << "Number of matches = " << numberOfMatches << endl;
    cout << setw(10) << "Number of mismatches = " << numberOfMismatches << endl;
    cout << setw(10) << setprecision(3) << "Matched percentage = " << 100.0f * float(numberOfMatches) / float( removedSelection.size() ) << endl;
    cout << setw(10) << setprecision(3) << "Mismatched percentage = " << 100.0f * numberOfMismatches/removedSelection.size() << endl;
    
    return 0;
}
std::vector<int> generateSampleSpace(const size_t SAMPLE_SIZE, size_t minRange, size_t maxRange)
{
    std::vector<int> array(maxRange - minRange,0);
    for(int i = (int)minRange; i < maxRange; ++i){
        array[i] = i;
    }
    std::random_shuffle(array.begin(), array.end());
    std::vector<int> result(SAMPLE_SIZE);
    std::copy(array.begin(), array.begin() + SAMPLE_SIZE, result.begin());
    return result;
}
GeneList removeSelection(GeneList& list, const std::vector<int>& sampleSpace)
{
    GeneList removedList;
    GeneList newList;
    
    //add the elements to be removed
    for(int i = 0; i < sampleSpace.size(); ++i){
        removedList.push_back(list[sampleSpace[i]]);
    }
    //add all elements except for the ones to be removed
    for(int i = 0; i < list.size(); ++i){
        if(std::count(sampleSpace.begin(),sampleSpace.end(),i) == 0){
            newList.push_back(list[i]);
        }
    }
    
    list.swap(newList);
    return removedList;
}
GeneList readGeneList(const std::string& filePath){
    ifstream fileInput(filePath.c_str());
    GeneElement tmpGene;
    GeneList geneList;
    
    fileInput >> tmpGene; //ignore metadata;
    while(fileInput >> tmpGene){
        geneList.push_back(tmpGene);
    }
    return geneList;
}

template<typename Heuristic>
GeneElement matchGene(const GeneList& geneList, const GeneElement& element, const Heuristic& dist)
{
    long maxScore = -1;
    GeneElement maxGeneElement;
    //find the maximum score
    for(int i = 0; i < geneList.size(); ++i)
    {
        int score = dist(geneList[i].nuclotides,element.nuclotides);
        if(score > maxScore){
            maxScore = score;
            maxGeneElement = geneList[i];
        }
    }
    
    return maxGeneElement;
}
