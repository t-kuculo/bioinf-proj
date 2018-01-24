#pragma once
#include <iostream>
#include <list>
#include <set>
#include <fstream>
#include <sstream>
#include <map>
#include <string.h>
#include <cstdlib>
#include <string>
#include <algorithm>
#include "sys/types.h"
#include "sys/sysinfo.h"

using namespace std;


long long getRAM(){
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    
    long long totalVirtualMem = memInfo.totalram;
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;
    
    long long totalPhysMem = memInfo.totalram;
    totalPhysMem *= memInfo.mem_unit;
    return totalPhysMem ;
}

// remove a character from string
string removeChar(string s, char a){
    string new_string;
    for(int i=0; i<s.size(); i++)
        if (s[i]!= a) new_string += s[i];
    return new_string;
}

// convert quality string to float
float getQuality(string quality){
    double weight = 0;
    for (int i = 0; i < quality.size(); i++)
        weight += int(quality[i])-int('!');
    weight /= quality.size();
    return weight;
}

// print progress
void printProgress(int i, int j){
    cout<<"[";
    int position = 30*((float)i)/j;
    for(int j=0; j<30; ++j){
        if(j<position) cout<<"#";
        else cout <<" ";
    }
    cout<<"]"<<(int)(((float)i/j)*100)<<"%   ("<<i<<"/"<<j<<")\r";
    cout.flush();
}

// gets an appropriate length(q+inserted bases) sequence for next edge
string getEdge(string sequence, int q){
    int edge_counter = 0;
    string edge;
    for(int i=0; i<sequence.size(); i++){
        if(islower(sequence[i]))
            edge += toupper(sequence[i]);
        else{
            edge += sequence[i];
            edge_counter++;
        }
        if(edge_counter == q) break;
    }
    return edge;
}


