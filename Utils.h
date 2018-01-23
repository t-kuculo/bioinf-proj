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

using namespace std;

// remove a character from string
string remove_char(string s, char a){
    string new_s;
    for(int i=0; i<s.size(); i++)
        if (s[i]!= a) new_s += s[i];
    return new_s;
}

// convert quality string to float
float get_quality(string quality){
    float w = 0;
    for (int i = 0; i < quality.size(); i++)
        w += int(quality[i])-int('!');
    w /= quality.size();
    return w;
}

// print progress
void print_progress(int i, int j){
    cout<<"[";
    int pos = 30*((float)i)/j;
    for(int j=0; j<30; ++j){
        if(j<pos) cout<<"#";
        else cout <<" ";
    }
    cout<<"]"<<(int)(((float)i/j)*100)<<"%\r";
    cout.flush();
}

// gets q-length sequence for next edge
string get_edge(string sequence, int q){
    int b = 0;
    string edge;
    for(int i=0; i<sequence.size(); i++){
        if(isupper(sequence[i])){
            b++;
            edge += sequence[i];
        }
        else if(islower(sequence[i])){
            edge += toupper(sequence[i]);
        }
        else if(sequence.at(i) == '_'){
            edge += sequence[i];
        }
        if(b == q) break;
    }
    return edge;
}


