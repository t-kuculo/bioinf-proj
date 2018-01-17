#include "Utils.h"


// helper function that returns file lines
list<string> read_file(string path) {
	list<string> lines;
	ifstream infile(path);
	for (string line; getline(infile, line); )
		lines.push_back(line);
	return lines;
}

// helper function that cleans up lines (remove new line at end)
string clean(string str) {
	if (str.back() == '\n')
		str = str.substr(0, str.size() - 1);
	return str;
}


/*
 * This class contains all relevant data read from files.
 * 
 * Attributes:
 *      backbone - string containing backbone sequence
 *      mappings - maps sequences and qualities to the index in the backbone where they align
 * 
 * Author: Tin Kuculo
 */ 
class Data {
    map<string, string> quality;
    map<string, string> sequence;
    
    // read layout file, return sequence object
    string get_backbone(string path) {
        list<string> lines = read_file(path);
        return clean(lines.back());
    }

    // get list of reads
    void get_reads(string path) {
        map<string, tuple<string, string>> reads;
        string name;
        list<string> lines = read_file(path);
        while (!lines.empty())
        {
            name = clean(lines.front());
            name = name.substr(1, name.size()); //ignore first char '@'
            lines.pop_front();
            sequence[name] = clean(lines.front());
            lines.pop_front();
            lines.pop_front();
            quality[name] = clean(lines.front());
            lines.pop_front();
            //printf("read: %s %s.. %s..\n", name.c_str(),sequence[name].substr(0,20).c_str(),quality[name].substr(0,20).c_str());
        }
    }

    // TODO: Do not use new_quality
    // Read mapping file, expects .paf format i
    map<int, list<tuple<string, string>>> get_mappings(string path){
        string read_id, cigar, temp, read, q, new_read, new_quality;
        int read_start, read_end, target_start, j, k;
        map<int, list<tuple<string, string>>> mappings;
       
        ifstream mfile(path);
        
        // TODO: sometimes there is extra elements in line: get proper tag
        // Read from file: we are only interested in a few fields, ignore all others.
        while(mfile >> read_id >> temp >> read_start >> read_end >> temp >> temp >> temp >> target_start >> temp){
            while (temp.substr(0,2) != "cg") mfile >> temp;
            cigar = temp.substr(5,-1);
            //printf("read: %s %d %d %d %s..\n", read_id.c_str(), read_start, read_end, target_start, cigar.substr(0,30).c_str());
            read = sequence[read_id].substr(read_start, read_end-read_start);
            q = quality[read_id].substr(read_start, read_end-read_start);
            //printf("cigar: %s..\nsequence:     %s..\nquality:      %s..\n", cigar.substr(0,15).c_str(), read.substr(0,50).c_str(), q.substr(0,50).c_str());
            temp = "";
            new_read = "";
            new_quality = "";
            j=0;
            // Modify sequence to include insertions, deletions
            for(int i=0; i<cigar.size(); i++){
                if(isdigit(cigar[i])) 
                    temp += cigar[i];
                else{
                    k = atoi(temp.c_str());
                    switch(cigar[i]){
                        case 'M':
                            new_read += read.substr(j, k);
                            new_quality += q.substr(j, k);
                            j += k;
                            break;
                            
                        case 'I':
                            for(int g=0; g<k; g++) new_read += tolower(read[j+g]);
                            new_quality += q.substr(j, k);
                            j += k;
                            break;
                            
                        case 'D':
                            new_read += string(k,'_');
                            new_quality += string(k,'(');
                    }
                temp = "";
                }  
            }
            
            // Add new sequence to dictionary
            if(!mappings.count(target_start))
                mappings[target_start] = *new list<tuple<string, string>>();
            mappings[target_start].push_front(make_tuple(new_read, new_quality));
            //printf("new sequence: %s\nnew quality:  %s\n\n", new_read.substr(0,50).c_str(), new_quality.substr(0,50).c_str());
        }
        return mappings;
    }

public:
	string backbone;
	map<int, list<tuple<string, string>>> mappings;
 
    void prepare_data(string backbone_path, string reads_path, string mappings_path) {
        
        this->backbone = get_backbone(backbone_path);
        get_reads(reads_path);
        mappings = get_mappings(mappings_path);
    }

};
