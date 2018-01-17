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
 *      quality  - confidence of sequence base reads (FASTQ), keys sequence name
 *      mappings - maps sequences to the index in the backbone where they align
 * 
 * Author: Tin Kuculo
 */ 
class Data {
public:
	string backbone;
	map<int, list<string>> mappings;
    void prepare_data(string, string, string);
    
private:
    map<string, string> quality;
    map<string, string> sequence;
    string get_backbone(string);
    map<string, tuple<string, string>> get_reads(string path);
    map<int, list<string>>get_mappings(string);

public:     
    void prepare_data(string backbone_path, string reads_path, string mappings_path) {
        
        this->backbone = get_backbone(backbone_path);
        map<string, tuple<string, string>> reads = get_reads(reads_path);

        // get sequences and qualities
        for (std::map<string, tuple<string, string>>::iterator it = reads.begin(); it != reads.end(); ++it) {
            string name = it->first;
            tuple<string, string> data = it->second;
            this->sequence[name] = get<0>(data);
            this->quality[name] = get<1>(data);
        }
        
        mappings = get_mappings(mappings_path);
    }

private:
    // read layout file, return sequence object
    string get_backbone(string path) {
        list<string> lines = read_file(path);
        return clean(lines.back());
    }

    // get list of reads
    map<string, tuple<string, string>> get_reads(string path) {
        map<string, tuple<string, string>> reads;
        string name, sequence, quality;
        list<string> lines = read_file(path);
        while (!lines.empty())
        {
            name = clean(lines.front());
            name = name.substr(1, name.size()); //ignore first char '@'
            lines.pop_front();
            sequence = clean(lines.front());
            lines.pop_front();
            lines.pop_front();
            quality = clean(lines.front());
            lines.pop_front();
            reads[name] = make_tuple(sequence, quality);
        }
        return reads;
    }


    // Read mapping file, expects .paf format i
    map<int, list<string>> get_mappings(string path){
        string read_id, cigar, temp, read, q, new_read, new_quality;
        map<string, tuple<string, string> reads = this->sequence;
        int read_start, read_end, target_start, j, k;
        map<int, list<string>> mappings;
        
        ifstream m(path);
        
        // Read from file: we are only interested in a few fields, ignore all others.
        while(infile >> read_id >> temp >> read_start >> read_end >> temp >> temp >> temp >> target_start >> temp
            >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> cigar){
                
            read = reads[read_id].substr(read_start, read_end-read_start);
            q = quality[read_id].substr(read_start, read_end-read_start);
            temp = "";
            new_read = "";
            j=0;
            // Modify sequence to include insertions, deletions
            for(int i=5; i<cigar.size(); i++){
                if(isdigit(cigar[i])) 
                    temp += cigar[i];
                else{
                    k = atoi(temp);
                    switch(cigar[i]){
                        case 'M':
                            new_read += read.substr(j, k)
                            new_quality += q.substr(j, k)
                        case 'I':
                            for(int g=0; g<k; g++) new_read += tolower(read[j+g]);
                            new_quality += q.substr(j, k)
                            
                        case 'D':
                            new_read += string(k,'_');
                            new_quality += string(k,'(');
                    }
                    temp = ""
                    j += k;
                }  
            }
            // Add new sequence to dictionary
            if(mappings.count(target_start))
                mappings[target_start].push_front(make_tuple(new_read, new_quality));
            else:
                mappings[target_start] = *new list<string>(make_tuple(new_read, new_quality));
        }
        return mappings;
    }
};



