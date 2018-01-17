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


// Contains all used data
class Data {
public:
	string backbone;
	map<string, string> sequence;
	map<string, string> quality;
	map<int, list<Mapping>> index_to_mapping;

    map<int, list<string>>get_mappings(string);
    
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


		list<Mapping> mappings = get_mappings(mappings_path);
		// create map where: key = index on backbone, value = list of Mappings that start at that index
		// main programs iterates over backbone and checks if there are mapped queries that start at that index.
		for (list<Mapping>::iterator it = mappings.begin(); it != mappings.end(); ++it) {
			int index = it->t_start;
			list<Mapping> temp = index_to_mapping[index];
			temp.push_back(*it);
			this->index_to_mapping[index] = temp;
		}

	}
    
    // Read mapping file, expects .paf format i
    map<int, list<string>> get_mappings(string path){
        string read_id, cigar, temp, read, new_read;
        map<string, string> reads = this->sequence;
        int read_start, read_end, target_start, j, k;
        map<int, list<string>> mappings;
        
        ifstream m(path);
        
        // Read from file: we are only interested in a few fields, ignore all others.
        while(infile >> read_id >> temp >> read_start >> read_end >> temp >> temp >> temp >> target_start >> temp
            >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> cigar){
                
            read = reads[read_id].substr(read_start, read_end-read_start);
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
                            
                        case 'I':
                            for(int g=0; g<k; g++) new_read += tolower(read[j+g]);
                            
                        case 'D':
                            new_read += string(k,'_');
                    }
                    temp = ""
                    j += k;
                }  
            }
            // Add new sequence to dictionary
            if(mappings.count(target_start))
                mappings[target_start].push_front(new_read);
            else:
                mappings[target_start] = *new list<string>(new_read);
        }
        return mappings;
    }
};



