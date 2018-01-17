#pragma once
#include <iostream>
#include <list>
#include <set>
#include <fstream>
#include <sstream>
#include <map>
#include <string.h>
#include <cstdlib>

using namespace std;
extern int k, g;

class Node;
Node *find_node(Node *, int, string);

// Edge class, contains g-length sequence and a weight
class Edge {
public:
	string *seq;
	int weight;
	Node *next;
};

// Node of k-mer graph, contains a k-length sequence and a list of edges
class Node {
public:
	string * seq;
	list<Edge> * edges;
	int index; //position in backbone

    bool operator==(const Node &n1){
        return(seq->data()==n1.seq->data() && index == n1.index);
    }
   /*
   Update current node in k-mer graph.
   If matching edge exists, add weight. Otherwise, add new edge.
   Returns: next node
   */
	void update(string seq, string quality, int index) {
		int q = 0;
		// convert quality string to int
		for (int i = 0; i < quality.size(); i++)
			q += int(quality[i])-int('!');
        
        for (list<Edge>::iterator it = this->edges->begin(); it != this->edges->end(); ++it) {
            // if matching edge found, add quality to weight
            if (seq.compare(it->seq->c_str()) == 0){
                it->weight += q;
                return;
            }
        }
            
        // else, create new edge, create new node if necessary 
        string node_seq = seq.substr(seq.size() - k, k);
        string edge_seq = seq.substr(0,g);
        
        Edge * new_edge = (Edge *)malloc(sizeof(Edge));
        new_edge->seq = new string(edge_seq);
        new_edge->weight = q;
        void * temp = find_node(this, index, node_seq);
        
        if(temp == nullptr){
            Node * new_node = (Node *)malloc(sizeof(Node));
            new_node->seq = new string(node_seq);
            new_node->index = index;
            new_node->edges = new list<Edge>();
            new_edge->next = new_node;
        }
        else
            new_edge->next = (Node *)temp;
    
        edges->push_back(*new_edge);
        
    }

    // Function to advance forward, takes edge sequence
    Node *next(string edge){
        for (list<Edge>::iterator it = this->edges->begin(); it != this->edges->end(); ++it) {
            if (edge.compare(it->seq->c_str()) == 0){
                return it->next;
            }
        }
        printf("Error! Edge not found.");
        return 0;
    }

};

// BFS from current node: find correct index and sequence
Node *find_node(Node *current, int i, string seq){
    Node *node;
    list<Node *> queue = *new list<Node *>();
    set<Node *> visited = *new set<Node *>();
    
    queue.push_back(current);
    while(!queue.empty()){
        node = queue.front();
        queue.pop_front();
        if(node->index==i && node->seq->compare(seq)==0){
            return node;
        }
        for(list<Edge>::iterator it = current->edges->begin(); it != current->edges->end(); ++it) {
            if(!visited.count(it->next))
                queue.push_back(it->next);
        }
        visited.insert(node);
    }
    return nullptr;
}

// Print whole tree from root (BFS)
void print_tree(Node *root, int maxindex){
    Node *node;
    list<Node *> queue = *new list<Node *>();
    set<Node *> visited = *new set<Node *>();
    
    queue.push_back(root);
    while(!queue.empty()){
        node = queue.front();
        queue.pop_front();
        if(node->index>=maxindex)
            continue;
        if(!visited.count(node))
            printf("\n%s[%d. %s]\n",string(node->index,' ').c_str(), node->index, node->seq->c_str());
        for(list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it) {
            if(!visited.count(node))
                printf("%s  --%s(%d)-->[%d. %s]\n",string(node->index,' ').c_str(),it->seq->c_str(), 
                                            it->weight, node->next(it->seq->data())->index, node->next(it->seq->data())->seq->c_str());
                
            if(!visited.count(it->next))
                queue.push_back(it->next);
        }
        visited.insert(node);
    }
}
    
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

// TODO: NORMAL FSCANF
// get mappings
class Mapping {
public:
	string q_name;
	int q_start;
	int q_end;
	string t_name;
	int t_start;
	int t_end;
};

list<Mapping> get_mappings(string path) {
	list<Mapping> mappings;
	list<string> lines;

	string pythonCommand = "python read_mappings.py " + path;
	system(pythonCommand.c_str()); // run python script that changes mapping data format
	lines = read_file("temp");
	//remove("temp"); //delete the temp file after reading it
    int i=0;
	while (!lines.empty()) {
		Mapping mapping;
		mapping.q_name = clean(lines.front());
		lines.pop_front();
		mapping.q_start = stoi(clean(lines.front()));
		lines.pop_front();
		mapping.q_end = stoi(clean(lines.front()));
		lines.pop_front();
		mapping.t_name = clean(lines.front());
		lines.pop_front();
		mapping.t_start = stoi(clean(lines.front()));
		lines.pop_front();
		mapping.t_end = stoi(clean(lines.front()));
		lines.pop_front();
		mappings.push_back(mapping);
        i++;
	}
	return mappings;
}

// Contains all used data
class Data {
public:
	string backbone;
	map<string, string> sequence;
	map<string, string> quality;
	map<int, list<Mapping>> index_to_mapping;

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
};

