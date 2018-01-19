#include "Utils.h"
#include <string>
#include <algorithm>
#include "Data.h"

// TODO: clean up!!

Node * create_node(string *seq, int index, list<Edge> *edges){
    Node * node = (Node *)malloc(sizeof(Node));
    node->seq = seq;
    node->index = index;
    node->edges = edges;
    return node;
}
Edge create_edge(string *seq, int weight, Node *next){
    Edge *edge = (Edge *)malloc(sizeof(Edge));
    edge->seq = seq;
    edge->weight = weight;
    if(next)
        edge->next = next;
    return *edge;
}

void test_graph(Node * current, int limit) {
    int i=0;
	while (i<=limit && current->edges->size()) {
		printf("[%d. %s]---%s(%d)-->", current->index, current->seq->c_str(), 
        current->edges->front().seq->c_str(), current->edges->front().weight);
		current = current->next(current->edges->front().seq->c_str());
        i++;
	}
    if(i>=limit)
        printf(" .. (cont.)\n");
}


// used in debug. prints node name and edges + weights
void print_node(Node * node){
    printf("%d. %s -> [ ", node->index, node->seq->c_str());
    for(list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it){
        printf("%s-%d ",it->seq->c_str(), it->weight);
    }
    printf("]\n");
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


// TODO: move to utils
string remove_char(string s, char a){
    string new_s;
    for(int i=0; i<s.size(); i++)
        if (s[i]!= a) new_s += s[i];
    return new_s;
}


void insert(Node *current, string sequence, string quality){
    Node *root = current;
    string seq,new_edge, q;
    bool found=false;
    int i=0;
    
    printf("Inserting: %s at index %d \n", sequence.substr(0,10).c_str(), current->index+g);
    
    while (i<sequence.size()-g) {
		seq = get_edge(sequence.substr(i,-1), g);
		q = quality.substr(i, seq.size());
        int t = i;
        i += seq.size();
        new_edge = remove_char(seq, '_');
        //printf("New edge: %s (from: %s..)\n", seq.c_str(), sequence.substr(t,10).c_str());
		current->update(new_edge, q, current->index+g);
        current = current->next(new_edge);
	}
    //print_tree(graph, root->index+3*g);
    //printf("Inserted:  %s at index %d\n", sequence.substr(0,50).c_str(), root->index);
}

 /*
  * Initializes graph: convert backbone sequence to k-mer graph.
  * 
  * input:   backbone string
  * returns: pointer to root node
  * 
  * author: Ana Brassard
 */
Node *init_graph(string backbone) {
	string seq;
    Node *root, *first;
    Edge start;
	int i = k+g;
    
    // graph starts with empty node with edge that points to first k of backbone:
    //  []--AA-->[AA]   i.e. [root]--start-->[first]
    first = create_node(new string(backbone.substr(0, k)), 0, new list<Edge>());
    start = create_edge(new string(backbone.substr(0, k)), 0, first);
    list<Edge> *temp = new list<Edge>();
    temp->push_front(start);
    root = create_node(new string("*"), -g, temp);
    graph = root;
    //convert backbone to k-mer graph
    insert(first, backbone.substr(k,-1), string(backbone.size()-k, ')')); // mid-low confidence to backbone
    
    return root;
}

string get_sequence(Node *root){
    // list of paths: (path, weight, next node)
    list<tuple<string, int, Node *>> paths, new_paths;
    list<Node *> next;
    string path, max_path;;
    Node * current;
    int weight, last;
    int max_weight = 0;
    
    // get index of last node
    current = root;
    while(current->edges->size())
        current = current->next(current->edges->front().seq->data());
    last = current->index;
    
    // init
    current = root;
    for(list<Edge>::iterator it = current->edges->begin(); it != current->edges->end(); ++it){
        paths.push_back(make_tuple(it->seq->data(), it->weight, it->next));
        next.push_back(it->next);
    }
    
    while(!next.empty()){
        current = next.front();
        next.pop_front();
        new_paths.clear();
        for(list<tuple<string, int, Node *>>::iterator p = paths.begin(); p != paths.end(); ++p){
            if(!(*get<2>(*p)==*current)){
                new_paths.push_front(*p);
                continue;
            }
            
            // if node has no further edges and weight > current max and path goes to last node, add string path to found paths
            if(current->index == last && get<1>(*p) > max_weight){
                max_weight = get<1>(*p);
                max_path = get<0>(*p);
                continue;
            }
            
            // otherwise, update path with following edges
            for(list<Edge>::iterator it = current->edges->begin(); it != current->edges->end(); ++it){
                path = get<0>(*p) + it->seq->data();
                weight = get<1>(*p) + it->weight;
                new_paths.push_back(make_tuple(path, weight, it->next));
                if(find(next.begin(), next.end(),it->next) == next.end())
                    next.push_back(it->next);
                
            }
        }
        paths = new_paths;
    }
    printf("Found max path: %s ... (%d)\n", (max_path.substr(0,30)).c_str(), max_weight);
    return max_path;
}

// return backbone node at index i
Node *get_backbone_node(int i){
    Node * current = graph;
    while(current->index != i)
        current = current->edges->front().next;
    return current;
}

// iterate over data, update graph, find best path, write to file
/*
 * Main program of Sparc algorithm. 
 * Creates backbone k-mer graph, updates graph for each matched read, searches and finds heaviest path.
 * New path is consensus.
 * 
 * Needs input:
 *   -g : g (int)
 *   -k : k (int)
 *   -b : backbone file (.fasta)
 *   -r : reads file (.fastq)
 *   -m : mapping file (.paf)
 *   -o : output file
 * 
 * Author: Ana Brassard
*/
//global variables
int g, k; 
Node *graph;
int main(int argc, char* argv[])
{
	Data data;
	Node *current;
    Edge edge;
	string backbone_path, reads_path, mappings_path, output_path;
    string sequence, quality, final_sequence, new_edge;
    list<tuple<string, string>> mappings;
    int prev, index, maxweight = 0;;

	// get input parameters
	for (int i = 1; i < argc; ++i)
	{
        // TODO: add help option
        // TODO: check file formats
		if (strcmp(argv[i], "-g")==0) {
			i++;
			g = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "-k")==0) {
			i++;
			k = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "-b")==0) {
			i++;
			backbone_path = argv[i];
			continue;
		}
		if (strcmp(argv[i], "-r")==0) {
			i++;
			reads_path = argv[i];
			continue;
		}
		if (strcmp(argv[i], "-m")==0) {
			i++;
			mappings_path = argv[i];
			continue;
		}
        if (strcmp(argv[i], "-o")==0) {
			i++;
			output_path = argv[i];
			continue;
		}
	}
	printf("k=%d, g=%d\n",k,g);
	data.prepare_data(backbone_path, reads_path, mappings_path);
	printf("data ready\n");
    
    // initialize graph: 
	graph = init_graph(data.backbone);
	printf("graph initialized\n");
    //print_tree(graph, 30);
    
    int strt,cnt = 0;
    //go through mappings, insert to graph
    for(map<int, list<tuple<string, string>>>::iterator m_it=data.mappings.begin(); m_it!= data.mappings.end(); ++m_it){
        index = m_it->first;
        if(index>=12) return 0;
        mappings = m_it->second;
        // find node indices before and after matched index (e.g. if g=3, for index 2 get 0, 3)
        for(int i=0; i*g<=index; i++)
            prev = i*g;
        if(index==0) prev = -g;  // special case: adding to root
        //printf("Adding sequence %d to index %d, starting from %d\n", index, prev+g,prev);
        for(list<tuple<string, string>>::iterator it = mappings.begin(); it != mappings.end(); ++it){
            sequence = get<0>(*it);
            quality = get<1>(*it);
            //trim sequence start to align with node if not adding to root
            if(index>0){
                sequence = sequence.substr(prev+g-index, -1);
                quality = quality.substr(prev+g-index, -1);
            }
           
            // while sequence start is deletion of length k, skip to next node index
            while(sequence.substr(0,k)==string(k,'_')){
                sequence = sequence.substr(g, -1);
                quality = sequence.substr(g, -1);
                prev += g;
            }
            //printf("mapping:\n    %s\n    %s\n", sequence.substr(0,25).c_str(), data.backbone.substr(prev+g,25).c_str());
            // find first node match
            current = find_node(graph, prev+g, sequence.substr(0,k)); //TODO: why doesn't it find correct node?
            print_tree(graph, prev+3*g);
            if(current){
                // if node is found, insert rest of sequence
                printf("\nMatching root %s found at %d\n", current->seq->c_str(), prev+g);
                insert(current, sequence.substr(k,-1), quality.substr(k,-1));
            }
            else{
                // if node is not found, create new edge from previous backbone node, then insert rest
                current = get_backbone_node(prev);
                new_edge = get_edge(sequence, k);
                current->update(new_edge, quality.substr(0, new_edge.size()), prev+g);
                current = current->next(new_edge);
                printf("\nRoot of %s not found at %d, created new node: %s \n", 
                                                sequence.substr(0, 15).c_str(), prev+g, current->seq->c_str());
                insert(current, sequence.substr(new_edge.size(), -1), quality.substr(new_edge.size(), -1));
            }
            cnt++;
        }
    }
    
    /*
    // go through backbone indices, add new sequence to graph if matched 
    current = graph;
    do{
        if(data.mappings.count(current->index)){
            //printf("At index: %d, %d mappings\n",current->index, data.mappings[current->index].size());
            mappings = data.mappings[current->index];
            for(list<tuple<string, string>>::iterator it = mappings.begin(); it != mappings.end(); ++it){
                sequence = get<0>(*it);
                quality = get<1>(*it);
                insert(get_backbone_node(prev), sequence, quality);
            }
        }
        current = current->next(current->edges->front().seq->c_str()); //edge to next backbone node is always at head of edge list
        //printf("Next node: [%d. %s]\n", current->index, current->seq->substr(0,25).c_str());
    }while(current->edges->size());
    */
    printf("graph constructed, %d sequences inserted\n", cnt);
    print_tree(graph, 30);
        
    final_sequence = get_sequence(graph);
    
    ofstream output;
    output.open(output_path);
    output << ">" << "consensus_output\n";
    output << final_sequence << "\n\n";
    output.close();
    
    printf("new sequence written to %s\n", output_path.c_str());
	return 0;
}


