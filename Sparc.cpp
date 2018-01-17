#include "Utils.h"
#include <string>
#include <algorithm>
#include "Data.h"

// TODO: clean up!!

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
    Node * current;
    Node * root;
	int i = k+g;
    
    // Create first node (root)
	seq = backbone.substr(0, k);
	root = (Node *)malloc(sizeof(Node));
    root->seq = new string(seq);
    root->index = 0;
    root->edges = new list<Edge>();
    
    seq = backbone.substr(k, g);
    string quality(g, '$'); // backbone quality is 
    root->update(seq, quality, g);
    current = root->next(seq);
    
    // Add remaining nodes and edges
	while (i<backbone.size()-k) {
		seq = backbone.substr(i, g);
		string quality(g, '(');
		current->update(seq, quality, i+g-k);
		i += g;
        current = current->next(seq);
	}
    
    return root;
}
// used in debug. prints node name and edges + weights
void print_node(Node * node){
    printf("%d. %s -> [ ", node->index, node->seq->c_str());
    for(list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it){
        printf("%s-%d ",it->seq->c_str(), it->weight);
    }
    printf("]\n");
}

// TODO: quality without added deletions
// gets sequence for next edge
string get_edge(string sequence){
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
        if(b == g) break;
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
// insert new read to k-mer graph
void insert(Node *current, string sequence, string quality){
    Node *root = current;
    string seq, q;
    int j, i = k;
    
    while (i<sequence.size()-g) {
		seq = get_edge(sequence.substr(i,-1));
		q = quality.substr(i, seq.size());
        seq = remove_char(seq, '_');
		current->update(seq, q, current->index+g);
		i += g;
        current = current->next(seq);
	}
    printf("Inserted: %s at index %d\n", sequence.substr(0,50).c_str(), root->index);
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
    printf("Found max path: %s ... (%d)\n", (root->seq->data()+max_path.substr(0,30)).c_str(), max_weight);
    return root->seq->data()+max_path;
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

int g, k;
int main(int argc, char* argv[])
{
	Data data;
	Node *graph, *current;
    Edge max;
	string backbone_path, reads_path, mappings_path, output_path;
    string sequence, quality, final_sequence;
    list<tuple<string, string>> mappings;
    int maxweight = 0;
    int last;

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
    
    // initialize graph: convert backbone to k-mer graph
	graph = (Node *)malloc(sizeof(Node));
    graph->seq = new string(data.backbone.substr(0, k));
    graph->index = 0;
    graph->edges = new list<Edge>();
    insert(graph, data.backbone, string(data.backbone.size(), ')')); // mid-low confidence to backbone
	printf("graph initialized\n");
    
    // go through backbone indices, add new sequence to graph if matched 
    current = graph;
    while(true){
        if(data.mappings.count(current->index)){
            mappings = data.mappings[current->index];
            for(list<tuple<string, string>>::iterator it = mappings.begin(); it != mappings.end(); ++it){
                sequence = get<0>(*it);
                quality = get<1>(*it);
                insert(current, sequence, quality);
            }
        }
        if(!current->edges->size())
            break;
        current = current->next(current->edges->front().seq->c_str());
        //printf("Next node: %d. %s\n",current->index, current->seq->c_str());
    }
    
    printf("graph constructed\n");
    
    final_sequence = get_sequence(graph);
    //print_tree(graph, 30);
    
    ofstream output;
    output.open(output_path);
    output << ">" << "consensus_output\n";
    output << final_sequence << "\n";
    output.close();
    
    printf("new sequence written to %s\n", output_path.c_str());
	return 0;
}


