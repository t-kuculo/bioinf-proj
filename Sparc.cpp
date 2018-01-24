/* 
 * This file contains the main program: 
 *      1. construct k-mer graph
 *      2. add reads
 *      3. find heaviest path
 * 
 *  author: Ana Brassard
 */
#include "Data.h"
#include <ctime>

// global variables:
int g, k; 
Node *graph;


// Initializes graph: convert backbone sequence to k-mer graph.
void init_graph(string backbone) {
    Node *first;
    Edge first_edge;
    
    // graph starts with empty node with edge that points to first k of backbone:
    //  []--AA-->[AA]   i.e. [graph]--first_edge-->[first]
    graph = CreateNode("*", -g);
    first = CreateNode(backbone.substr(0, k), 0);
    first_edge = CreateEdge(backbone.substr(0, k), getQuality(string(k, ')')), first);
    graph->edges->push_front(first_edge);
    
    // add rest of backbone
    Node *node, *current;
    Edge edge;
    int i=k, index=g;
    string seq;
    current = first;
    while(i<backbone.size()){
        // get next edge
        seq = backbone.substr(i,g);
        edge = CreateEdge(seq, getQuality(string(g, ')')), nullptr);
        
        // get next node from edge
        if(seq.size()>=k)
            seq = seq.substr(g-k, k);
        node = CreateNode(seq, index);
        
        // connect: current --> edge --> node
        edge.next = node;
        current->edges->push_front(edge);
        
        index += g;
        i+=g;
        current = node;
    } 
    // add empty node to end of backbone
    node = CreateNode("*", index);
    edge = CreateEdge("",0,node);
    current->edges->push_back(edge);
}


// Traverses graph and finds heaviest sequence (string).
string get_best_sequence(){
    
    // get index of last node
    Node * current = graph;
    int last;
    while(current->edges->size())
        current = current->edges->front().next;
    last = current->index;

    // initialize path lists
    list<tuple<string, float, Node *>> paths, new_paths;    // list of paths: [(path, weight, next node)]
    current = graph;
    for(list<Edge>::iterator it = current->edges->begin(); it != current->edges->end(); ++it)
        paths.push_back(make_tuple(it->seq->data(), it->weight, it->next));
    
    // traverse graph from root to last (not including) 
    double weight;
    string path;
    while(true){
        new_paths.clear();
        for(list<tuple<string, float, Node *>>::iterator p = paths.begin(); p != paths.end(); ++p){
            for(list<Edge>::iterator it = get<2>(*p)->edges->begin(); it != get<2>(*p)->edges->end(); ++it){
                path = get<0>(*p) + it->seq->data();
                weight = get<1>(*p) + it->weight;
                // check if new path to node is heavier than last found path to that node
                bool found = false;
                for(list<tuple<string, float, Node *>>::iterator n_p = new_paths.begin(); n_p != new_paths.end(); n_p++){
                    if(get<2>(*n_p)->seq->compare(it->next->seq->data()) == 0){
                        found = true;    
                        if(weight > get<1>(*n_p))
                            *n_p = make_tuple(path, weight, it->next);
                    }
                }
                if(!found)
                    new_paths.push_back(make_tuple(path, weight, it->next));
            }
        }
        if(new_paths.empty())
            break;
        printProgress(get<2>(new_paths.front())->index, last);
        paths = new_paths;
    }
    
    //printf("\nFound max paths:\n");
    // Iterate over found paths to last node, get heaviest
    string max_path;
    double max_weight = 0;
    if(paths.size()>1){
        for(list<tuple<string, float, Node *>>::iterator p = paths.begin(); p != paths.end(); ++p){
            if(get<1>(*p) > max_weight){
                max_weight = get<1>(*p);
                max_path = get<0>(*p);
    }}}
    else{
        max_path = get<0>(paths.front());
        max_weight = get<1>(paths.front());
    }
   // printf("\nFound max path: %s ... (%f)\n", (max_path.substr(0,30)).c_str(), max_weight);
    return max_path;
}


// iterate from beginning of backbone to last index of mappings
void add_sequences(map<int, list<tuple<string, string>>> mappings){
    // add matched reads to graph    
    list<Node *> current_nodes, previous_nodes;
    set<Node *> visited;
        
    previous_nodes.push_back(graph);
    current_nodes.push_back(graph->edges->front().next);
    
    for(int i=0; i<=mappings.rbegin()->first; i++){
        // if i passed current node, advance forward: get set of nodes at next index
        if(i > current_nodes.front()->index){
            visited.clear();
            previous_nodes = current_nodes;
            current_nodes.clear();
            for(list<Node *>::iterator n = previous_nodes.begin(); n!= previous_nodes.end(); ++n){
                if(!(*n)->edges->size()) continue;
                for(list<Edge>::iterator e = (*n)->edges->begin(); e!= (*n)->edges->end(); ++e){
                    if(visited.count(e->next)==0){
                        current_nodes.push_back(e->next);
                        visited.insert(e->next);
            }}}
        }
        Node *current, *previous, *new_current;

        Edge edge;
        string sequence, quality, new_edge;
        int prev, index, maxweight = 0;
    
        // if no mappings at current index, continue
        if(!mappings.count(i)) continue;
        
        // else: add mapped sequences to graph
        for(list<tuple<string, string>>::iterator it = mappings[i].begin(); it != mappings[i].end(); ++it){
            sequence = get<0>(*it);
            quality = get<1>(*it);
            
            //trim sequence start to align with current node
            if(i%g != 0){
                sequence = sequence.substr(i%g, -1);
                quality = quality.substr(i%g, -1);
            }
            
            // look for matching node in current nodes
            bool found = false;
            if(i==0)
                new_edge = getEdge(sequence, k); // special case: adding to root
            else
                new_edge = getEdge(sequence, g);
            for(list<Node *>::iterator n = current_nodes.begin(); n != current_nodes.end(); ++n){
                if((*n)->seq->compare(new_edge.substr(new_edge.size()-k, k))==0){
                    found = true;
                    current = *n;
                    break;
                }
            }
            
            // if sequence root does not match, create new edge from previous backbone node and add new node
            if(!found){
                current = CreateNode(new_edge.substr(new_edge.size()-k, k).data(), current_nodes.front()->index);
                edge = CreateEdge(removeChar(new_edge, '_'), getQuality(quality.substr(0, new_edge.size())), current);
                previous_nodes.front()->edges->push_back(edge);
                current_nodes.push_back(current);
            }
            
            // insert rest of sequence
            sequence = sequence.substr(new_edge.size(), -1);
            quality = quality.substr(new_edge.size(), -1);
            Insert(current, sequence, quality, previous_nodes, current_nodes); //see: Node.h
        }
        printProgress(i, mappings.rbegin()->first);
         
    }
}


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
*/
int main(int argc, char* argv[])
{
    clock_t begin = clock();
	// get input parameters
    string backbone_path, reads_path, mappings_path, output_path;
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
    
    Data data;
	data.prepare_data(backbone_path, reads_path, mappings_path);
	printf("data ready\n");
    
	init_graph(data.backbone);
    data.backbone = "";
	printf("graph initialized\n");
    
    printf("adding sequences...\n");
    add_sequences(data.mappings);    
    data.mappings.clear();
    //print_tree(graph, 27);
    
    printf("\nfinding best path...\n");
    string best_sequence = get_best_sequence();
    
    // write to file
    ofstream output;
    output.open(output_path);
    output << ">" << "consensus_output\n";
    output << best_sequence << "\n\n";
    output.close();
    
    clock_t end = clock();        
    printf("\nnew sequence written to %s\n", output_path.c_str());
    printf("\nCPU time elapsed: %fs\n", double(end-begin)/CLOCKS_PER_SEC);
    long long RAM = getRAM();
    printf("total RAM used: %lldB (~%fGB)\n", RAM, (double)(RAM/1024/1024)/1024);
    
	return 0;
}


// idea: when inserting index > i, delete nodes before i, only keep heaviest paths