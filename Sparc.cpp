#include "Data.h"


 /*
  * Initializes graph: convert backbone sequence to k-mer graph.
  * 
  * input:   backbone string
  * 
  * author: Ana Brassard
 */
void init_graph(string backbone) {
    Node *first;
    Edge first_edge;
    
    // graph starts with empty node with edge that points to first k of backbone:
    //  []--AA-->[AA]   i.e. [root]--start-->[first]
    graph = create_node("*", -g,  new list<Edge>());
    first = create_node(backbone.substr(0, k), 0, new list<Edge>());
    first_edge = create_edge(backbone.substr(0, k), get_quality(string(k, ')')), first);
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
        edge = create_edge(seq, get_quality(string(g, ')')), nullptr);
        
        // get next node from edge
        if(seq.size()>=k)
            seq = seq.substr(g-k, k);
        node = create_node(seq, index, new list<Edge>());
        
        // connect: current --> edge --> node
        edge.next = node;
        current->edges->push_front(edge);
        
        index += g;
        i+=g;
        current = node;
    } 
    // add empty node to end of backbone
    node = create_node("*", index, new list<Edge>());
    edge = create_edge("",0,node);
    current->edges->push_back(edge);
}

string get_sequence(Node *root){
    // list of paths: (path, weight, next node)
    list<tuple<string, float, Node *>> paths, new_paths;
    list<Node *> next;
    string path, max_path;;
    Node * current;
    int weight, last;
    float max_weight = 0;
    
    // get index of last node
    current = root;
    while(current->edges->size())
        current = current->edges->front().next;
    last = current->index;
    printf("Last index: %d\n", last);
    
    // init
    current = root;
    for(list<Edge>::iterator it = current->edges->begin(); it != current->edges->end(); ++it){
        paths.push_back(make_tuple(it->seq->data(), it->weight, it->next));
        next.push_back(it->next);
    }
    
    while(true){
        new_paths.clear();
        for(list<tuple<string, float, Node *>>::iterator p = paths.begin(); p != paths.end(); ++p){
            // if weight > current max and path goes to last node, add string path to found paths
            if(get<1>(*p) > max_weight){
                max_weight = get<1>(*p);
                max_path = get<0>(*p);
                //printf("New max: %f\n", max_weight);
            }
            
            // otherwise, update path with following edges
            if(!get<2>(*p)->edges) continue;
            for(list<Edge>::iterator it = get<2>(*p)->edges->begin(); it != get<2>(*p)->edges->end(); ++it){
                path = get<0>(*p) + it->seq->data();
                weight = get<1>(*p) + it->weight;
                // Add new path to node only if heaviest path to it
                bool found = false;
                for(list<tuple<string, float, Node *>>::iterator p_ = new_paths.begin(); p_ != new_paths.end(); ++p_){
                    if(!get<2>(*p_)->seq->compare(it->next->seq->data())){
                        found = true;
                        if(weight>get<1>(*p_)){
                            *p_ = make_tuple(path, weight, it->next);
                            break;
                        }
                    }
                }
                if(!found)
                    new_paths.push_back(make_tuple(path, weight, it->next));
            }
        }
        if(new_paths.empty()) break;
        cout<<"[";
        int pos = 30*(float)(get<2>(new_paths.front())->index)/last;
        for(int j=0; j<30; ++j){
            if(j<pos) cout<<"#";
            else cout <<" ";
        }
        cout<<"]"<< (int) (pos/30.0*100)<<"%\r";
        cout.flush();
        paths = new_paths;
    }
    printf("\nFound max path: %s ... (%f)\n", (max_path.substr(0,30)).c_str(), max_weight);
    return max_path;
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
 * 
 * Author: Ana Brassard
*/
int g, k; 
Node *graph;
int main(int argc, char* argv[])
{
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
    
    // get data
    Data data;
	data.prepare_data(backbone_path, reads_path, mappings_path);
	printf("data ready\n");
    
    // initialize graph
	init_graph(data.backbone);
	printf("graph initialized\n");

    // add matched reads to graph    
    list<Node *> current_nodes, previous_nodes;
    set<Node *> visited;
        
    previous_nodes.push_back(graph);
    current_nodes.push_back(graph->edges->front().next);
    printf("adding sequences...\n");
    
    // iterate from beginning of backbone to last index of mappings
    for(int i=0; i<=data.mappings.rbegin()->first; i++){
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
        if(!data.mappings.count(i)) continue;
        
        // else: add mapped sequences to graph
        for(list<tuple<string, string>>::iterator it = data.mappings[i].begin(); it != data.mappings[i].end(); ++it){
            sequence = get<0>(*it);
            quality = get<1>(*it);
            
            //trim sequence start to align with current node
            if(i%g != 0){
                sequence = sequence.substr(i%g, -1);
                quality = quality.substr(i%g, -1);
            }
            
            // look for matching node in current nodes
            bool found = false;
            new_edge = get_edge(sequence, g);
            for(list<Node *>::iterator n = current_nodes.begin(); n != current_nodes.end(); ++n){
                if((*n)->seq->compare(new_edge.substr(new_edge.size()-k, k))==0){
                    found = true;
                    current = *n;
                    break;
                }
            }
            
            // if sequence root does not match, create new edge from previous backbone node and add new node
            if(!found){
                current = create_node(new_edge.substr(new_edge.size()-k, k).data(), current_nodes.front()->index, new list<Edge>());
                edge = create_edge(remove_char(new_edge, '_'), get_quality(quality.substr(0, new_edge.size())), current);
                previous_nodes.front()->edges->push_back(edge);
                current_nodes.push_back(current);
            }
            
            // insert rest of sequence
            sequence = sequence.substr(new_edge.size(), -1);
            quality = quality.substr(new_edge.size(), -1);
            insert(current, sequence, quality, previous_nodes, current_nodes);
        }
        cout<<"[";
        int pos = 30*((float)(i))/data.mappings.rbegin()->first;
        for(int j=0; j<30; ++j){
            if(j<pos) cout<<"#";
            else cout <<" ";
        }
        cout<<"]"<<(int)(((float)i)/(data.mappings.rbegin()->first)*100)<<"%\r";
        cout.flush();
         
    }
    cout<<endl;
    
    printf("graph constructed\n");
    print_tree(graph, 27);
        
    printf("finding best path...\n");
    string final_sequence = get_sequence(graph);
    
    ofstream output;
    output.open(output_path);
    output << ">" << "consensus_output\n";
    output << final_sequence << "\n\n";
    output.close();
    
    printf("new sequence written to %s\n", output_path.c_str());
	return 0;
}


