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
Edge create_edge(string *seq, float weight, Node *next){
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
		printf("[%d. %s]---%s(%f)-->", current->index, current->seq->c_str(), 
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
        printf("%s-%f ",it->seq->c_str(), it->weight);
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

// convert quality string to float
float get_quality(string quality){
    float w = 0;
    for (int i = 0; i < quality.size(); i++)
        w += int(quality[i])-int('!');
    w /= quality.size();
    return w;
}

/*
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
*/

void insert(Node *current, string sequence, string quality, list<Node *> previous_nodes, list<Node *> current_nodes ){
    Node *node *next_backbone;
    list<Node *> queue;// = *new list<Node *>();
    set<Node *> visited;// = *new set<Node *>();
    map<tuple<int, string>, Node *> nodes;
    map<int, Edge> edges;
    string seq, q, new_edge, new_node_seq;
    int index = current->index;
    int i=0;
    
    nodes[make_tuple(index, current->seq->data())] = current;
    
    
    // create map of all nodes to add (new nodes)
    while(i<sequence.size()-g){
        seq = get_edge(sequence.substr(i,-1), g);
		q = quality.substr(i, seq.size());
        
        new_edge = remove_char(seq, '_');
        edges[index] = create_edge(new string(new_edge), get_quality(q), nullptr);
        
        new_node_seq = seq.substr(seq.size()-k, k);
        node = new_node(new_node_seq, index+g);
        nodes[make_tuple(node->index, new_node_seq)] = node;
       
        index += g;
        i+=seq.size();
    }
    
    //printf("Nodes and edges to add ready\n");
    while(true){
        //printf("current index: %d\n", current_nodes.front()->index);
        visited.clear();
        previous_nodes = current_nodes;
        current_nodes.clear();
        for(list<Node *>::iterator n = previous_nodes.begin(); n!= previous_nodes.end(); ++n){
            //printf("node: %s\n", (*n)->seq->c_str());
            if(nodes.count(make_tuple((*n)->index, (*n)->seq->data())))
                nodes[make_tuple((*n)->index, (*n)->seq->data())] = *n;
            
            if(!(*n)->edges->size()) continue;
            for(list<Edge>::iterator e = (*n)->edges->begin(); e!= (*n)->edges->end(); ++e){
                //printf("edge: %s\n", (*e).seq->c_str());
                if(!visited.count(e->next)){
                    current_nodes.push_back(e->next);
                    visited.insert(e->next);
                }
            }
        }
        if(visited.empty()) break;
    }
    /*

    
    // DFS tree, if some node already exists, change pointer to existing node
    // skip nodes with higher index than possible with current sequence
    queue.push_front(graph);
    while(!queue.empty()){
        node = queue.front();
        queue.pop_front();
        if(nodes.count(make_tuple(node->index, node->seq->data()))){
            nodes[make_tuple(node->index, node->seq->data())] = node;
            //printf("Found node: %d. %s\n", node->index, node->seq->data());
        }
        for(list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it) {
            if(!visited.count(it->next) && it->next->index < current->index+sequence.size())
                queue.push_front(it->next);
        }
        visited.insert(node);
    }
    */
    // convert node map to node array
    Node *node_list[edges.size()];
    i = 0;
    for(map<tuple<int, string>, Node *>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
        node_list[i++] = it->second;

    //printf("\nInserting %d elements at index %d from %s\n", (int)edges.size(), node_list[0]->index, sequence.substr(0,25).c_str());
    //printf("Starting from node: [%d. %s]...\n", node_list[0]->index, node_list[0]->seq->c_str());
    
    Edge edge;
    bool found;
    //For found and created nodes...
    for(i=0; i<edges.size(); i++){
         node = node_list[i];
         //printf("Current node: %d. %s (%d edges)\n", node->index, node->seq->c_str(), (int)node->edges->size());
         edge = edges[node->index];
         found=false;
         for (list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it) {
            // if matching edge found, add quality to weight
            if (edge.seq->compare(it->seq->data()) == 0){
                //printf("Found edge: %s (%f)\n", edge.seq->c_str(), edge.weight);
                it->weight += edge.weight;
                found=true;
                break;
            }
        }
        // else, create new edge, connect to previous and next node
        if(!found){
            edge.next = node_list[i+1];
            node->edges->push_back(edge);
            //printf("Added edge %s to node %d, next = [%d. %s]\n", edge.seq->c_str(), node->index,
                                                              //  node->edges->back().next->index,node->edges->back().next->seq->c_str());
        }
    }
    // get next backbone node
    index = current->index;
    while(current->edges->front()->index <= index+sequence.size())
        current = current->edges->front().next;
    new_edge = create_edge("",0,current);
    node_list[i]->edges->push_back(new_edge);
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
    start = create_edge(new string(backbone.substr(0, k)), get_quality(string(k, ')')), first);
    list<Edge> *temp = new list<Edge>();
    temp->push_front(start);
    root = create_node(new string("*"), -g, temp);
    graph = root;
    
    //convert backbone to k-mer graph
    list<Node *> current;
    current.push_back(first);
    list<Node *> previous;
    previous.push_back(root);
    insert(first, backbone.substr(k,-1), string(backbone.size()-k, ')'), previous, current); // mid-low confidence to backbone
    
    return root;
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
	Node *current, *previous, *new_current;
    Edge edge;
	string backbone_path, reads_path, mappings_path, output_path;
    string sequence, quality, final_sequence, new_edge;
    list<tuple<string, string>> mappings;
    list<Node *> next_nodes, current_nodes, previous_nodes;
    set<Node *> visited;
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
    
    printf("adding sequences...\n");
    int strt,cnt = 0;
    previous_nodes.push_back(graph);
    current_nodes.push_back(graph->edges->front().next);
    for(int i=0; i<data.backbone.size()-g; i++){
        // advance forward: get set of nodes at next index
        if(i > current_nodes.front()->index){
            //printf("Advancing: [%d. %s] ->", previous_nodes.front()->index,  previous_nodes.front()->seq->c_str()); 
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
            //printf(" [%d. %s]\n", current_nodes.front()->index, current_nodes.front()->seq->c_str());
        }
        
        // if no mappings at current index, continue
        if(!data.mappings.count(i)) continue;
        
        mappings = data.mappings[i];
        for(list<tuple<string, string>>::iterator it = mappings.begin(); it != mappings.end(); ++it){
            
            //printf("Currently %d nodes at index %d\n", (int)current_nodes.size(), current_nodes.front()->index);
            sequence = get<0>(*it);
            quality = get<1>(*it);
            //printf("i = %d, current = [%d. %s]\n", i, current_nodes.front()->index, current_nodes.front()->seq->c_str());
            
            //trim sequence start to align with node if not adding to root
            if(i%g != 0){
                sequence = sequence.substr(g-i%g, -1);
                quality = quality.substr(g-i%g, -1);
            }
            
            //printf("mapping:\n    %s\n    %s\n", sequence.substr(0,25).c_str(), data.backbone.substr(i,25).c_str());
            //print_tree(graph, 20);
            
            // look for matching root node at current index
            bool found = false;
            for(list<Node *>::iterator n = current_nodes.begin(); n != current_nodes.end(); ++n){
                //printf("Searching: %s\n", (*n)->seq->c_str());
                if((*n)->seq->compare(sequence.substr(0,k))==0){
                    found = true;
                    current = *n;
                    break;
                }
            }
            
            // if sequence root matches, insert rest of sequence
            if(found){
                //printf("\nMatching root %s found at %d\n", current->seq->c_str(), current->index);
                insert(current, sequence.substr(k,-1), quality.substr(k,-1), previous_nodes, current_nodes);
            }
            // if sequence root does not match, create new edge from previous backbone node, then insert rest
            else{
                new_edge = get_edge(sequence, k);
                current = new_node(new_edge, current_nodes.front()->index);
                edge = create_edge(new string(new_edge), get_quality(quality.substr(0, new_edge.size())), current);
                previous_nodes.front()->edges->push_back(edge);
                current_nodes.push_back(current);
                //printf("\nRoot of %s not found at %d, created new node: %s \n", sequence.substr(0, 15).c_str(), current->index, current->seq->c_str());
                insert(current, sequence.substr(new_edge.size(), -1), quality.substr(new_edge.size(), -1), previous_nodes, current_nodes);
            }
        }
        cnt++;
        cout<<"[";
        int pos = 30*((float)(cnt))/data.mappings.size();
        for(int j=0; j<30; ++j){
            if(j<pos) cout<<"#";
            else cout <<" ";
        }
        cout<<"]"<<(int)(((float)cnt)/(data.mappings.size())*100)<<"%\r";
        cout.flush();
         
    }
    cout<<endl;
    
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
    printf("graph constructed\n");
    print_tree(graph, 27);
        
    printf("finding best path...\n");
    final_sequence = get_sequence(graph);
    
    ofstream output;
    output.open(output_path);
    output << ">" << "consensus_output\n";
    output << final_sequence << "\n\n";
    output.close();
    
    printf("new sequence written to %s\n", output_path.c_str());
	return 0;
}


