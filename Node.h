/* 
 * This file contains classes that make the k-mer graph
 * 
 * Author: Tin Kuculo
 */
 
#include "Utils.h"

class Node;
extern int k, g;
extern Node *graph;

/*
 * Edge class:
 *      seq = pointer to g-length (or more/less) string sequence
 *      weight = weight of edge (double)
 *      next = pointer to next Node
 */
class Edge {
public:
	string *seq;
	double weight;
	Node *next;
};

 /*
  * Node class: 
  *     seq = pointer to k-length string sequence
  *     index = position of node in backbone (int)
  *     edges = pointer to list of edges, points to next edges
  */
class Node {
public:
	string * seq;
	list<Edge> * edges;
	int index; //position in backbone

    bool operator==(const Node &n1){
        return(seq->data()==n1.seq->data() && index == n1.index);
    }
};

Node *create_node(string seq, int index, list<Edge> *edges){
    Node * node = (Node *)malloc(sizeof(Node));
    node->seq = new string(seq);
    node->index = index;
    node->edges = edges;
    return node;
}

Edge create_edge(string seq, double weight, Node *next){
    Edge *edge = (Edge *)malloc(sizeof(Edge));
    edge->seq = new string(seq);
    edge->weight = weight;
    if(next)
        edge->next = next;
    return *edge;
}   

/* Adds new sequence to graph. Inserts new nodes aligned to backbone.
 * 
 * Input: current node, sequence to input, quality of sequence, 
 *        list of nodes at previous index, list of nodes at current index
 * 
 * Authors: Ana Brassard & Tin Kuculo
 */
void insert(Node *current, string sequence, string quality, list<Node *> previous_nodes, list<Node *> current_nodes ){
    Node *node, *next_backbone;
    set<Node *> visited;// = *new set<Node *>();
    map<tuple<int, string>, Node *> nodes;
    map<int, Edge> edges;
    string seq, q, new_edge, new_node_seq;
    int index = current->index;
    int i=0;
    
    // prepare map of all nodes to add, fill with new nodes, stop when less than g bases left to add
    nodes[make_tuple(index, current->seq->data())] = current;
    while(get_edge(sequence.substr(i,-1), g).size() >= g){
        seq = get_edge(sequence.substr(i,-1), g);
		q = quality.substr(i, seq.size());
        
        new_edge = remove_char(seq, '_');
        edges[index] = create_edge(new_edge, get_quality(q), nullptr);
        
        new_node_seq = seq.substr(seq.size()-k, k);
        node = create_node(new_node_seq, index+g, new list<Edge>());
        nodes[make_tuple(node->index, new_node_seq)] = node;
       
        index += g;
        i+=seq.size();
    }
    
    // go through graph, if node to add found, replace new node with found node in map
    while(true){
        visited.clear();
        previous_nodes = current_nodes;
        current_nodes.clear();
        for(list<Node *>::iterator n = previous_nodes.begin(); n!= previous_nodes.end(); ++n){
            
            // if current node is in map, replace new node with current node
            if(nodes.count(make_tuple((*n)->index, (*n)->seq->data())))
                nodes[make_tuple((*n)->index, (*n)->seq->data())] = *n;
            
            // stop search if last index, also remember current node as next backbone node
            if((*n)->index == index){
                next_backbone = (*n);
                visited.clear();
                break;
            }
            
            // add next nodes as current nodes (continue search)
            for(list<Edge>::iterator e = (*n)->edges->begin(); e!= (*n)->edges->end(); ++e){
                if(!visited.count(e->next)){
                    current_nodes.push_back(e->next);
                    visited.insert(e->next);
                }
            }
        }
        if(visited.empty()) break;
    }

    // convert node map to node array
    Node *node_list[nodes.size()];
    int j=0;
    for(map<tuple<int, string>, Node *>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
        node_list[j++] = it->second;
        
    // if last node is new (not found), add edge that connects to backbone
    if(node_list[j-1]->edges->empty())
        edges[index] = create_edge("", 0, next_backbone);
    
    /* commented: ignore tail of sequence
    else{
        seq = sequence.substr(i, -1);
		q = quality.substr(i, -1);
        new_edge = remove_char(seq, '_');
        edges[index] = create_edge(new_edge, get_quality(q), next_backbone);
    }
    */
    
    Edge edge;
    bool found;
    
    // connect nodes and edges
    for(i=0; i<nodes.size()-1; i++){
        node = node_list[i];
        edge = edges[node->index];
        found=false;
        for (list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it) {
            // if matching edge found, add quality to weight
            if (edge.seq->compare(it->seq->data()) == 0){
               it->weight += edge.weight;
               found=true;
               break;
            }
        }
        // else, create new edge, connect to previous and next node
        if(!found){
            edge.next = node_list[i+1];
            node->edges->push_back(edge);
        }
    }
    // if there is an extra edge (if a new edge that leads to backbone was added)
    if(edges.count(node_list[i]->index))
        node_list[i]->edges->push_back(edges[node_list[i]->index]);

}

// Prints tree from given root to max index (BFS)
void print_tree(Node *root, int maxindex){
    Node *node;
    list<Node *> queue = *new list<Node *>();
    set<Node *> visited = *new set<Node *>();
    set<Node *> added = *new set<Node *>();
    
    queue.push_back(root);
    while(!queue.empty()){
        node = queue.front();
        queue.pop_front();
        if(node->index>maxindex)
            continue;
        printf("\n%s[%d. %s]\n",string(node->index+g,' ').c_str(), node->index, node->seq->c_str());
        //printf("\n[%d. %s]\n", node->index, node->seq->c_str());
        for(list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it) {
            if(!visited.count(node))
                printf("%s  --%s(%f)-->[%d. %s]\n",string(node->index+g,' ').c_str(),it->seq->c_str(), 
                                            it->weight, it->next->index, it->next->seq->c_str());
                // printf(" --%s(%f)-->[%d. %s]\n",it->seq->c_str(), it->weight, it->next->index, it->next->seq->c_str());
                
            if(!visited.count(it->next) && !added.count(it->next)){
                queue.push_back(it->next);
                added.insert(it->next);
            }
                
        }
        visited.insert(node);
    }
}

