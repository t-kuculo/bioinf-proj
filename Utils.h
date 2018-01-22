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
class Node;
extern int k, g;
extern Node *graph;
Node *find_node(Node *, int, string);

// Edge class, contains g-length sequence and a weight
class Edge {
public:
	string *seq;
	float weight;
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
		float q = 0;
		// convert quality string to int
		for (int i = 0; i < quality.size(); i++)
			q += int(quality[i])-int('!');
        q /= quality.size();
        
        for (list<Edge>::iterator it = this->edges->begin(); it != this->edges->end(); ++it) {
            // if matching edge found, add quality to weight
            if (seq.compare(it->seq->c_str()) == 0){
                it->weight += q;
                return;
            }
        }
            
        // else, create new edge, create new node if necessary 
        string node_seq = seq.substr(seq.size() - k, k);
        string edge_seq = seq;
        
        Edge * new_edge = (Edge *)malloc(sizeof(Edge));
        new_edge->seq = new string(edge_seq);
        new_edge->weight = q;
        void * temp = find_node(graph, index, node_seq);
        
        if(temp == nullptr){
            Node * new_node = (Node *)malloc(sizeof(Node));
            new_node->seq = new string(node_seq);
            new_node->index = index;
            new_node->edges = new list<Edge>();
            new_edge->next = new_node;
        }
        else
            new_edge->next = (Node *)temp;
    
        //printf("new edge: [%s]--%s-->[%s]\n", this->seq->c_str(), new_edge->seq->c_str(), new_edge->next->seq->c_str());
        edges->push_back(*new_edge);
        
    }

    // Function to advance forward, takes edge sequence
    Node *next(string edge){
        for (list<Edge>::iterator it = this->edges->begin(); it != this->edges->end(); ++it) {
            if (edge.compare(it->seq->c_str()) == 0){
                return it->next;
            }
        }
        printf("Error! Edge not found.\n");
        return 0;
    }

};

Node * new_node(string seq, int index){
    Node * new_node = (Node *)malloc(sizeof(Node));
    new_node->seq = new string(seq);
    new_node->index = index;
    new_node->edges = new list<Edge>();
    return new_node;
}
    
    
// DFS from given node: find correct index and sequence
Node *find_node(Node *current, int i, string seq){
    Node *node;
    list<Node *> queue = *new list<Node *>();
    set<Node *> visited = *new set<Node *>();
    
    queue.push_front(current);
    while(!queue.empty()){
        node = queue.front();
        queue.pop_front();
        if(node->index==i && node->seq->compare(seq)==0){
            return node;
        }
        for(list<Edge>::iterator it = node->edges->begin(); it != node->edges->end(); ++it) {
            if(!visited.count(it->next) && it->next->index <= i)
                queue.push_front(it->next);
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
 