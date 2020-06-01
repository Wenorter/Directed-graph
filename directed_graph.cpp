//------------------------------
//DATA STRUCTURES AND ALGORITHMS
//ASSIGNMENT 1: DIRECTED GRAPGH
//------------------------------

//---------CONTENTS-------------
//PART 1: ADD, CHECK, REMOVE, RETRIEVE VERTEX
//PART 2: ADD, CHECK, REMOVE EDGE
//PART 3: DEGREE
//PART 5: CKECKING FACTORS
//PART 6: PRE-ORDER, IN-ORDER, AND POST-ORDER TRAVERSALS OF THE MST 
//------------------------------

#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

//---------
//LIBRARIES
//---------

#include <bits/stdc++.h>
#include <cstdlib>
#include <sstream>
#include <ctime>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <limits>
#include <utility> 
#include <map> 

//---------
//NAMESPACE
//---------

using namespace std;

//---------------------------------------------
//------------------CLASSES--------------------
//---------------------------------------------

//------------
//VERTEX CLASS
//------------

template <typename T> 
class vertex {
public:
	int id; //id of the vertex
	T weight; //weight of the vertex

	vertex(int vertID, T vertWeight) : id(vertID), weight(vertWeight){} // constructor for the vertex class
};

//--------------------
//DIRECTED GRAPH CLASS
//--------------------

template <typename T>
class directed_graph{

//---------------
//PRIVATE CLASSES
//---------------

private:
	//----------
	// vertex_list stores all vertices in the graph, as well as the vertices' weights
	//----------
	unordered_map<int, T> vertex_idweight; // each element is a pair(id, weight) for a vertex
	
	//----------
	// adj_list stores all edges in the graph, as well as the edges' weights.
	//----------
	unordered_map<int, unordered_map<int, T>> adj_list; // each element is a pair(vertex, the neighbours of this vertex)
															// each neighbour is also a pair (neighbour_vertex, weight for edge from vertex to neighbour_vertex)

//---------------
//PUBLIC CLASSES
//---------------

public:
	directed_graph(); //A constructor for directed_graph. The graph should start empty.
	~directed_graph(); //A destructor. Depending on how you do things, this may not be necessary.

	bool contains(const int&) const; //Returns true if the graph contains the given vertex_id, false otherwise.
	bool adjacent(const int&, const int&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

	void add_vertex(const vertex<T>&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const int&, const int&, const T&); //Adds a weighted edge from the first vertex to the second.

	void remove_vertex(const int&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const int&, const int&); //Removes the edge between the two vertices, if it exists.

	size_t in_degree(const int&) const; //Returns number of edges coming in to a vertex.
	size_t out_degree(const int&) const; //Returns the number of edges leaving a vertex.
	size_t degree(const int&) const; //Returns the degree of the vertex (both in edges and out edges).

	size_t num_vertices() const; //Returns the total number of vertices in the graph.
	size_t num_edges() const; //Returns the total number of edges in the graph.

	vector<vertex<T>> get_vertices(); //Returns a vector containing all the vertices.
	vector<vertex<T>> get_neighbours(const int&); //Returns a vector containing all the vertices reachable from the given vertex. The vertex is not considered a neighbour of itself.
	vector<vertex<T>> get_second_order_neighbours(const int&); // Returns a vector containing all the second_order_neighbours (i.e., neighbours of neighbours) of the given vertex.
															  // A vector cannot be considered a second_order_neighbor of itself.
	bool reachable(const int&, const int&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
	bool contain_cycles() const; // Return true if the graph contains cycles (there is a path from any vertices directly/indirectly to itself), false otherwise.

	vector<vertex<T>> depth_first(const int&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
	vector<vertex<T>> breadth_first(const int&) const; //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

	directed_graph<T> out_tree(const int&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges. This means every vertex in the tree is reachable from the root.

	vector<vertex<T>> pre_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of a pre-order traversal of the minimum spanning tree starting at the given vertex.
	vector<vertex<T>> in_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of an in-order traversal of the minimum spanning tree starting at the given vertex.
	vector<vertex<T>> post_order_traversal(const int&, directed_graph<T>&); // returns the vertices in ther visitig order of a post-order traversal of the minimum spanning tree starting at the given vertex.

	vector<vertex<T>> significance_sorting(); // Return a vector containing a sorted list of the vertices in descending order of their significance.
};

//---------------------------------------------
//------------------METHODS--------------------
//---------------------------------------------

//CONSTRUCTOR
template <typename T> 
directed_graph<T>::directed_graph() {}

//DESTRUCTOR
template <typename T> 
directed_graph<T>::~directed_graph() {}

//---------------------------------------------
//	PART 1: ADD, CHECK, REMOVE, RETRIEVE VERTEX
//---------------------------------------------

//CHECKING VERTEX
template <typename T> 
bool directed_graph<T>::contains(const int& u_id) const{
                                                            //NOTE TO SELF:
	if(vertex_idweight.find(u_id) != vertex_idweight.end()){//In C++ end() returns a value when the searching element is not found.
                                                            //Basically if the this sign does not appear the element has been found, this is why we use != and not ==
                                                            //If we write == instead of =! the if statement will result in infinite loop and search for noneexistent element forever.	
		return true;
	}

	return false;
}



//CHECKING IF VERTEX IS NEXT TO ANOTHER VERTEX
template <typename T> 
bool directed_graph<T>::adjacent(const int& u_id, const int& v_id) const{ 
	
	if (contains(u_id) && contains(v_id)){					//if both vertices are found
		auto vertex = adj_list.find(u_id);
		if(vertex != adj_list.end()){ 						//neighbour list checking while using the u_id 
			auto neighbours = vertex->second;
			if(neighbours.find(v_id) != neighbours.end()){ //checking the neighbour list in the list of their neighbours 
            
				return true;
			} 
		}
	}

	return false;
}

//ADDING NEW VERTEX
template <typename T>  
void directed_graph<T>::add_vertex(const vertex<T>& u) {

	if(!contains(u.id)){

		vertex_idweight.insert({u.id, u.weight}); 
		adj_list.insert({u.id, unordered_map<int, T>()});		// Step 1: add the new vertex to all_vertices
																// Step 2: add an entry for this vertex in adj_list but add no edge
		 		
	}

}

//---------------------------------------------
//  PART 2: ADD, CHECK, REMOVE EDGE
//---------------------------------------------

//ADDING NEW EDGE
template <typename T> 
void directed_graph<T>::add_edge(const int& u_id, const int& v_id, const T& uv_weight) {

	if(contains(u_id) && contains(v_id) && adj_list[u_id].find(v_id)==adj_list[u_id].end()){ 	// Step 1: make sure both vertices are in the graph	and the neighbour vertex is found																												// Step 2: make sure the edge is not already in the graph
		adj_list[u_id].insert({v_id, uv_weight}); 			                                    // Step 2: add this edge to adj_list
	}

}

//REMOVING VERTEX
template <typename T> 
void directed_graph<T>::remove_vertex(const int& u_id){ // remove the vertex, as well as all the incident edges
	
	vertex_idweight.erase(u_id); 			// Step 1: remove the vertex from all_vertices
	
	adj_list.erase(u_id); 					// Step 2: remove all edges starting from this vertex
	
	for (auto& x: adj_list){ 				// Step 3: iterate adj_list to remove all edges ending at this vertex
		x.second.erase(u_id);
	}
}

//REMOVING EDGE
template <typename T> 
void directed_graph<T>::remove_edge(const int& u_id, const int& v_id){								
		
		if(adjacent(u_id,v_id)){           //Step 1: checking if the vertices are adjacent.

			adj_list[u_id].erase(v_id);     //Step 1: removing the edge.
		}
	
}

//---------------------------------------------
//  PART 3: DEGREE
//---------------------------------------------

//FINDING THE AMOUNT OF EDGES POINTING TO A VERTEX
template <typename T> 
size_t directed_graph<T>::in_degree(const int& u_id) const { 

	size_t count = 0; //counter
	for(auto& x: adj_list){ // vertex id's
		if (adjacent(x.first, u_id)){ //checks if current id is adjacent to cpecified id
			count++;
		}		
	}
	return count;
}

//FINDING THE AMOUNT OF EDGES POINTING AWAY TO A VERTEX
template <typename T> 
size_t directed_graph<T>::out_degree(const int& u_id) const { 
	
	if(!contains(u_id)){ //if it doesn't contain 

		return 0; 
	}
	
	return adj_list.at(u_id).size(); //return number of edges
}

//FINDING THE AMOUNT OF EDGES POINTING TO THE VERTEX AND AWAY FROM THE VERTEX
template <typename T>
size_t directed_graph<T>::degree(const int& u_id) const { 
	
	return in_degree(u_id) + out_degree(u_id); 
}

//---------------------------------------------
//  PART 3: NUMS
//---------------------------------------------

//COUNTS THE NUMBER OF EGES IN A VERTEX
template <typename T> 
size_t directed_graph<T>::num_vertices() const { 

	int count = 0;
	for (int i=0; i < vertex_idweight.size(); i++) {
		for (int i=0; i < adj_list.size(); i++) {
			count++;
		}

	 return count;

	}

	return false; //don't change that stupid!
}


//COUNTS THE NUMBER OF EGES IN A VERTEX
template <typename T> 
size_t directed_graph<T>::num_edges() const { 

	size_t count = 0;
	for (auto& x: adj_list){ // x == pair<int, unordered_map<int,T>>
		count += x.second.size(); // x.second == unordered_map<int, T>
	}

	return count;
}

//---------------------------------------------
//  PART 4: GET METHODS
//---------------------------------------------

//GETS THE ID AND WEIGHT OF EACH VERTICE
template <typename T> 
vector<vertex<T>> directed_graph<T>::get_vertices() {

	vector<vertex<T>> vert;

	for(auto& x: vertex_idweight){ 								// Step 1: iterate vertex_list to get all vertex_ids
		vert.push_back(vertex<T>(x.first, x.second));           // Step 2: build a vertex class for each vertex_id
	}

	return vert;                                                // Step 3: return a vector of the vertex classes for all vertex id's

}

//GETS THE NEIGHBOURS BY THE VERTEX ID KEY
template <typename T> 
vector<vertex<T>> directed_graph<T>::get_neighbours(const int& u_id){

	vector<vertex<T>> vert;

	if(contains(u_id)){ 													// Step 1: make sure the vertex is in the graph
		for (auto x: adj_list[u_id]){ 										// Step 2: find all edges starting from the vertex of u_id
			vert.push_back(vertex<T>(x.first, vertex_idweight[x.first]));	// Step 3: add the end_node of each edge to the result
		}
	}

	return vert;

}

//GETS THE NEIGHBOURS OF THE NEIGHBOURS BY THE VERTEX ID KEY
template <typename T> 
vector<vertex<T>> directed_graph<T>::get_second_order_neighbours(const int& u_id) { 
	
	vector<vertex<T>> vert;
	map<int,vertex<T>> list; //controls that you only have unique keys

	if(contains(u_id)){ 														
		for (auto x: adj_list[u_id]){ 											
			for (auto f: get_neighbours(x.first)){
				if(f.id != u_id){
				 	list.insert({f.id, f}); //insert doesn't allow duplicates.
				}			
			}																			
		}

		for(auto i: list){ //copying from map to vector
			vert.push_back(i.second);
		}
	}
	
	return vert;
}

//---------------------------------------------
//  PART 5: CKECKING FACTORS
//---------------------------------------------

//CHECKS IF SECOND VERTEX IS REACHABLE FROM THE FIRST
template <typename T> 
bool directed_graph<T>::reachable(const int& u_id, const int& v_id) const { 
	
	if (adjacent(u_id, v_id)) return true; //for optimization.

	for (auto v: breadth_first(u_id)){ //uses bft to traverse through the graph
		if(v.id == v_id){ //if id's equal
			return true;
		}
	}
	return false; 
}

//CHECKS IF THE GRAPH CONTAINS CYCLE EDGE
template <typename T> 
bool directed_graph<T>::contain_cycles() const { 
	for(auto x: vertex_idweight){
		auto found = adj_list.find(x.first);
		if (found == adj_list.end()) //if id is found continues to travel
			continue;
		for (auto v: found -> second ){
			
			if(x.first != v.first && reachable(v.first, x.first)){//checks of vertex id's do not have each other in their neighbour list

				return true;
			}
		}

	}
	return false; 
}

//DEPTH FIRST TRAVELSAL
template <typename T> 
vector<vertex<T>> directed_graph<T>::depth_first(const int& u_id) { 
	//ITERATIVE METHOD
    //NOTE TO SELF: traversal is only complete when stack is empty in the end
    // travels -> o -> o -> o -> o straight away without 'looking around' first in comparison to BFT
	vector<vertex<T>> dft_vect;
	vector<int> list;
	
	if (contains(u_id)){

		stack<int> s; 
		s.push(u_id); //starting node pushed to stack
		
		while (!s.empty()){
			auto id = s.top();
			s.pop();
			bool visited = find(list.begin(),list.end(), id) != list.end(); 
			if (visited) continue; //marking with visited
			list.push_back(id); //put id's into list
			for (auto v: get_neighbours(id)){

				s.push(v.id);					 							
			}
		}
	}
	
	for(auto i: list) //copying from map to vector
		{ 
			dft_vect.push_back({i,vertex_idweight.find(i) -> second});
		}


	return dft_vect;
}

//BREADTH FIRST TRAVELSAL
template <typename T> 
vector<vertex<T>> directed_graph<T>::breadth_first(const int& u_id) const{ 
    
	//ITERATIVE METHOD
    // is stationary, remains on the same node, it's looking for other nodes around and visits them.
	vector<vertex<T>> bft_vect;
	vector<int> list;
	
	if (contains(u_id)){

		queue<int> q;   
		q.push(u_id); //starting node pushed to quere.
		
		while (!q.empty()){
			auto id = q.front(); //gets the first id from a queue into a variable
			q.pop(); //removes it from queue
			bool visited = find(list.begin(),list.end(), id) != list.end(); //marking with visited
			if (visited) continue; //continues to visit
			list.push_back(id); //records id when visits
			auto found = adj_list.find(id); 
			if (found == adj_list.end()) //duplicate check
				continue;
			for (auto v: found -> second ){
				q.push(v.first);

			}
		}
	}
	
	for(auto i: list) //copying from map to vector
		{ 
			bft_vect.push_back({i,vertex_idweight.find(i) -> second});
		}

	return bft_vect; 
}

//SPANNING TREE FROM A SPECIFIED ROOT POINT
template <typename T> 
directed_graph<T> directed_graph<T>::out_tree(const int& u_id) {
//NOTE TO SELF: CAN'T HAVE EDGE DEGREE MORE THAN TWO.
// this tree is produced by getting the neighbours of the vertex by the key u_id.
// it's important to consider the less weightened edge, from EVERY visited node.
// you cant have edge cycles.
//NOTE: I used reverse deleta algorithm

	directed_graph<T> tree; //initial tree size.
	//vector<int> unprocessed; //stack of unprocessed nodes.

	if (contains(u_id)){

		//queue<int> s;
		//s.push(u_id); //starting node pushed to stack
		multimap<T, pair<int,int>> edges; //multimap for storing edges

		for (auto i : adj_list){ //copying the graph
			if(reachable(u_id, i.first)){ //checks if the id is reachable

				
				tree.add_vertex(vertex(i.first, vertex_idweight[i.first])); //adds the root

				for (auto v: get_neighbours(i.first)){
					
					edges.insert({adj_list[i.first][v.id], {i.first, v.id}}); //insert edge id and weight into edge list 
					tree.add_vertex(vertex(v.id, vertex_idweight[v.id])); //adds vertices in the tree
					tree.add_edge(i.first, v.id, adj_list[i.first][v.id]); //adds edge in a tree

				}
			}
		}

		auto size = tree.depth_first(u_id).size(); //gets the number of verdices in the graph
	
		for(auto it = edges.rbegin(); it != edges.rend(); ++it){ //iterates through edge list
			int uID = it -> second.first; //first id
			int vID = it -> second.second; // second id
			
			tree.remove_edge(uID, vID); //removes the edge between first id and second id
			//if the number of vertices connected now are less then connected before, means that the graph is disconnected
			if(size != tree.depth_first(u_id).size()){

				tree.add_edge(uID, vID, it -> first); //adds edge back to reconnect the graph

			}
		}
	}

	return tree;	 

}


//-------------------------------------------------------------------------------------
//	PART 6: pre-order, in-order, and post-order traversals of the Minimum Spanning Tree
//-------------------------------------------------------------------------------------

//PRE-ORDER TRAVERSAL
template <typename T> 
vector<vertex<T>> directed_graph<T>::pre_order_traversal(const int& u_id, directed_graph<T>& mst){ 
//1. Visit Node
//2. Traverse left
//3. Traverse right
//4. O(n) time complexity
//----------------
//PSEUDOCODE
//----------------
//1. if node == null then return
//2. visit(node)
//3. preorder(node.left)
//4. preorder(node.right)

	return mst.depth_first(u_id); //i just return the values using dft because it has idetical steps when it visits and traverses.
}

//IN-ORDER TRAVERSAL
template <typename T> 
vector<vertex<T>> directed_graph<T>::in_order_traversal(const int& u_id, directed_graph<T>& mst) { 
//1. Traverse left
//2. Visit node
//3. Traverse right
//In-order visits the vertex between the left traversal and the right traversal.
//----------------
//PSEUDOCODE
//----------------
//1. if node == null then return
//2. inorder(node.left)
//3. visit(node)
//4. inorder(node.right)
	
	vector<vertex<T>> iot_vect;
	auto tree = mst.out_tree(u_id);
	if (contains(u_id)){
		
		multimap<T, int> m;
		for (auto vec : tree.get_neighbours(u_id)){
			
			m.insert({tree.adj_list[u_id][vec.id], vec.id});
		}	
		auto v = m.begin();	
		if (v != m.end()){

			auto vert = in_order_traversal(v->second, mst);
			iot_vect.insert(iot_vect.end(), vert.begin(), vert.end());
				
			++v;				 							
		}
		iot_vect.push_back({u_id, tree.vertex_idweight[u_id]});

		for (; v != m.end(); ++v){
			
			auto vert = in_order_traversal(v->second, mst);
			iot_vect.insert(iot_vect.end(), vert.begin(), vert.end());
				
						 							
		}
	}

	return iot_vect;

}

//POST-ORDER TRAVERSAL
template <typename T> //IN PROGRESS
vector<vertex<T>> directed_graph<T>::post_order_traversal(const int& u_id, directed_graph<T>& mst) { 
//1. Traverse left
//2. Traverse right
//3. Visit node
//Records node at the end of the branch first
//----------------
//PSEUDOCODE
//----------------
//1. if node == null then return
//2. postorder(node.left)
//3. postorder(node.right)
//4. visit(node)
vector<vertex<T>> pot_vect;
	auto tree = mst.out_tree(u_id); //tree mst
	if (contains(u_id)){
		multimap<T, int> edges;
		for (auto vec : tree.get_neighbours(u_id)){
			
			edges.insert({tree.adj_list[u_id][vec.id], vec.id});
		}		
		for (auto v = edges.begin(); v != edges.end(); ++v){
			
			auto vert = post_order_traversal(v->second, mst);
			pot_vect.insert(pot_vect.end(), vert.begin(), vert.end());
						 							
		}
		pot_vect.push_back({u_id, tree.vertex_idweight[u_id]}); //3. visit nodes
	}

	return iot_vect; 
}

//SORTING ALGORITHM
template <typename T> //DONE
vector<vertex<T>> directed_graph<T>::significance_sorting() { 
	
	//https://www.geeksforgeeks.org/sorting-algorithms/
	//Bubble sort

	vector<vertex<T>> vert;

	for (auto v : vertex_idweight){
		vert.push_back({v.first, v.second}); //gets all the vertices
	}

	for (size_t i = 0; i < vert.size() - 1; i++) //for all vertices in a descending order

		for (size_t j = 0; j < vert.size() - 1 - i; j++){
			if (vert[j].id < vert[j+1].id){ //if first id is less than second id

				swap(vert[j].id,vert[j+1].id); //id swap
				swap(vert[j].weight,vert[j+1].weight); //weight swap

			}
		}

	

	return vert; 
}
#endif
//----------------------END OF THE ASSIGNMENT 1---------------------