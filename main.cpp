#include "directed_graph.hpp"

int main() {

	directed_graph<int> g1;

	//Assignment vertices 1-5
	vertex<int> v0(1, 800); //Node A
	vertex<int> v1(2, 3000); //Node B
	vertex<int> v2(3, 400); //Node C
	vertex<int> v3(4, 710); //Node D
	vertex<int> v4(5, 221); //Node E
	//BFS example vertices

	g1.add_vertex(v0); // add vertex class
	g1.add_vertex(v1);
	g1.add_vertex(v2);
	g1.add_vertex(v3);
	g1.add_vertex(v4);


	//g1.remove_vertex(1); // delete vertex by id
	cout << "SUCCESS!" << endl;
	cout << "---------------" << endl;
	cout << "OUTPUT START" << endl;
	cout << "---------------" << endl;
	vector<vertex<int>> vertex_list = g1.get_vertices();
	cout << "All vertices: ";
	for (vertex<int> vt : vertex_list) {
	 	cout << "(" << vt.id << ", " << vt.weight << ") ";
	}
	cout << endl;

	//Assignment edges
	// g1.add_edge(1, 2, 6); //A->B
	// g1.add_edge(1, 3, 9); //A->C
	// g1.add_edge(2, 5, 3); //B->E
	// g1.add_edge(3, 4, 4); //C->D
	// g1.add_edge(4, 3, 7); //D->C
	// g1.add_edge(4, 1, 1); //D->A
	// g1.add_edge(4, 5, 5); //D->E

	g1.add_edge(1, 2, 6); 
   	g1.add_edge(1, 3, 9);
   	g1.add_edge(2, 5, 3);
   	g1.add_edge(3, 4, 4);
	g1.add_edge(4, 1, 1);
	g1.add_edge(4, 3, 7);
	g1.add_edge(4, 5, 5);
   	
	
	//g1.remove_vertex(4);

	cout << "All neighbours of 1: ";
	vector<vertex<int>> neighbour_list = g1.get_neighbours(1);
	for (vertex<int> nb : neighbour_list) {
	 	cout << "(" << nb.id << ", " << nb.weight << ") ";
	}
	cout << endl;

	cout << "Second Order Neighbours of 1: ";
	vector<vertex<int>> second_neighbour_list = g1.get_second_order_neighbours(1);
	for (vertex<int> snb : second_neighbour_list) {
	 	cout << "(" << snb.id << ", " << snb.weight << ") ";
	}
	cout << endl;

	cout << "Total number of edges: ";
	cout << g1.num_edges() << endl;
	cout << endl;

	cout << "DFS: ";
	vector<vertex<int>> DFS = g1.depth_first(1);
	for (vertex<int> dfs : DFS) {
	 	cout << "(" << dfs.id << ", " << dfs.weight << ") ";
	}
	cout << endl;

	cout << "BFS: ";
	vector<vertex<int>> BFS = g1.breadth_first(1);
	for (vertex<int> bfs : BFS) {
	 	cout << "(" << bfs.id << ", " << bfs.weight << ") ";
	}
	cout << endl;

	cout << "Contains Cycles? ";
    if (g1.contain_cycles() == 0){
        cout << "No" << endl;
    }
    else {
        cout << "Yes" << endl;
    }

	cout << "Pre Order Travelsal: ";
	auto mst = g1.out_tree(1);
	vector<vertex<int>> PROT = g1.pre_order_traversal(1, mst);
	
	for (vertex<int> prot : PROT) {
	 	cout << "(" << prot.id << ", " << prot.weight << ") ";
	}

	cout << endl;
	cout << "In Order Travelsal: ";
	vector<vertex<int>> IOT = g1.in_order_traversal(1, mst);
	
	for (vertex<int> iot : IOT) {
	 	cout << "(" << iot.id << ", " << iot.weight << ") ";
	}
	cout << endl;

	cout << "Post Order Travelsal: ";
	vector<vertex<int>> POT = g1.post_order_traversal(1, mst);
	
	for (vertex<int> pot : POT) {
	 	cout << "(" << pot.id << ", " << pot.weight << ") ";
	}
	cout << endl;

	cout << "Significance Sorting: ";
	vector<vertex<int>> SORT = g1.significance_sorting();
	for (vertex<int> sort : SORT) {
	 	cout << "(" << sort.id << ", " << sort.weight << ") ";
	}
	cout << endl;
	cout << "---------------" << endl;
	cout << "OUTPUT END" << endl;
}
