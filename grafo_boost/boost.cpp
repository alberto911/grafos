#include <iostream>
#include <chrono>
#include "graph.h"

using namespace std;
using namespace std::chrono;

int main() {
	// Directed graph
	MyGraph<DirectedGraph> dg;
	int V = 14, E = 24;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 1; i <= V; ++i)
		dg.addVertex(i);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Add vertices: " << duration << endl;

	t1 = high_resolution_clock::now();	
	dg.loadEdges("edges.txt", E);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Add edges: " << duration << endl;

	t1 = high_resolution_clock::now();
	dg.DFS(1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "DFS: " << duration << endl;

	t1 = high_resolution_clock::now();
	dg.BFS(1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "BFS: " << duration << endl;

	t1 = high_resolution_clock::now();
	dg.Dijkstra(1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Dijkstra: " << duration << endl;

	t1 = high_resolution_clock::now();
	dg.FloydWarshall();
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Floyd-Warshall: " << duration << endl;
	
	//Undirected graph
	MyGraph<UndirectedGraph> ug;
	int E2 = 19;

	for (int i = 1; i <= V; ++i)
		ug.addVertex(i);
	ug.loadEdges("edges2.txt", E2);
	
	t1 = high_resolution_clock::now();
	ug.Prim(1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Prim: " << duration << endl;

	t1 = high_resolution_clock::now();
	ug.Kruskal();
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Kruskal: " << duration << endl;
	
	return 0;
}

