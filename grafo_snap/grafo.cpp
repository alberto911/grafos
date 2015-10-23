#include <iostream>
#include <string>
#include <fstream>
#include <stack>
#include <queue>
#include <algorithm>
#include <utility>
#include <limits>
#include <iomanip>
#include <chrono>
#include "Snap.h"

#undef max

using namespace std::chrono;

typedef TPt<TNodeEDatNet<TInt, TInt> > DGraph;

void addVertex(DGraph g, int id) {
	g->AddNode(id);
}

void addEdge(DGraph g, int source, int destination, int weight) {
	g->AddEdge(source, destination, weight);
}

void removeVertex(DGraph g, int id) {
	g->DelNode(id);
}

void removeEdge(DGraph g, int src, int dst) {
	g->DelEdge(src, dst);
}

void loadDirected(DGraph g, std::string filename, int E) {
	int u, v, weight;
	
	std::ifstream file(filename);
	if (file.is_open()) {
		for (int i = 0; i < E; ++i) {
			file >> u >> v >> weight;
			addEdge(g, u, v, weight);
		}
	}
}

void loadUndirected(DGraph g, std::string filename, int E) {
	int u, v, weight;
	
	std::ifstream file(filename);
	if (file.is_open()) {
		for (int i = 0; i < E; ++i) {
			file >> u >> v >> weight;
			addEdge(g, u, v, weight);
			addEdge(g, v, u, weight);
		}
	}
}

void DFS(DGraph g, int source) {
	int n = g->GetNodes();
	bool visited[n];
	std::stack<DGraph::TObj::TNodeI> pila;

	for (int i = 0; i < n; ++i) {
		visited[i] = false;
	}

	DGraph::TObj::TNodeI u = g->GetNI(source);
	pila.push(u);

	while (!pila.empty()) {
		DGraph::TObj::TNodeI v = pila.top();
		pila.pop();
		
		std::cout << v.GetId() << " ";
		for (int e = v.GetOutDeg()-1; e >= 0; e--) {
			int id = v.GetOutNId(e);
			if (!visited[id-1]) {
				visited[id-1] = true;
				pila.push(g->GetNI(id));
			}
		}
	}
	std::cout << std::endl;
}

void BFS(DGraph g, int source) {
	int n = g->GetNodes();
	bool visited[n];
	std::queue<DGraph::TObj::TNodeI> cola;

	for (int i = 0; i < n; ++i) {
		visited[i] = false;
	}

	DGraph::TObj::TNodeI u = g->GetNI(source);
	std::cout << source << " ";
	visited[source-1] = true;
	cola.push(u);

	while (!cola.empty()) {
		DGraph::TObj::TNodeI v = cola.front();
		cola.pop();
		
		for (int e = 0; e < v.GetOutDeg(); e++) {
			int id = v.GetOutNId(e);
			if (!visited[id-1]) {
				std::cout << id << " ";
				visited[id-1] = true;
				cola.push(g->GetNI(id));
			}
		}
	}
	std::cout << std::endl;
}

void FloydWarshall(DGraph g) {	
	int n = g->GetNodes();
	int distances[n][n];

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j)		
				distances[i][j] = 0;
			else
				distances[i][j] = std::numeric_limits<int>::max() / 2;
		}	
	}

	for (DGraph::TObj::TEdgeI ei = g->BegEI(); ei < g->EndEI(); ei++) {
		distances[ei.GetSrcNId()-1][ei.GetDstNId()-1] = ei.GetDat();
	}

	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (distances[i][j] > distances[i][k] + distances[k][j])
					distances[i][j] = distances[i][k] + distances[k][j];
			}
		}
	}

	// Imprimir
	std::cout << "All-pairs shortest paths\n";
	std::cout << "       ";
	for (DGraph::TObj::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++)
		std::cout << std::setw(5) << NI.GetId();
	std::cout << std::endl;
	
	DGraph::TObj::TNodeI NI = g->BegNI();
	for (int i = 0; i < n; ++i, NI++) {
		std::cout << std::setw(3) << NI.GetId() << " -> ";
		for (int j = 0; j < n; ++j) {
			if (distances[i][j] == (std::numeric_limits<int>::max)() / 2)
				std::cout << std::setw(5) << "inf";
			else
				std::cout << std::setw(5) << distances[i][j];
		}
		std::cout << std::endl;
	}
}

void Dijkstra(DGraph g, int source) {
	int n = g->GetNodes();
	int distances[n];
	int parent[n];
	std::priority_queue<std::pair<int, DGraph::TObj::TNodeI> > q;
	
	distances[source-1] = 0;
	parent[source-1] = source;

	for (DGraph::TObj::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
		int id = NI.GetId();
		if (id != source) {
			distances[id-1] = std::numeric_limits<int>::max() / 2;
			parent[id-1] = -1;
		}
		q.push(std::make_pair(distances[id-1], NI));
	}

	while (!q.empty()) {
		DGraph::TObj::TNodeI u = q.top().second;
		int dist = q.top().first;
		q.pop();

		for (int e = 0; e < u.GetOutDeg(); e++) {
			int idV = u.GetOutNId(e);
			int alt = dist + g->GetEI(u.GetId(), idV).GetDat();
			if (alt < distances[idV-1]) {
				distances[idV-1] = alt;
				parent[idV-1] = u.GetId();
				q.push(std::make_pair(alt, g->GetNI(idV)));
			}
		}
	}

	std::cout << "Shortest paths from " << source << ":\n";
	for (DGraph::TObj::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
		int id = NI.GetId();		
		if (parent[id-1] == -1)
			std::cout << "Vertex " << id << " can't be reached\n";
		else {
			std::cout << "Distance to " << id << " = " << distances[id-1] << ", ";
			std::cout << "parent = " << parent[id-1] << std::endl;
		}
	}
}

void Prim(DGraph g, int source) {
	int n = g->GetNodes();
	int distances[n];
	int parent[n];
	std::priority_queue<std::pair<int, DGraph::TObj::TNodeI>, std::vector< std::pair<int, DGraph::TObj::TNodeI> >, std::greater< std::pair<int, DGraph::TObj::TNodeI> > > q;

	distances[source-1] = 0;
	parent[source-1] = source;

	for (DGraph::TObj::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
		int id = NI.GetId();
		if (id != source) {
			distances[id-1] = std::numeric_limits<int>::max() / 2;
			parent[id-1] = -1;
		}
		q.push(std::make_pair(distances[id-1], NI));
	}

	while (!q.empty()) {
		DGraph::TObj::TNodeI u = q.top().second;
		q.pop();
		distances[u.GetId()-1] = 0;

		for (int e = 0; e < u.GetOutDeg(); e++) {
			int idV = u.GetOutNId(e);
			int alt = g->GetEI(u.GetId(), idV).GetDat();
			if (distances[idV-1] != 0 && alt < distances[idV-1]) {
				distances[idV-1] = alt;
				parent[idV-1] = u.GetId();
				q.push(std::make_pair(alt, g->GetNI(idV)));
			}
		}
	}

	std::cout << "Minimum spanning tree:" << std::endl;
	for (DGraph::TObj::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
		int id = NI.GetId();
		std::cout << "Parent of " << id << " = ";
		if (parent[id-1] != id)
			std::cout << parent[id-1] << std::endl;
		else
			std::cout << "no parent" << std::endl;
	}
}

int findSet(int x, int parent[]) {
	if (parent[x-1] != x)
		parent[x-1] = findSet(parent[x-1], parent);
	return parent[x-1];
}

void setUnion(int u, int v, int parent[]) {
	u = findSet(u, parent);
	v = findSet(v, parent);
	parent[u-1] = v;
}

void Kruskal(DGraph g) {
	int n = g->GetNodes();
	int parent[n];
	std::vector<std::pair<int, DGraph::TObj::TEdgeI> > edges;
	std::queue<DGraph::TObj::TEdgeI> tree;

	for (DGraph::TObj::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
		int id = NI.GetId();
		parent[id-1] = id;
	}

	for (DGraph::TObj::TEdgeI ei = g->BegEI(); ei < g->EndEI(); ei++)
		edges.push_back(std::make_pair(ei.GetDat(), ei));

	std::sort(edges.begin(), edges.end());

	for (unsigned int i = 0; i < edges.size(); ++i) {
		int u = edges[i].second.GetSrcNId();
		int v = edges[i].second.GetDstNId();
		if (findSet(u, parent) != findSet(v, parent)) {
			tree.push(edges[i].second);
			setUnion(u, v, parent);
		}
	}

	while (!tree.empty()) {
		std::cout << "Edge between " << tree.front().GetSrcNId() << " and " << tree.front().GetDstNId();
		std::cout << ", weight = " << tree.front().GetDat() << std::endl;
		tree.pop();
	}
}

int main() {
	DGraph dg = DGraph::TObj::New();
	int V = 14, E = 24, E2 = 19;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 1; i <= V; ++i)
		addVertex(dg, i);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "Add vertices: " << duration << std::endl;

	t1 = high_resolution_clock::now();
	loadDirected(dg, "edges.txt", E);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "Add edges: " << duration << std::endl;
	
	t1 = high_resolution_clock::now();
	DFS(dg, 1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "DFS: " << duration << std::endl;

	t1 = high_resolution_clock::now();
	BFS(dg, 1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "BFS: " << duration << std::endl;

	t1 = high_resolution_clock::now();
	Dijkstra(dg, 1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "Dijkstra: " << duration << std::endl;

	t1 = high_resolution_clock::now();
	FloydWarshall(dg);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "Floyd-Warshall: " << duration << std::endl;

	DGraph ug = DGraph::TObj::New();
	for (int i = 1; i <= V; ++i)
		addVertex(ug, i);
	loadUndirected(ug, "edges2.txt", E2);

	t1 = high_resolution_clock::now();
	Prim(ug, 1);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "Prim: " << duration << std::endl;

	DGraph dg2 = DGraph::TObj::New();
	for (int i = 1; i <= V; ++i)
		addVertex(dg2, i);
	loadDirected(dg2, "edges2.txt", E2);

	t1 = high_resolution_clock::now();
	Kruskal(dg2);
	t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "Kruskal: " << duration << std::endl;

  return 0;
}
