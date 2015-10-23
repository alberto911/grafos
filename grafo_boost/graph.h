#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/property_map/property_map.hpp>

using namespace boost;

struct VData {
	int id;
	VData() { }
	VData(int i) : id(i) { }
};

typedef adjacency_list < listS, vecS, directedS, VData, property < edge_weight_t, float > > DirectedGraph;
typedef adjacency_list < listS, vecS, undirectedS, VData, property < edge_weight_t, float > > UndirectedGraph;

template <class Graph>
class DFSVisitor : public default_dfs_visitor {
public:
	void discover_vertex(typename Graph::vertex_descriptor u, const Graph& g) {
		std::cout << g[u].id << " ";
	}
};

template <class Graph>
class BFSVisitor : public default_bfs_visitor {
public:
	void discover_vertex(typename Graph::vertex_descriptor u, const Graph& g) {
		std::cout << g[u].id << " ";
	}
};

template <class Graph>
class MyGraph {
public:
	Graph graph;
	std::map<int, typename Graph::vertex_descriptor> hash;
	
	MyGraph();

	typename Graph::vertex_descriptor addVertex(int id);
	typename Graph::vertex_descriptor* getDescriptor(int id);
	bool addEdge(int idU, int idV, float weight);
	void loadEdges(std::string filename, int E);
	void removeVertex(int id);
	void removeEdge(int idU, int idV);
	void DFS(int id);
	void BFS(int id);
	void Prim(int id);
	void Kruskal();
	void Dijkstra(int source);
	void FloydWarshall();
};

template <class Graph>
MyGraph<Graph>::MyGraph() {
	Graph graph(0);
}

template <class Graph>
typename Graph::vertex_descriptor MyGraph<Graph>::addVertex(int id) {
	typename Graph::vertex_descriptor v = add_vertex(VData(id), graph);
	hash.insert(std::make_pair(id, v));
	return v;
}

template <class Graph>
typename Graph::vertex_descriptor* MyGraph<Graph>::getDescriptor(int id) {
	typename std::map<int, typename Graph::vertex_descriptor>::iterator it = hash.find(id);
	if (it != hash.end())
		return &it->second;
	return nullptr;
}

template <class Graph>
bool MyGraph<Graph>::addEdge(int idU, int idV, float weight) {
	typename Graph::vertex_descriptor* u = getDescriptor(idU);
	typename Graph::vertex_descriptor* v = getDescriptor(idV);
	if (u && v) {
		add_edge(*u, *v, weight, graph);
		return true;
	}
	return false;
}

template <class Graph>
void MyGraph<Graph>::loadEdges(std::string filename, int E) {
	int u, v;
	float weight;
	
	std::ifstream file(filename);
	if (file.is_open()) {
		for (int i = 0; i < E; ++i) {
			file >> u >> v >> weight;
			addEdge(u, v, weight);
		}
	}
}

template <class Graph>
void MyGraph<Graph>::removeVertex(int id) {
	typename Graph::vertex_descriptor* v = getDescriptor(id);
	if (v) {
		clear_vertex(*v, graph);
		remove_vertex(*v, graph);
		hash.erase(id);
	}
}

template <class Graph>
void MyGraph<Graph>::removeEdge(int idU, int idV) {
	typename Graph::vertex_descriptor* u = getDescriptor(idU);
	typename Graph::vertex_descriptor* v = getDescriptor(idV);
	if (u && v) {
		remove_edge(*u, *v, graph);
	}
}

template <class Graph>
void MyGraph<Graph>::DFS(int id) {
	DFSVisitor<Graph> vis;
	typename Graph::vertex_descriptor* v = getDescriptor(id);
	if (v) {
		std::cout << "Visited vertices:\n";
		depth_first_search(graph, visitor(vis).root_vertex(*v));
	}
	else
		std::cout << "Root vertex does not exist";
	std::cout << std::endl;
}

template <class Graph>
void MyGraph<Graph>::BFS(int id) {
	BFSVisitor<Graph> vis;
	typename Graph::vertex_descriptor* v = getDescriptor(id);
	if (v) {
		std::cout << "Visited vertices:\n";
		breadth_first_search(graph, *v, visitor(vis));
	}
	else
		std::cout << "Source vertex does not exist";
	std::cout << std::endl;
}

template <class Graph>
void MyGraph<Graph>::Prim(int id) {
	typename Graph::vertex_descriptor* v = getDescriptor(id);

	if (v) {
		std::vector<typename Graph::vertex_descriptor> p (num_vertices(graph));
		prim_minimum_spanning_tree(graph, &p[0], root_vertex(*v));

		std::cout << "Minimum spanning tree:" << std::endl;
		for (std::size_t i = 0; i < p.size(); ++i) {
			std::cout << "Parent of " << graph[vertex(i, graph)].id << " = ";
			if (p[i] != i)
			  std::cout << graph[vertex(p[i], graph)].id << std::endl;
			else
			  std::cout << "no parent" << std::endl;
		}
	}
	else
		std::cout << "Source vertex does not exist\n";
}

template <class Graph>
void MyGraph<Graph>::Kruskal() {
	typename property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);
	std::vector<typename Graph::edge_descriptor> treeEdges;

	kruskal_minimum_spanning_tree(graph, std::back_inserter(treeEdges));

	std::cout << "Minimum spanning tree:" << std::endl;
	for (typename std::vector<typename Graph::edge_descriptor>::iterator it = treeEdges.begin(); it != treeEdges.end(); ++it) {
		std::cout << "Edge between " << graph[source(*it, graph)].id << " and " << graph[target(*it, graph)].id;
		std::cout << ", weight = " << weight[*it] << std::endl;
	}
}

template <class Graph>
void MyGraph<Graph>::Dijkstra(int source) {
	typename Graph::vertex_descriptor* s = getDescriptor(source);
	if (s) {
		std::vector<typename Graph::vertex_descriptor> p (num_vertices(graph));
  		std::vector<int> d (num_vertices(graph));
		dijkstra_shortest_paths(graph, *s, predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, graph))).
								distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, graph))));

		std::cout << "Shortest paths from " << graph[*s].id << ":\n";
		typename Graph::vertex_iterator it, end;
		for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
			if (p[*it] == *it && *it != *s)
				std::cout << "Vertex " << graph[*it].id << " can't be reached\n";
			else {
    			std::cout << "Distance to " << graph[*it].id << " = " << d[*it] << ", ";
    			std::cout << "parent = " << graph[p[*it]].id << std::endl;
			}
  		}
	}
	else
		std::cout << "Source vertex does not exist";
	std::cout << std::endl;
}

template <class Graph>
void MyGraph<Graph>::FloydWarshall() {
	int V = num_vertices(graph);
	float** D = new float*[V];
	for (int i = 0; i < V; ++i)
		D[i] = new float[V];

  	floyd_warshall_all_pairs_shortest_paths(graph, D);

	std::cout << "All-pairs shortest paths\n";
	std::cout << "       ";
	typename Graph::vertex_iterator it, end;
	for (boost::tie(it, end) = vertices(graph); it != end; ++it)
		std::cout << std::setw(5) << graph[*it].id;
	std::cout << std::endl;

	boost::tie(it, end) = vertices(graph);
	for (int i = 0; i < V; ++i, ++it) {
		std::cout << std::setw(3) << graph[*it].id << " -> ";
		for (int j = 0; j < V; ++j) {
			if (D[i][j] == (std::numeric_limits<float>::max)())
				std::cout << std::setw(5) << "inf";
			else
				std::cout << std::setw(5) << D[i][j];
		}
		std::cout << std::endl;
	}

	for (int i = 0; i < V; ++i)
		delete [] D[i];
	delete [] D;
}

