include <Pair>




class Vertex {
    Public: 
    Vertex(int vertex_id, int ports_num) : vertex_id(vertex_id), ports_num(ports_num) {}
    
    void addPort()
    {
        ports_num++;
    }

    int vertex_id;
    int ports_num;
};

Typedef pair<Vertex, Port> vport;
Typedef pair<int, int> vport_id;


class Edge {
    pair<vport_id, vport_id> edge_id;
};

template<class T>
class PortGraph {
    private:
	// maps vport vp (vertix and port) to all the edges that has an vp as a source
	map<vport_id, set<Edge>> adjacency_list;
    // maps vport_id (vertix and port) to its vport
	map<vport_id, vport> vport_map;


    public:
    // Build empty PortGraph
    PortGraph();

    // Build PortGraph with 
    PortGraph(int n_vertices, vector<int> ports_num, vector<pair<pair<int, int>, pair<int, int>>> edges_list)
    {
        for (i = 0; i < n_vertices; ++i) {
            addVertix(i, ports_num[i]);
        }        
        for (auto pr : edges_list) {
            addEdge(pr);
        }
    }

    addVertix(int vertex_id, int ports_num) {
        Vertex v(vertex_id, ports_num);
        for (i = 0; i < ports_num; ++i) {
            adjacency_list[pair<int, int>(vertex_id, i)] = set();
            vport_map[pair<int, int>(vertex_id, i)] = vport(v, Port(i));
        }
    }

    addEdge(pair<vport_id, vport_id> p) {
        adjacency_list[p.first].add(Edge(p));
    }


    /*
    BFS
    DFS
    min_spanning_tree
    shortestpath(weight_function)
    transpose_graph()
    strongly_connected_components
    topological_sort
    max_flow
    min_cut()
    colouring
    induced_graph(vports/vertices/edges) // or filter
    is_subgraph(subgraph)
    netlistToPortGraph // maybe constructor
    findPathCost(vport, vport, cost_function)
    make_connected // maybe with min edges
    is_bipartite()
    is_reachable(vport source, vport dest)
    is_reachable(port src, port src)
    */
};

