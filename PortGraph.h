#include "Utilities.h"
using namespace std;

template<class P>
class Port {
    public:
    Port(int portId, const P& attr = P()) : port_id(portId) {}

    int getPortId()  { return port_id; }

    P getAttribute()  { return attribute; }

    bool operator==(Port p) { return p.port_id == port_id; }

    private:
    int port_id;
    P attribute;
};

template<class V, class P>
class Vertex {
    public: 
    Vertex(int vertex_id, int ports_num, V vertex_attr = V(), PortsAttributes ports_attr = PortsAttributes(ports_num))
     : vertex_id(vertex_id), attribute(vertex_attr) {
        if (ports_attr.empty()) {
            for (int i = 0; i < ports_num)
                ports.insert(Port(i));
        }
        else {
            for (int i = 0; i < ports_num)
                ports.insert(Port(i), ports_attr[i]);
        }
    }
    
    void addPort(P attr = P())
    {
        ports.insert(Port(ports.size(), attr))
    }

    int vertexId() { return vertex_id; }

    P getAttribute() { return attribute; }

    bool operator==(Vertex v) { return v.vertex_id == vertex_id; }

    private:
    int vertex_id;
    V attribute;
    PortSet ports;
};


template<class V, class P, class E>
class Edge {
    public:
    Edge(const vport& src, const vport& dst, E attr = E()) 
    : source(src), dest(dst), edge_id(src.fst.vertexId(), dst.fst.vertexId()), attribute(attr) {}

    vport getSource() { return source; }

    vport getDest() { return dest; }

    edge_id getEdgeId() { return id; }

    bool operator==(Edge e) { return e.id == id; }

    E getAttribute() { return attribute; }

    private:
    vport source;
    vport dest
    edge_id id;
    E attribute;
};

template<class V, class P, class E>
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

