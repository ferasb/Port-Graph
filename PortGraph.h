#include "Utilities.h"
using namespace std;

template <class P>
class Port {
    public:
    Port(int portId, P attr = P()) : port_id(portId) {}

    int getPortId()  { return port_id; }

    P getAttribute()  { return attribute; }

    bool operator==(Port p) { return p.port_id == port_id; }

    private:
    int port_id;
    P attribute;
};

template <class V, class P>
class Vertex {
    public: 
    Vertex(int vertex_id, int ports_num, V vertex_attr = V(), PortsAttributes ports_attr = PortsAttributes())
     : vertex_id(vertex_id), attribute(vertex_attr) {
        if (ports_attr.empty()) {
            for (int i = 0; i < ports_num; i++)
                ports.insert(Port<P>(i));
        }
        else {
            for (int i = 0; i < ports_num; i++)
                ports.insert(Port<P>(i, ports_attr[i]));
        }
    }
    
    void addPort(P attr = P()) { ports.insert(Port<P>(ports.size(), attr)); }

    int vertexId() { return vertex_id; }

    V getAttribute() { return attribute; }

    PortSet getPorts() { return ports; }

    int getPortsNum() { return ports.size(); }

    bool operator==(Vertex v) { return v.vertex_id == vertex_id; }

    void print() {
        ostringstream s;
        s << "Vertrex " << vertex_id << " with attribute " << attribute << ", with " << ports.size() << " ports." << endl;
        fprintf(stderr, s.str().c_str());
    }

    private:
    int vertex_id;
    V attribute;
    PortSet ports;
};


template <class V, class P, class E>
class Edge {
    public:
    Edge(vport src, vport dst, E attr = E()) 
    : source(src), dest(dst), edge_id(src.first.vertexId(), dst.first.vertexId()), attribute(attr) {}

    vport getSource() { return source; }

    vport getDest() { return dest; }

    edge_id getEdgeId() { return id; }

    bool operator==(Edge e) { return e.id == id; }

    E getAttribute() { return attribute; }

    // void print() {
    //     ostringstream s;
    //     s << "Edge (" << id.first << ", " << id.second << ") " <<  " with attribute " << attribute << "." << endl;
    //     fprintf(stderr, s.str().c_str());
    // }

    private:
    vport source;
    vport dest;
    edge_id id;
    E attribute;
};

template <class V = int, class P = int, class E = int>
class PortGraph {
    struct cmpVport {
        bool operator()(const vport_id& a, const vport_id& b) const {
            return a.first < b.first || (a.first == b.first && a.second < b.second);
        }
    };
    struct cmpEdge {
        bool operator()(const Edge<V, P, E>& a, const Edge<V, P, E>& b) const {
            return a.getEdgeId().first.first < b.getEdgeId().first.first;
        }
    };
    
    private:
	// maps vport vp (vertix and port) to all the edges that has an vp as a source
	map<vport_id, set<Edge<V, P, E>, cmpEdge>, cmpVport> adjacency_list;
    // maps vport_id (vertix and port) to its vport
	map<vport_id, vport, cmpVport> vport_map;


    public:
    // Build PortGraph with 
    PortGraph(int n_vertices, vector<int> ports_num, vector<edge_id> edges_list,
    VerticesAttributes verticesAttributes =  VerticesAttributes(),
    vector<PortsAttributes> portsAttributes = vector<PortsAttributes>(),
    EdgesAttributes edgesAttributes = EdgesAttributes())
    {
        for (int i = 0; i < n_vertices; ++i) {
            addVertix(i, ports_num[i], verticesAttributes[i], portsAttributes[i]);
        }
        for (int i = 0; i < n_vertices; ++i) {
            addEdge(edges_list[i], edgesAttributes[i]);
        }
    }

    void addVertix(int vertex_id, int ports_num, V attr = V(), PortsAttributes ports_attr = PortsAttributes())
    {
        Vertex<V, P> v(vertex_id, ports_num, attr, ports_attr);
        for (int i = 0; i < ports_num; ++i) {
            adjacency_list[pair<int, int>(vertex_id, i)] = set<Edge<E> >();
            vport_map[pair<int, int>(vertex_id, i)] 
            = vport(Vertex<V, P>(vertex_id, ports_num, attr, ports_attr), Port<P>(i, ports_attr[i]));
        }
    }

    void addEdge(edge_id id, E attr = E()) 
    {
        adjacency_list[id.first].insert(Edge<V, P, E>(vport_map[id.first], vport_map[id.second], attr));
    }

    void Print() 
    {
        for (auto x : vport_map) {
            Vertex<V, P>& v = x.second.first;
            v.print();
        }
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

