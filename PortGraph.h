//
// Created by Mario Barbara and Feras Bisharat
//

#ifndef PORJECT_OFF_PORTGRAPH_H
#define PORJECT_OFF_PORTGRAPH_H

#include "Utilities.h"
using namespace std;

/**************** Port class *****************/
/**
    includes :
     int port_id - the id for the port which is unique per port
     P attribute - the attribute for the port of type "P"
    Port Class requires the class P to have :
     1) default constructor
     2) copy constructor
    Port Class methods :
     1) getPortId : gets port id
     2) getAttribute : gets port attribute
     3) print : prints the info of the Port
 **/
template <class P>
class Port {
public:
    /* Port Class constructor , create new port with 'portId' id and 'attr' attribute */
    explicit Port(int portId = -1, P attr = P()) : port_id(portId) , attribute(attr) {}

    /* returns the id of 'this' Port */
    int getPortId() { return port_id; }

    /*returns the attribute of 'this' Port */
    P getAttribute()  { return attribute; }

    /* Comparison operator of 'Port' Class
      - two 'Port' are equal if and only if there id's are the same
    */
    bool operator==(const Port& that) { return that.getPortId() == port_id; }

    /* Smaller operator of 'Port' Class
      - this 'Port' is smaller if the only if 'this' id is smaller the 'that' id
    */
    bool operator<(const Port& that) {
        return this->port_id < that.getPortId();
    }

    /* assign new attribute value for 'this' port */
    void setAttribute(P& attr) {
        attribute = attr;
    }

    /* print method for Port Class , prints fields in the following format
        Port <Port id> with attribute <Port's attribute>
    */
    void print() {
        ostringstream s;
        s << "Port " << port_id << " with attribute " << attribute << endl;
        fprintf(stderr, s.str().c_str());
    }

private:
    int port_id;
    P attribute;
};
/**********************************************/

/**************** Vertex class ****************/
/**
    includes :
      1) int vertex_i - the id for the vertex which is unique per vertex
      2) V attribute - the attribute for the vertex of type "V"
      3) int n_outgoing_edges - number of out going edges for this vertex , dosen't count multi-edges between vertecies
      4) int n_ingoing_edges - number of in going edges to this vertex , dosen't count multi-edges between vertecies
      5) PortMap ports - maping from port id to port Class that has the same id,the ports aren't unique per vertex and can be shared
    Vertex Class requires the class P to have :
      1) default constructor
      2) copy constructor
    Vertex Class requires the class V to have :
      1) default constructor
      2) copy constructor
    Vertex Class methods :
      1) incOut/Ingoing edges
      2) decOut/Ingoing edges
      3) vertexId
      4) getAttribute
      5) getPorts
      6) getPortsNum
      7) getPort
      8) addPortAttr
      9) print
 **/
template <class V, class P>
#define PortMap map<int, Port<P>>
class Vertex {
public:

    /* Vertex Class defualt constructor only init  n_ingoing_edges and  n_outgoing_edges to zero */
    Vertex () {
        n_ingoing_edges = 0;
        n_outgoing_edges = 0;
    }

    /* Vertex Class constuctor , creates new vertex
      with vertex id and vertex_attr , init ports_num 'ports' with ascending id starting form 0 , init the ith port
      with the  attribute ith attribute in ports_attr array if present,
      else uses the defualt constructor and maps them accordingly
    */
    Vertex(int vertex_id , int ports_num , V vertex_attr = V(), PortsAttributes ports_attr = PortsAttributes())
            : vertex_id(vertex_id), attribute(vertex_attr) {
        if (ports_attr.empty()) {
            for (int i = 0; i < ports_num; i++)
                ports.insert(pair<int,Port<P>>(i,Port<P>(i)));
        }
        else {
            for (int i = 0; i < ports_num; i++)
                ports.insert(pair<int,Port<P>>(i, Port<P>(i,ports_attr[i])));
        }
        // new version
        n_ingoing_edges = 0;
        n_outgoing_edges = 0;
    }

    /* increasing the number of out going edges only when adding a new out going edge to 'this' vertex
       - NOT including multi edges
    */
    void incOutgoingEdges() {
        n_outgoing_edges++;
    }

    /* increasing the number of in going edges only when adding a new in going edge to 'this' vertex
     - NOT including multi edges
    */
    void incIngoingEdges() {
        n_ingoing_edges++;
    }

    /* decreasing the number of out going edges only when removing the last out going ede to 'this' edge*/
    void decOutgoingEdges() {
        if(n_outgoing_edges >0)
            n_outgoing_edges--;
    }

    /* decreasing the number of in going edges only when removing the last in going ede to 'this' edge*/
    void decIngoingEdges() {
        if(n_ingoing_edges >0)
            n_ingoing_edges--;
    }

    /* returns the number of out going edges from 'this' vertex */
    int outGoingEdges() { return n_outgoing_edges;}

    /* returns the number of in going edges to 'this' vertex */
    int inGoingEdges() { return n_ingoing_edges;}

    /* retruns the id of 'this' vertex :- vectex_id */
    int vertexId() { return vertex_id; }

    /* retruns the attribute of 'this' vertex :- attribute */
    V getAttribute() { return attribute; }

    /* return a map of type <int,Port> that maps the ports id to the corresponding Port Class */
    PortMap getPorts() { return ports; }

    /* returns the number of ports that 'this' vertex has */
    int getPortsNum() { return ports.size(); }

    /* Comparison operator of 'Vertex' Class
       - two 'Vetex' are equal if and only if there id's are the same
    */
    bool operator==(Vertex that_v) { return that_v.vertexId() == vertex_id; }

    /* Samller operator of 'Vertex' Class
       - this 'Vertex' is smaller if the only if 'this' id is smaller the 'that' id
    */
    bool operator<(Vertex that_v) {
        return  that_v.vertexId() > vertex_id;
    }

    /* returns 'port' Class with port_id id if present else error */
    Port<P> getPort(int port_id) {
        //check -
        return ports[port_id];
    }

    /* assgin 'this' vertex a new vertex_id with id */
    void setVertexId(int id) {
        vertex_id = id;
    }

    /* assgin 'this' vertex a new attribute with attr */
    void setAttribute(V attr) {
        attribute = attr;
    }

    /* adds a new port with port_id id and attr attribute to 'this' vertex if not present
       else updates a new attribute to the corresponding port with "port_id" id
    */
    void addPortWithAttr(int port_id,P attr) {
        Port<P> p = Port<P>(port_id, attr);
        ports[port_id] = p;
    }

    /* print method for Vertex Class , prints fields in the following format :
       vertex <Vertex_id> with attribute <Vertex's attribute> with <vertex's port number>
       and Foreach port "in" this vertex prints :
       Port <Port id> with attribute <Port's attribute>
    */
    void print() {
        ostringstream s;
        s << "Vertex " << vertex_id << " with attribute " << attribute << ", with " << ports.size() << " ports." << endl;
        s << "Vertex " << vertex_id << " with Ports  : " << endl ;
        fprintf(stderr, s.str().c_str());
        for (auto p : ports)
            p.second.print();
    }

    /* removes the port with "port_id" id if present Else no changes happens */
    void removePort(int port_id) {
        if(ports.find(port_id) == ports.end())
            return;
        ports.erase(ports[port_id]);
    }

private:
    int vertex_id;
    V attribute;
    int n_outgoing_edges;
    int n_ingoing_edges;
    PortMap ports;
};
/*********************************************/


/**************** Edge class ****************/
/**
    includes :
      1)vport source - the sourse of 'this edge
      2)vport dest - the destination of 'this' edge
      3)edge_id id - the id of 'this' edge
      4)E attribute - the attribute of 'this' edge
    Edge Class requires the class P to have :
      1) default constructor
      2) copy constructor
    Edge Class requires the class V to have :
      1) default constructor
      2) copy constructor
    Edge Class requires the class E to have :
      1) default constructor
      2) copy constructor
    Edge Class methods :
      1) getSource
      2) getDest
      3) EdgeId
      4) getAttribute
      5) printIds
 **/
template <class V, class P, class E>
class Edge {
public:

    /*Edge Class defualt constructor*/
    Edge() = default ;

    /*Edge Class constructor , init new directed 'edge' from src to dst with optional attribute att */
    Edge(vport src , vport dst , E attr = E()): source(src), dest(dst), attribute(attr) {
        vport_id p1 = std::make_pair(src.first.vertexId(),src.second.getPortId());
        vport_id p2 = std::make_pair(dst.first.vertexId(),dst.second.getPortId());
        id = std::make_pair(p1, p2) ;
    }

    /* return the vport source of 'this' edge - vertex and port*/
    vport getSource() { return source; }

    /* return the vport destination of 'this' edge - vertex and port*/
    vport getDest() { return dest; }

    /* returns the edge id*/
    edge_id EdgeId() { return id; }

    /* Comparison operator of 'Edge' Class
       - two 'Port' are equal if and only if there id's are the same
    */
    bool operator==(Edge e) { return e.id == id; }

    /* return the attribute of 'this' edge*/
    E getAttribute() { return attribute; }

protected:
    ostringstream toString() {
        ostringstream s;
        s << "Edge (" << id.first.first << ", " << id.first.second << ") " << "-- (" << id.second.first << ", " << id.second.second << ") " ;
        s << "with attribute " << attribute << endl ;
        return s;
    }

    ostringstream toStringIds() {
        ostringstream s;
        s << "Edge (" << id.first.first << ", " << id.first.second << ") " << "-- (" << id.second.first << ", " << id.second.second << ") " ;
        return s;
    }

public:
    /* print method for Edge Class , prints fields in the following format :
       Edge <source vport> <destination vport> with attribute <Edge's attribute>
    */
    void print() {
        ostringstream s = toString();
        fprintf(stdout, s.str().c_str());
    }

    void printIds() {
        ostringstream s = toStringIds();
        fprintf(stdout, s.str().c_str());
    }

private:
    vport source;
    vport dest;
    edge_id id;
    E attribute;
};
/**********************************************/

/** vport is a pair of <Vertex Class> and <Port Class> **/
/** this vport is smaller than that vport only and only if
    this vertex id is smaller than that vertex id OR
    both vertecies are EQUAL and this port id is smaller than that port id
**/
struct cmpVport {
    bool operator()(const vport_id& a,const vport_id& b) const  {
        return (a.first < b.first) || (a.first == b.first && a.second < b.second);
    }
};
/**********************************************/

/** this edge is smaller than that edge only and only if
    this edge id is smaller than that edge id
**/
template <class V = int, class P = int, class E = int>
struct cmpEdge {
    bool operator()(const Edge<V, P, E>& a,const Edge<V, P, E>& b) const {
        return a.EdgeId() < b.EdgeId() ;
    }
};
/*************************************************/

/* Forward declarations */
template <class V , class P , class E >
class DFSIterator;

template <class V , class P , class E >
class BFSIterator;

template <class V , class P , class E >
class DFSVertexIterator;

template <class V , class P , class E >
class BFSVertexIterator;


/**************************************************/

/**************** Port Graph class ****************/
/**
    includes :
      1) map adjacency_list for vports and it's reverse
      2) bool is_transpose : true if the graph is reversed Else false
      3) map vertex_map : maps vport id to vports
      4) map vertex_neighbors for vertecies and it's reverse
      4) map shortest_paths and shortest_paths_weights : maps the path from 'src' to 'dst' and it's weight
    PortGraph Class requires the class P to have :
      1) default constructor
      2) copy constructor
    PortGraph Class requires the class V to have :
      1) default constructor
      2) copy constructor
    PortGraph Class requires the class E to have :
      1) default constructor
      2) copy constructor
 **/
template <class V = int, class P = int, class E = int>
class PortGraph {
private:

    // maps vport vp (vertix and port) to all the edges that has an vp as a source
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adjacency_list;
    // maps vport vp (vertix and port) to all the reverse edges that has an vp as a source
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> backwards_adjacency_list;

    // flag for the reverse Graph
    bool is_reversed ;

    // maps vport_id (vertix and port) to its vport
    // vertex id -> vertex
    // maps vertex_id (int) to its vertex
    map<vertex_id, Vertex<V,P>> vertex_map;

    //maps vertex_id -> neighbor vertices
    map<int,map<int,bool>> vertex_neighbors;
    //maps vertex_id -> neighbors vertices of the reverse graph
    map<int,map<int,bool>> reverse_vertex_neighbors;

    // params for SCC
    const int UNSEEN = -1;
    const int SEEN = 1;

    // for Algo's in case of reverse
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& AdjacencyList() {
        return is_reversed ? backwards_adjacency_list : adjacency_list ;
    }

    // for Algo's in case of reverse
    map<int,map<int,bool>>& VertexNeighbors() {
        return is_reversed ? reverse_vertex_neighbors : vertex_neighbors ;
    }

    // for shortest paths
    typedef vector<edge_id> Path;
    map<pair<vport_id, vport_id>, double> shortest_paths_weights;
    map<pair<vport_id, vport_id>, Path> shortest_paths;
    const vport_id  END_VPORT = vport_id(-1,-1);
    const int  END_VERTEX = -1;

public:

    /* Port graph defualt constructor - empty graph */
    PortGraph() = default;

    /* create new Port Graph and consder it reversed or not if 'isReversed' false or Not*/
    PortGraph(bool isReversed) :is_reversed(isReversed) {}

    /* init new port graph with n_vertices vertices which the ith vertex has verticesAttributes[ith] attribute
       and ports_num[ith] ports ,each kth port has portsAttributes[ith][kth] attribute .
       Foreach ith_edge in edges_list :
       init new edges with edgesAttributes[ith_edge] attribute
    */
    PortGraph(int n_vertices, vector<int> ports_num, vector<edge_id> edges_list,
              VerticesAttributes verticesAttributes =  VerticesAttributes(),
              vector<PortsAttributes> portsAttributes = vector<PortsAttributes>(),
              EdgesAttributes edgesAttributes = EdgesAttributes())
    {
        is_reversed = false;
        for (int i = 0; i < n_vertices; ++i) {
            if(verticesAttributes.empty() && portsAttributes .empty())
                addVertex(i, ports_num[i]);
            else if(portsAttributes .empty())
                addVertex(i, ports_num[i],verticesAttributes[i]);
            else if(verticesAttributes.empty())
                addVertex(i, ports_num[i],V(),portsAttributes[i]);
            else {
                addVertex(i, ports_num[i],verticesAttributes[i],portsAttributes[i]);
            }
        }

        //add reverse edge
        for (int i = 0; i < edges_list.size(); ++i) {
            if(edgesAttributes.empty())
                addEdge(edges_list[i]);
            else {
                addEdge(edges_list[i], edgesAttributes[i]);
            }
        }

    }

    /* reverse every directed edge and
       reverse each vertex neighbor
    */
    void reverseGraph() {
        //check
        is_reversed = !is_reversed;
        auto tmp1 = adjacency_list;
        adjacency_list = backwards_adjacency_list;
        backwards_adjacency_list = tmp1;
        auto tmp2 = vertex_neighbors;
        vertex_neighbors = backwards_adjacency_list;
        backwards_adjacency_list = tmp2;
    }

    /* True if the graph is reversed Else False */
    bool isReversed() {
        return is_reversed;
    }

    /* add new vertex with 'vertex_id' id and 'ports_num' port with 'attr' attribute
       the ith port has the ports_attr[ith] attribute if present
       Else the port has the defualt value
    */
    void addVertex(int vertex_id, int ports_num, V attr = V(), PortsAttributes ports_attr = PortsAttributes()) {
        for (int i = 0; i < ports_num; ++i) {
            adjacency_list[vport_id(vertex_id, i)] = set<Edge<V, P, E>, cmpEdge<V,P,E>>();
            // reverse Graph
            backwards_adjacency_list[vport_id(vertex_id, i)] = set<Edge<V, P, E>, cmpEdge<V,P,E>>();
            P portAttr =  ports_attr.empty() ? P() : ports_attr[i] ;
        }
        vertex_map[vertex_id] = Vertex<V, P>(vertex_id, ports_num, attr, ports_attr);
        if (vertex_neighbors.find(vertex_id) == vertex_neighbors.end())
            vertex_neighbors[vertex_id] = map<int,bool>();
        reverse_vertex_neighbors[vertex_id] = map<int,bool>();
    }

    /* add new edge with id == <src,dst> with 'attr' attibute if Not present
       Else updates the Edge attribute with 'attr'
    */
    void addEdge(edge_id id, E attr = E()) {
        // check -
        // add edge
        Vertex<V,P> v1 = vertex_map[id.first.first];
        Vertex<V,P> v2 = vertex_map[id.second.first];
        vport vp1 = vport(v1,v1.getPort(id.first.second));
        vport vp2 = vport(v2,v2.getPort(id.second.second));
        // check
        adjacency_list[id.first].insert(Edge<V, P, E>(vp1, vp2, attr));
        // add reverse Edge
        backwards_adjacency_list[id.second].insert(Edge<V, P, E>(vp2, vp1, attr));
        //check
        if (vertex_neighbors[id.first.first].find(id.second.first) == vertex_neighbors[id.first.first].end()
           && id.first.first != id.second.first) {
            if (reverse_vertex_neighbors[id.second.first].find(id.first.first) == reverse_vertex_neighbors[id.second.first].end()
               && id.first.first != id.second.first) {
                reverse_vertex_neighbors[id.second.first].insert(pair<int, bool>(id.first.first, true));
            }
            vertex_map[id.first.first].incOutgoingEdges();
            vertex_map[id.second.first].incIngoingEdges();
            vertex_neighbors[id.first.first].insert(pair<int, bool>(id.second.first, true));
        }
    }

    /*Print the vertices and ports per vertex */
    void Print() {
        for (auto x : vertex_map) {
            x.second.print();
        }
    }

    /*Print the edges*/
    void PrintEdges() {
        auto adj_list = AdjacencyList();
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it)
            for(auto e : (*it).second)
                e.print();
    }

    /*Print the edges without attributes*/
    void PrintEdgesIds() {
        auto adj_list = AdjacencyList();
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it)
            for(auto e : (*it).second)
                e.printIds();
    }

    /* return the outGoing edges id's form the vport 'id' if present
       Else error
    */
    vector<edge_id> getOutgoingEdges(vport_id id) {
        auto adj_list = AdjacencyList();
        vector<edge_id> to_return;
        for(auto e : adj_list[id]) {
            to_return.push_back(e.EdgeId());
        }
        return to_return;
    }

/********** Topological Sort **********/

// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: is g a DAG (return value), a topological ordering of g (order).
// comment: order is valid only if g is a DAG.
    vector<vport_id> topological_sort() {
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adj_list = AdjacencyList();
        // compute indegree of all nodes
        map<vport_id,int> vport_indegree = map<vport_id,int>();
        for(auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            set<Edge<V, P, E>,cmpEdge<V,P,E>> edge_set = (*it).second;
            vport_id id = (*it).first;
            for(auto edge : edge_set) {
                vport_indegree[edge.EdgeId().second]++ ;
            }
        }

        // order sources first
        vector<vport_id> vport_order = vector<vport_id>();
        for(auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if (vport_indegree[id] == 0)
                vport_order.push_back(id);
        }
        // go over the ordered nodes and remove outgoing edges,
        // add new sources to the ordering
        for(auto id : vport_order) {
            for (auto edge : adj_list[id]) {
                vport_id neighbor = edge.EdgeId().second;
                vport_indegree[neighbor] -- ;
                if (vport_indegree[neighbor] == 0)
                    vport_order.push_back(neighbor);
            }
        }
        // check -
        return vport_order;
    }

/********** Strongly Connected Components **********/

protected:
    void KosarajuDFS(const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list
            ,vport_id id, map<vport_id,vport_id>& S, map<vport_id,int>& colorMap, int color) {

        colorMap[id] = color;
        for(auto edge : adj_list[id]){
            vport_id neighbor = edge.EdgeId().second;
            if(colorMap[neighbor] == UNSEEN)
                KosarajuDFS(adj_list,neighbor, S,colorMap,color);
        }
        S.insert({id,id});
    }

// Compute the number of SCCs and maps nodes to their corresponding SCCs.
// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: the number of SCCs (return value), a mapping from node to SCC color (components).
    int findSCC(map<vport_id ,int>& components) {
        // first pass: record the `post-order' of original graph
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adj_list = AdjacencyList();
        map<vport_id ,vport_id> postOrder , dummy ;
        map<vport_id ,int>  seen ;
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it){
            vport_id id = (*it).first;
            seen[id] = UNSEEN;
            components[id] = UNSEEN;
        }
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if(seen[id] == UNSEEN)
                KosarajuDFS(adj_list , id,postOrder,seen,SEEN);
        }
        // reverse iterator
        // second pass: explore the SCCs based on first pass result
        int numSCC = 0;
        for (auto rit = adj_list.rbegin(); rit != adj_list.rend(); ++rit) {
            vport_id id = (*rit).first;
            reverseGraph();
            map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> rev_adj_list = AdjacencyList();
            reverseGraph();
            if (components[postOrder[id]] == UNSEEN)
                KosarajuDFS(rev_adj_list, postOrder[id], dummy, components, numSCC++);
        }
        return numSCC;
    }

public:
// Computes the SCC graph of a given digraph.
// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: strongly connected components graph of g (sccg).
    void findSCCgraph(vector<set<int>>& sccg) {
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adj_list = AdjacencyList();
        map<vport_id, int> component;
        int n = findSCC(component);
    }

/********** Min Spanning Tree **********/
private:
    /* Union Find Class where every 'group' id is a 'vport' id
       and two groups are merged if the "leaders" of both groups are neighbors
    */
    struct unionfindPG {
        map<vport_id ,int> rank;
        map<vport_id ,vport_id >parent;

        explicit unionfindPG (map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list) {
            for (auto it = adj_list.begin(); it != adj_list.end(); ++it){
                vport_id id = (*it).first;
                rank[id] = 0;
                parent[id] = id;
            }

        }

        vport_id find(vport_id x) {
            vport_id tmp = x;
            while(x!=parent[x]) x=parent[x];
            while(tmp!=x) {
                vport_id remember=parent[tmp];
                parent[tmp]=x;
                tmp=remember;
            }
            return x;
        }

        void unite(vport_id p, vport_id q) {
            p = find(p);
            q = find(q);
            if(q==p) return;
            if(rank[p] < rank[q]) parent[p] = q;
            else parent[q] = p;
            if(rank[p] == rank[q]) rank[p]++;
        }
    };

    /*find MST for the current port graph according to the 'weightFunc' function
      the function uses kruskal algorthim to find mst where each node in the graph from a 'vport' type
      and two node are neihbors if there is a directed edge between them
      comment : weight function input include the id of the edgs and it's attribute
    */
public:
    pair<double,PortGraph<V,P,E>> findMST(double(*weightFunc)(edge_id ,E attr)) {
        int n = AdjacencyList().size();
        PortGraph<V,P,E> mst = PortGraph<V,P,E>();
        // (weight , edge_id)
        vector<pair<double,edge_id>> edges;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it){
            set<Edge<V, P, E>,cmpEdge<V,P,E>> s_edge = (*it).second;
            for(auto e : s_edge)
                edges.push_back(pair<double,edge_id>(weightFunc(e.EdgeId(),e.getAttribute()) ,e.EdgeId()));
        }
        // with weightFunc
        sort(edges.begin(), edges.end());
        double mst_cost = 0;
        unionfindPG components(adj_list);
        for (pair<double,edge_id>& e : edges) {
            if (components.find(e.second.first) != components.find(e.second.second)) {
                // W function on edges
                edge_id id = e.second;
                mst.addVport(id.first,vertex_map[id.first.first].getAttribute(),vertex_map[id.first.first].getPort(id.first.second).getAttribute());
                mst.addVport(id.second,vertex_map[id.second.first].getAttribute(),vertex_map[id.second.first].getPort(id.second.second).getAttribute());
                mst.addEdge(id,getEdge(id).getAttribute());
                mst_cost += e.first;
                components.unite(id.first, id.second);
            }
        }
        return pair<double,PortGraph<V,P,E>>(mst_cost,mst);
    }

/*************** BFS/DFS ***************/

    /* used to indicate the end of 'vports' in pg and ther is no more */
    PGVportIterator vportEnd() {
        return PGVportIterator(END_VPORT);
    }

    /* used to indicate the end of 'vertices' in pg the there is no more */
    PGVertexIterator vertexEnd() {
        return PGVertexIterator(-1);
    }

    /* return a vector of the 'vports' ids in the graph */
    vector<vport_id> getVports() {
        vector<vport_id> res;
        for(auto ver : vertex_map){
            for(auto p : ver.second.getPorts()) {
                vport_id id = vport_id(ver.second.vertexId(), p.second.getPortId());
                res.push_back(id);
            }
        }
        return res;
    }

    /* return a vector of the 'vertices' ids in the graph */
    vector<int> getVertices() {
        vector<int> res;
        for(auto it : vertex_map)
            res.push_back(it.first);
        return res;
    }

    /* return a maping between vports id and the set of the out going edges from each vport */
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> getAdjList() {
        return AdjacencyList();
    }

    /* return a set of the out going edges "id's" vport */
    set<Edge<V,P,E>,cmpEdge<V,P,E>> getVportAdjList(vport_id id) {
        return AdjacencyList()[id];
    }

    /* return a vector of vertices ids that has in going edge from vertex with 'src' id
       two vertices are considered neighbores if they have one edge between them at least
    */
    vector<int> getVertexAdjList(int src) {
        vector<int> res;
        for(auto it : VertexNeighbors()[src])
            res.push_back(it.first);
        return res;
    }

/************** Bipartite **************/

    /* returns True if there a directed edge from 'src' to 'dst' Else False */
    bool isNeighbers(vport_id& src, vport_id& dst) {
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for(const Edge<V,P,E>& e : adj_list[src]){
            if(e.EdgeId().second == dst)
                return true;
        }
        return false;
    }

    /* return True if there is portA in vertex 'src' and port2 in vertex 'dst'
       such that there is a directed edge from < src,portA > to < dst,portB >
       Else False
    */
    bool isNeighbers(int src, int dst) {
        const PortMap mp_src = vertex_map[src].getPorts();
        for (auto it1 = mp_src.begin(); it1 != mp_src.end();it1++) {
            //Port<P> p1 = (*it1).second;
            vport_id id1 = vport_id(src,(*it1).second.getPortId());
            const PortMap mp_dst = vertex_map[dst].getPorts();
            for (auto it2 = mp_dst.begin(); it2 != mp_dst.end();it2++) {
                //Port<P> p2 = (*it2).second;
                vport_id id2 = vport_id(dst,(*it2).second.getPortId());
                if (isNeighbers(id1, id2))
                    return true;
            }
        }
        return false;
    }

private:
    // This function returns true if
    // port graph is Bipartite, else false
    bool isBipartiteUtil(vport_id src,map<vport_id,int>& colorArr) {
        colorArr[src] = 1;

        // Create a queue (FIFO) of vertex numbers a
        // nd enqueue source vertex for BFS traversal
        queue<vport_id> q;
        q.push(src);

        // Run while there are vertices in queue (Similar to
        // BFS)
        while (!q.empty()) {
            // Dequeue a vertex from queue ( Refer
            // http://goo.gl/35oz8 )
            vport_id u = q.front();
            q.pop();

            // Return false if there is a self-loop
            if (isNeighbers(u,u))
                return false;

            // Find all non-colored adjacent vertices
            for (auto v : getVports()) {
                // An edge from u to v exists and
                // destination v is not colored
                if (isNeighbers(u,v) && colorArr[v] == -1) {
                    // Assign alternate color to this
                    // adjacent v of u
                    colorArr[v] = 1 - colorArr[u];
                    q.push(v);
                }

                    // An edge from u to v exists and destination
                    // v is colored with same color as u
                else if (isNeighbers(u,v) && colorArr[v] == colorArr[u])
                    return false;
            }
        }

        // If we reach here, then all adjacent vertices can
        // be colored with alternate color
        return true;
    }

public:
    // Returns true if port garph is Bipartite, else false
    bool isBipartite() {
        // Create a color array to store colors assigned to all
        // veritces. Vertex/ number is used as index in this
        // array. The value '-1' of colorArr[i] is used to
        // ndicate that no color is assigned to vertex 'i'.
        // The value 1 is used to indicate first color is
        // assigned and value 0 indicates second color is
        // assigned.
        // if src not in graph ...
        map<vport_id,int> colorArr;

        for (auto id : getVports()) {
            colorArr.insert(pair<vport_id, int>(id,-1));
        }
        // This code is to handle disconnected graoh
        for (auto id : getVports())
            if (colorArr[id] == -1)
                if (!isBipartiteUtil(id,colorArr))
                    return false;
        return true;
    }

/********** Induced Graph **********/

    /* return the edge with 'id' id if present Else error*/
    Edge<V,P,E> getEdge(edge_id id) {
        // check -
        assert(isNeighbers(id.first,id.second));
        auto adj_list = this->AdjacencyList();
        for(auto e : adj_list[id.first]){
            if(e.EdgeId() == id)
                return e;
        }
        //won't get here
        return Edge<V,P,E>();
    }

    /* returns the vertex's attribute with id 'vertex_i'*/
    V getVertexAttr(int vertex_i) {
        return vertex_map[vertex_i].getAttribute();
    }

    /* returns the Port's attribute with 'port_id' id in 'vertex_id' id
       if present Else error
    */
    P getPortAttr(int vertex_id , int port_id) {
        return vertex_map[vertex_id].getPort(port_id).getAttribute();
    }

    /* adds a new vport with 'vertex_attr' vertex atttibute and
      'port_attr' port attribute to the graph , if the vertex isn't present
       creates new one and adds a new port
    */
    void addVport(vport_id id,V vertex_attr,P port_attr) {
        //check -
        if(vertex_map.find(id.first) == vertex_map.end()) {
            Vertex<V, P> v = Vertex<V, P>();
            v.setVertexId(id.first);
            v.setAttribute(vertex_attr);
            v.addPortWithAttr(id.second, port_attr);
            vertex_map[id.first] = v;
            return;
        }
        addPort(id,port_attr);
    }

    /* adds a new port with 'id.second' id and 'attr' attribute if not present
       else if vertex is peresnt updates the 'id.second' port in 'id.first' vertex
       with attr attribute Else error
    */
    void addPort(vport_id id,P attr = P()) {
        //check -
        vertex_map[id.first].addPortWithAttr(id.second,attr);
    }

    /* returns all the edges form vertex 'src' to vertex 'dst' */
    vector<edge_id> getNeighbers(int src, int dst) {
        vector<edge_id> res;
        const PortMap mp_src = vertex_map[src].getPorts();
        for (auto it1 = mp_src.begin(); it1 != mp_src.end();it1++) {
            //Port<P> p1 = (*it1).second;
            vport_id id1 = vport_id(src,(*it1).second.getPortId());
            const PortMap mp_dst = vertex_map[dst].getPorts();
            for (auto it2 = mp_dst.begin(); it2 != mp_dst.end();it2++) {
                //Port<P> p2 = (*it2).second;
                vport_id id2 = vport_id(dst,(*it2).second.getPortId());
                if (isNeighbers(id1, id2))
                    res.push_back(edge_id(id1, id2));
            }
        }
        return res;
    }

    ///  Induced Graph by vport sub set
    PortGraph inducedGraph(bool(*pred)(vport_id, V)) {
        map<vport_id,int> sub_vport_set;
        // new induced graph
        PortGraph<V,P,E> pg = PortGraph<V,P,E>(this->isReversed());
        map<int,int> hash;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        // get all vports that satisfies Pred
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if(pred(id,getPortAttr(id.first,id.second))){
                sub_vport_set.insert(pair<vport_id,int>(id,1));
            }
        }
        // get all the edges for the induced graph
        for (auto it1 = sub_vport_set.begin();it1 != sub_vport_set.end();it1++) {
            //check
            for (auto it2 = sub_vport_set.begin(); it2 != sub_vport_set.end();it2++) {
                if (it1 == it2)
                    continue;
                vport_id src = (*it1).first;
                vport_id dst = (*it2).first;
                // add new vport
                if (isNeighbers(src,dst)) {
                    // new
                    pg.addVport(src,vertex_map[src.first].getAttribute(),vertex_map[src.first].getPort(src.second).getAttribute());
                    pg.addVport(dst,vertex_map[dst.first].getAttribute(),vertex_map[dst.first].getPort(dst.second).getAttribute());
                    // new edge
                    edge_id e_id = edge_id(src,dst);
                    E e_attr = getEdge(e_id).getAttribute();
                    pg.addEdge(e_id,e_attr);
                }
            }
        }

        // new induced port graph
        return pg;
    }

    ///  Induced Graph by vertex sub set
    PortGraph inducedGraph(bool(*pred)(int, P)) {
        map<int,int> sub_vertex_set;
        // new induced graph
        PortGraph<V,P,E> pg = PortGraph<V,P,E>(this->isReversed());
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        // get all vports that satisfies Pred
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if (pred(id.first,vertex_map[id.first].getAttribute())) {
                sub_vertex_set.insert(pair<int,int>(id.first,1));
            }
        }
        // get all the edges for the induced graph
        for (auto it1 = sub_vertex_set.begin();it1 !=sub_vertex_set.end();it1++) {
            int src = (*it1).first;
            // add first vertex
            for (auto it : vertex_map[src].getPorts()) {
                auto p = it.second;
                pg.addVport(vport_id(src,p.getPortId()), vertex_map[src].getAttribute(),p.getAttribute());
            }

            for (auto it2 =  sub_vertex_set.begin(); it2 != sub_vertex_set.end();it2++) {
                if (it2 == it1)
                    continue;
                int dst = (*it2).first;
                vector<edge_id> vec = getNeighbers(src, dst);
                // add second vertex
                for (auto it : vertex_map[dst].getPorts()) {
                    auto p = it.second;
                    pg.addVport(vport_id(dst, p.getPortId()), vertex_map[dst].getAttribute(), p.getAttribute());
                }

                if (!vec.empty()) {
                    //new edge
                    for (auto e_id : vec) {
                        pg.addPort(e_id.first, vertex_map[src].getPort(e_id.first.second).getAttribute());
                        pg.addPort(e_id.second, vertex_map[dst].getPort(e_id.second.second).getAttribute());
                        Edge<V,P,E> edge = this->getEdge(e_id);
                        pg.addEdge(edge.EdgeId(), edge.getAttribute());
                    }
                }
            }
        }

        // new induced port graph
        return pg;
    }

    ///  Induced Graph by edge sub set
    PortGraph inducedGraph(bool(*pred)(vport_id, vport_id, E)) {
        PortGraph<V,P,E> pg = PortGraph<V,P,E>(this->isReversed());
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for (auto it = adj_list.begin(); it != adj_list.end(); it++) {
            for (auto e: (*it).second) {
                // edge satisfies Pred
                vport_id src = e.EdgeId().first;
                vport_id dst = e.EdgeId().second;
                if (pred(src, dst,e.getAttribute())) {
                    //add edge
                    pg.addVport(src, vertex_map[src.first].getAttribute(),vertex_map[src.first].getPorts()[src.second].getAttribute());
                    pg.addVport(dst, vertex_map[dst.first].getAttribute(),vertex_map[dst.first].getPorts()[dst.second].getAttribute());
                    pg.addEdge(e.EdgeId(),e.getAttribute());
                }
            }
        }
        // new induced port graph
        return pg;
    }


/********** is Reachable **********/

    /* Return true if dest is reachable from source
       Else return false
    */
    bool isReachable(vertex_id source, vertex_id dest) {
        BFSVertexIterator<V,P,E> itr = BFSVertexIterator<V,P,E>(this, source);

        for (; itr != vertexEnd(); itr = itr.next()) {
            if (*itr == dest)
                return true;
        }
        return false;
    }

    /* Return true if dest is reachable from source
       Else return false
    */
    bool isReachable(vport_id source, vport_id dest) {
        BFSIterator<V,P,E> itr = BFSIterator<V,P,E>(this, source);
        for (; itr != vportEnd(); itr = itr.next()) {
            if (*itr == dest)
                return true;
        }
        return false;
    }

/********** Shortest Paths **********/
    /*
     * Returns a vector of edges that represent the shortest path from source to dest (path)
     * If there's no path then return empty vector
     */
    Path shortestPath(WeightFunction wf, vport_id src, vport_id dst, bool newWeights = true) {
        vport_pair_id pair_id = pair<vport_id, vport_id>(src, dst);
        // if using a new weight function, clear cache
        if (newWeights) {
            shortest_paths_weights.clear();
            shortest_paths.clear();
        }
        try {
            Path p = shortest_paths.at(pair_id);
        } catch (...) {
            // failed, calculate all shortest paths from source
            vector<vport_id> vp_id_vec = getVports();
            map<vport_id, vport_id> prev;
            for (auto vp_id : vp_id_vec) {
                shortest_paths_weights.insert(pair<vport_pair_id, double>(vport_pair_id(src, vp_id), DBL_MAX));
                prev[vp_id] = vport_id(-2, -2); //undefined
            }
            shortest_paths_weights.at(vport_pair_id(src, src)) = 0;
            typedef pair<vport_id, double*> id_weight_pair;
            class id_weight_pair_comparator {
            public:
                bool operator()(const id_weight_pair& lhs, const id_weight_pair& rhs)
                {
                    return *(lhs.second) > *(rhs.second);
                }
            };
            priority_queue<id_weight_pair, vector<id_weight_pair>, id_weight_pair_comparator> pq;
            for (auto vp_id : vp_id_vec) {
                pq.push(id_weight_pair(vp_id, &(shortest_paths_weights.at(vport_pair_id(src, vp_id)))));
            }
            while (!pq.empty()) {
                id_weight_pair current = pq.top();
                pq.pop();
                auto outgoingEdges = getOutgoingEdges(current.first);
                for (auto e : outgoingEdges) {
                    double dist = *current.second + wf(e);
                    if (dist < shortest_paths_weights.at(vport_pair_id(src, e.second))) {
                        shortest_paths_weights.at(vport_pair_id(src, e.second)) = dist;
                        prev[e.second] = current.first;
                    }
                }
            }

            // save paths in cache
            for (auto vp_id : vp_id_vec) {
                Path p1;
                vport_id x = vp_id;
                while(x != src) {
                    vport_id y = prev[x];
                    if (y == vport_id(-2, -2)) {
                        p1.clear();
                        break;
                    }
                    p1.push_back(edge_id(y, x));
                    x = y;
                }
                std::reverse(p1.begin(), p1.end());
                shortest_paths[vport_pair_id(src, vp_id)] = p1;
            }
            return shortest_paths[pair_id];
        }
    }

    double shortestPathWeight(WeightFunction wf, vport_id src, vport_id dst, bool newWeights = true) {
        // if using a new weight function, clear cache
        if (newWeights) {
            shortest_paths_weights.clear();
            shortest_paths.clear();
        }
        vport_pair_id pair_id = vport_pair_id(src, dst);
        double shortest_path_weight = DBL_MAX;
        try {
            shortest_path_weight = shortest_paths_weights.at(pair_id);
        } catch (const std::out_of_range& oor) {
            shortestPath(wf, src, dst, newWeights);
            // this shouldn't fail
            try {
                shortest_path_weight = shortest_paths_weights.at(pair_id);
            }
            catch (const std::out_of_range& oor) {
                printf("Error\n");
                assert(0);
            }
        }
        return shortest_path_weight;
    }

/********** Max Flow **********/

   int maxFlowAux(map<edge_id, int>& capacity_map, map<vport_id, vport_id>& previous, vport_id src, vport_id dst) {
       previous[src] = vport_id(-1, -1); // undefined
       queue<pair<vport_id, int>> queue1;
       queue1.push(pair<vport_id, int>(src, std::numeric_limits<int>::max()));
       while (!queue1.empty()) {
           vport_id curr = queue1.front().first;
           int curr_flow = queue1.front().second;
           queue1.pop();
           for (auto neighbor : adjacency_list_undirected[curr]) {
               if (previous.count(neighbor) == 0 && capacity_map[edge_id(curr, neighbor)] != 0) {
                   previous[neighbor] = curr;
                   int flow = min(curr_flow, capacity_map[edge_id(curr, neighbor)]);
                   if (neighbor == dst)
                       return flow;
                   queue1.push(pair<vport_id, int>(neighbor, flow));
               }
           }
       }
       return 0;
   }

   int maxFlow(CapacityFunction cf, vport_id src, vport_id dst) {
       int max_flow = 0;
       map<edge_id, int> capacity_map;
       vector<edge_id> edges = getEdges();
       for (auto edge : edges)
           capacity_map[edge] = cf(edge);
       map<vport_id, vport_id> previous;
       int flow = maxFlowAux(capacity_map, previous, src, dst);
       while (flow != 0) {
           max_flow += flow;
           vport_id curr = dst;
           while (curr != src) {
               vport_id parent = previous[curr];
               capacity_map[edge_id(parent, curr)] -= flow;
               capacity_map[edge_id(curr, parent)] += flow;
               curr = parent;
           }
           flow = maxFlowAux(capacity_map, previous, src, dst);
       }
       return max_flow;
   }

/******************* Clique *******************/
    /* a method for inertnal use only , used to find vports or vertices Cliques
       with bruteforce method.
       input   : cadidates :- vector of vports or vertices ids with k-1 degree from 'this' port graph
                 current_set :- vector of vports or vertices ids that construct the current graph
                 K :- the size of the Clique
       compute : for each vport\vertex current from candidates assume its in the clique or Not
                 first check if current is a neighbor with all the nodes in 'current_set'
                 if so assume its in the clique ,aslo in 'current_set' and trie to find a Clique with size k-1
                 if the Clique is found then return the 'current_set' as a nodes of the Clique Else
                 assume the current node isnt in the Clique and trie find a Clique of size k
                 form the remaining node in 'candidates' if found return 'current_set' Else
                 return False
       comment : by finishing th computation 'current_set' has a Clique of size K if present
                 else unvalid set .
     */
protected:
    template <class T>
    bool findClique(vector<T>& candidates,vector<T>& current_set,int k) {
        if (current_set.size() == k)
            return true;
        if (candidates.empty())
            return false;
        T to_add = candidates.back();
        candidates.pop_back();

        // assume not in Clique
        bool res1 = findClique(candidates, current_set, k);
        if (res1)return true;

        //assume in clique
        bool in_Clique = true;
        for (T dst : current_set) {
            in_Clique = in_Clique && this->isNeighbers(to_add, dst) && this->isNeighbers(dst, to_add);
        }
        if (!in_Clique) {
            // not a solution
            candidates.push_back(to_add);
            return false;
        }
        // in Clique
        current_set.push_back(to_add);
        bool res2 = findClique(candidates, current_set, k);
        if (res2) return true;
        // not a solution
        current_set.pop_back();
        candidates.push_back(to_add);
        return false;
    }

public:
    /* returns true if the graph has no vertices */
    bool Empty() {
        return vertex_map.empty();
    }

    /// vport version for finding Clique
    /* returns a Clique of K vports as a port graph
       first it founds all the vports with k-1 dergree named as candidates
       then with bruteforce it tries to constuct a Clique of size K
       if found returns it as a port graph Else error
    */
    PortGraph<V,P,E> findVportClique(int k) {
        assert(k>=2);
        vector<vport_id> res;
        vector<vport_id> candidates;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for (auto it : adj_list) {
            vport_id id = it.first;
            if(adj_list[id].size() >= k-1)
                candidates.push_back(id);
        }
        findClique(candidates,res,k);
        PortGraph<V,P,E> pg = PortGraph<V,P,E>();
        for (auto id : res) {
            pg.addVport(id,vertex_map[id.first].getAttribute(),vertex_map[id.first].getPort(id.second).getAttribute());
        }
        // create new graph
        for (int i = 0; i < res.size()-1;i++)
            for (int j = i+1;j < res.size();j++) {
                edge_id e_id = edge_id(res[i],res[j]);
                edge_id rev_e_id = edge_id(res[j],res[i]);
                pg.addEdge(e_id,getEdge(e_id).getAttribute());
                pg.addEdge(rev_e_id,getEdge(rev_e_id).getAttribute());
            }
        // new Clique
        return pg;
    }

    /// vertex version for finding Clique
    /* returns a Clique of K vertices as a port graph
       first it founds all the vertices with k-1 dergree named as candidates
       then with bruteforce it tries to constuct a Clique of size K
       if found returns a vector of the vertices's ids Else error
    */
    vector<int> findVertexClique(int k) {
        assert(k>=2);
        vector<int> res;
        vector<int> candidates;
        for (auto it : vertex_map){
            if (it.second.outGoingEdges() >= k-1)
                candidates.push_back(it.first);
        }
        findClique(candidates,res,k);
        // check -
        PortGraph<V,P,E> pg;
        for (auto vertex_id : res){
            for (auto it : vertex_map[vertex_id].getPorts()) {
                vport_id id = vport_id(vertex_id,it.first);
                pg.addVport(id,vertex_map[vertex_id].getAttribute(),it.second.getAttribute());
            }
        }
        return res;
    }

/****************** SubGraph ******************/
    /* returns True if This graph is sub graph of that Graph
       check if all the vertcies in 'this' graph are a sub set of 'that' vertices also check the attributs if vertex_attr_check is true
       check if all the vports in 'this' graph are a sub set of 'that' vports also check the attributs if vports_attr_check is true
       check if all the edges in 'this' graph are a sub set of 'that' edges also check the attributs if edges_attr_check is true
       if all the above are sub set of that graph return True Else False
    */
    bool isSubGraph(PortGraph<V,P,E>& that,bool vertex_attr_check,bool ports_attr_check,bool edge_attr_check) {
        const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& sub_adj_list = that.getAdjList();

        // vertex - check
        for (int vertex_id : that.getVertices()) {
            // vertex id
            if(vertex_map[vertex_id] == vertex_map.end())
                return false;
            // vertex attribute
            if(vertex_attr_check && that.getVertexAttr(vertex_id) != vertex_map[vertex_id].getAttribute())
                return false;
        }

        // vports - check
        for(auto id : that.getVports()) {
            // vport id
            if(adj_list.find(id) == adj_list.end())
                return false;
            // port attribute
            if(ports_attr_check && that.getPortAttr(id.first,id.second) != getPortAttr(id.first,id.second))
                return false;
        }

        // edge - check
        for (auto it = sub_adj_list.begin(); it != sub_adj_list.end(); it++) {
            for(auto e : it->second){
                edge_id e_id = e.EdgeId();
                // edge id
                if(adj_list[e_id.first].find(e_id.second) == adj_list[e_id.first].end())
                    return false;
                // edge attribute
                if(edge_attr_check && that.getEdge(e_id).getAttribute() != getEdge(e_id).getAttribute())
                    return false;
            }
        }

        // is sub graph
        return true;
    }

/********* Difference between two PG *********/

    /* removes the edge with the id 'edge id' if present Else error
       first of all , remove the edge from the adjancecy list of vport with the id 'edge_id.first'
       then check if the src vertex and the dst vertex still neighbors if so finish Esle
       decrease the number of out doing edges from 'src' vertex and remove 'dst' vertex from its neighbors list
       decrease the number of in going edges to 'dst' vertex and remove 'src' vertex' from its backwards neihbors list
    */
    void removeEdge(edge_id id) {
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        if(adj_list[id.first].find(getEdge(id)) == adj_list[id.first].end())
            return;
        adj_list[id.first].erase(getEdge(id));
        reverseGraph();
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& rev_adj_list = AdjacencyList();
        rev_adj_list[id.second].erase(getEdge(edge_id(id.second,id.first)));
        reverseGraph();
        vport_id src = id.first;
        vport_id dst = id.second;
        if(isNeighbers(src.first,dst.first))
            return;
        vertex_map[src.first].decOutgoingEdges();
        vertex_map[dst.first].decIngoingEdges();
        vertex_neighbors[src.first].erase(dst.first);
        reverse_vertex_neighbors[dst.first].erase(src.first);
    }

    /* removes the vport with the id 'vport id' if present Else error
       first of all removes the port with the id 'vport id.second' from the set of ports
       in 'vport id.first' vertex and reomve each directed edge from or to 'vport id' port
    */
    void removePort(vport_id id) {
        if(vertex_map.find(id.first) == vertex_map.end())
            return;
        vertex_map[id.first].removePort(id.second);
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for(auto it1 = adj_list.begin(); it1 != adj_list.end() ; it1++){
            vport_id this_vport = it1.first;
            for(auto it2 = adj_list.begin(); it2 != adj_list.end() ; it2++){
                vport_id that_vport = it2.first;
                if(this_vport == id)
                    removeEdge(edge_id(id,that_vport));
                if(that_vport == id)
                    removeEdge(edge_id(this_vport,id));
            }
        }
    }

    /* reomve the vertex with the id 'vertex id' if present Else error
       first of all remove the vertex 'id' form the "vertex_map"
       and use removePort func to remove each port 'in' this vertex and thus including
       the removal of all out\in going edges from\to 'this' vertex
    */
    void removeVertex(int id) {
        if(vertex_map.find(id) == vertex_map.end())
            return;
        for(auto p :vertex_map[id].getPorts()){
            removePort(vport_id(id,p.getPortId()));
        }
        vertex_map.erase(id);
    }

    /* return the number of ports in 'id' vertex if present Else error */
    int getVertexPortNum(int id) {
        return vertex_map[id].getPorts().size();
    }

    /* return Port Graph without any vertices/vports/edges that appears in org_pg */
    PortGraph<V,P,E> diff(PortGraph<V,P,E>& org_pg) {
        PortGraph<V,P,E> diff_pg = *this;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = diff_pg.getAdjList();
        for(auto edge_set : org_pg.getAdjList())
            for(auto edge : edge_set)
                diff_pg.removeEdge(edge.EdgeId());

        for(auto it = adj_list.begin(); it != adj_list.end();it++) {
            vport_id id = it.first;
            if(adj_list[id].size() == 0)
                diff_pg.removeVport();
        }

        for(int id :  diff_pg.getVertices()){
            if(diff_pg.getVertexPortNum(id) == 0)
                diff_pg.removeVertex(id);
        }
        // diff pg
        return diff_pg;
    }

};

/********** Iterator Implementations **********/

/// DFS BY VPORT
template <class V , class P , class E >
class DFSIterator : public PGVportIterator {
private:
    vector<vport_id> path;
    PortGraph<V,P,E>* pg;
    map<vport_id,bool> visited;
    map<vport_id,bool> not_visited;

protected:
    DFSIterator(vector<vport_id> _path,PortGraph<V,P,E>* _pg,map<vport_id,bool>& _visited,map<vport_id,bool>& _not_visited,vport_id _current) {
        path = _path;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    DFSIterator(vport_id id) {
        current = id;
    }

    DFSIterator(PortGraph<V,P,E>* _pg,vport_id src) {
        this->pg = _pg;
        // if src not found in PG
        path.push_back(src);
        visited.insert(pair<vport_id,bool>(src,true));
        current = src;
        auto x = pg->getVports();
        for(auto p : x) {
            if(p != src)
                not_visited.insert(pair<vport_id, bool>(p, true));
        }
    }

    ~DFSIterator(){}

    DFSIterator operator++() {
        while(!not_visited.empty()) {
            vport_id current_src = current;
            for(Edge<V,P,E> e : pg->getVportAdjList(current_src)) {
                vport_id dst = e.EdgeId().second;
                if(visited.find(dst) == visited.end()){
                    //assert(not_visited.find(dst) != not_visited.end());
                    visited.insert(pair<vport_id,bool>(dst,true));
                    not_visited.erase(dst);
                    path.push_back(dst);
                    current = dst;
                    return (*this);
                }
            }
            path.pop_back();
            if (!path.empty())
                current = path.back();
            else {// get random new 'not visited' vport
                auto it = not_visited.begin();
                current = (*it).first;
                path.push_back(current);
                visited.insert(pair<vport_id,bool>(current,true));
                not_visited.erase(current);
                return (*this);
            }
        }
        current = this->END_VPORT;
        return *this;
    }

    DFSIterator<V,P,E> operator++(int) {
        DFSIterator<V,P,E> old = DFSIterator<V,P,E>(path,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    vport_id operator*() const {
        return current;
    }

    bool operator== (const PGVportIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGVportIterator& that)const {return !this->operator == (that) ;}
};

/// DFS BY VERTEX
template <class V , class P , class E >
class DFSVertexIterator : public PGVertexIterator {
private:
    vector<int> path;
    PortGraph<V,P,E>* pg;
    map<int,bool> visited;
    map<int,bool> not_visited;

protected:
    DFSVertexIterator(vector<int> _path,PortGraph<V,P,E>* _pg,map<int,bool>& _visited,map<int,bool>& _not_visited,int _current) {
        path = _path;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    DFSVertexIterator(int id) {
        current = id;
    }

    DFSVertexIterator(PortGraph<V,P,E>* _pg, int src) {
        this->pg = _pg;
        // if src not found in PG
        path.push_back(src);
        visited.insert(pair<int,bool>(src,true));
        current = src;
        for(auto p : pg->getVertices()) {
            if(p != src)
                not_visited.insert(pair<int, bool>(p, true));
        }
    }

    ~DFSVertexIterator() {}

    DFSVertexIterator operator++() {
        while (!not_visited.empty()) {
            int current_src = current;
            //done
            for(int dst : pg->getVertexAdjList(current_src)) {
                //vport_id dst = e.EdgeId().second;
                if(visited.find(dst) == visited.end()){
                    //assert(not_visited.find(dst) != not_visited.end());
                    visited.insert(pair<int,bool>(dst,true));
                    not_visited.erase(dst);
                    path.push_back(dst);
                    current = dst;
                    return (*this);
                }
            }
            //done
            path.pop_back();
            if (!path.empty())
                current = path.back();
            else {// get random new 'not visited' vertex
                auto it = not_visited.begin();
                current = (*it).first;
                path.push_back(current);
                visited.insert(pair<int,bool>(current,true));
                not_visited.erase(current);
                return (*this);
            }
        }
        current = this->END_VERTEX;
        return *this;
    }

    DFSVertexIterator<V,P,E> operator++(int) {
        DFSVertexIterator<V,P,E> old = DFSVertexIterator<V,P,E>(path,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    int operator*()const {
        return current;
    }

    bool operator== (const PGVertexIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGVertexIterator& that)const {return !this->operator == (that) ;}
};

/// BFS BY VPORT
template <class V , class P , class E >
class BFSIterator : public PGVportIterator {
private:
    vector<vport_id> queue;
    PortGraph<V,P,E>* pg;
    map<vport_id,bool> visited;
    map<vport_id,bool> not_visited;

protected:
    BFSIterator(vector<vport_id> _queue,PortGraph<V,P,E>* _pg,map<vport_id,bool>& _visited,map<vport_id,bool>& _not_visited,vport_id _current) {
        queue = _queue;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    BFSIterator(vport_id id) {
        current = id;
    }

    ~BFSIterator(){}

    BFSIterator(PortGraph<V,P,E>* _pg,vport_id src) {
        this->pg = _pg;
        // if src not found in PG
        queue.push_back(src);
        visited.insert(pair<vport_id,bool>(src,true));
        current = src;
        for(auto p : pg->getVports()) {
            if(p != src)
                not_visited.insert(pair<vport_id, bool>(p, true));
        }
    }

    BFSIterator operator++() {
        if(queue.empty()){
            current = this->END_VPORT;
            return *this;
        }
        vport_id current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for(Edge<V,P,E> e : pg->getVportAdjList(current_src)) {
            vport_id dst = e.EdgeId().second;
            if (visited.find(dst) == visited.end()) {
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<vport_id,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if (!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else {// get new random 'not visited' vport
            if (not_visited.empty()) {
                current = this->END_VPORT;
                return *this;
            }
            auto it = not_visited.begin();
            current = (*it).first;
            queue.push_back(current);
            visited.insert(pair<vport_id,bool>(current,true));
            not_visited.erase(current);
        }
        return (*this);
    }

    BFSIterator next() {
        if(queue.empty()){
            current = this->END_VPORT;
            return *this;
        }
        vport_id current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for (Edge<V,P,E> e : pg->getVportAdjList(current_src)){
            vport_id dst = e.EdgeId().second;
            if (visited.find(dst) == visited.end()){
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<vport_id,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if (!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else { // return end
            current = this->END_VPORT;
            return *this;
        }
        return (*this);
    }

    BFSIterator<V,P,E> operator++(int) {
        BFSIterator<V,P,E> old = BFSIterator<V,P,E>(queue,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    vport_id operator*()const {
        return current;
    }

    bool operator== (const PGVportIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGVportIterator& that)const {return !this->operator==(that);}
};

/// BFS BY VERTEX
template <class V , class P , class E >
class BFSVertexIterator : public PGVertexIterator {
private:
    vector<int> queue;
    PortGraph<V,P,E>* pg;
    map<int,bool> visited;
    map<int,bool> not_visited;

protected:
    BFSVertexIterator(vector<int> _queue,PortGraph<V,P,E>* _pg,map<int,bool>& _visited,map<int,bool>& _not_visited,int _current) {
        queue = _queue;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }

public:

    BFSVertexIterator(int id) {
        current = id;
    }

    ~BFSVertexIterator(){}

    BFSVertexIterator(PortGraph<V,P,E>* _pg,int src) {
        this->pg = _pg;
        // if src not found in PG
        queue.push_back(src);
        visited.insert(pair<int,bool>(src,true));
        current = src;
        for(auto p : pg->getVertices()) {
            if(p != src)
                not_visited.insert(pair<int, bool>(p, true));
        }
    }

    BFSVertexIterator operator++() {
        if(queue.empty()){
            current = this->END_VERTEX;
            return *this;
        }
        int current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for (int dst : pg->getVertexAdjList(current_src)) {
            if (visited.find(dst) == visited.end()) {
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<int,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if (!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else {// get new random 'not visited' vertex
            if(not_visited.empty()) {
                current = this->END_VERTEX;
                return *this;
            }
            auto it = not_visited.begin();
            current = (*it).first;
            queue.push_back(current);
            visited.insert(pair<int,bool>(current,true));
            not_visited.erase(current);
        }
        return (*this);
    }

    BFSVertexIterator next() {
        if(queue.empty()){
            current = this->END_VERTEX;
            return *this;
        }
        int current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for (int  dst : pg->getVertexAdjList(current_src)) {
            if (visited.find(dst) == visited.end()) {
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<int,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if (!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else { // return end
            current = this->END_VERTEX;
            return *this;
        }
        return (*this);
    }

    BFSVertexIterator<V,P,E> operator++(int) {
        BFSVertexIterator<V,P,E> old = BFSVertexIterator<V,P,E>(queue,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    int operator*() const {
        return current;
    }

    bool operator== (const PGVertexIterator& that) const {return (this->current) == (*that) ;}
    bool operator!= (const PGVertexIterator& that) const {return !this->operator == (that) ;}
};

#endif
