//
// Created by Mario Barbara and Feras Bisharat
//

#ifndef PROJECT1_PORTGRAPH_H
#define PROJECT1_PORTGRAPH_H
#include "Utilities.h"
using namespace std;

template <class P>
#define PortMap map<int,Port<P>>
//struct cmpPort{
//    bool operator () (const Port<P>& a,const Port<P>& b){
//        return a.getPortId() < b.getPortId() ;
//    }
//};

class Port{
public:

    explicit Port(int portId = -1, P attr = P()) : port_id(portId) , attribute(attr) {}

    int getPortId() const { return port_id; }

    const P getAttribute()  { return attribute; }

    bool operator==(const Port& p) { return p.getPortId() == port_id; }

    bool operator<(const Port& p){
        return this->port_id < p.getPortId();
    }

    void setAttribute(P& attr){
        attribute = attr;
    }

    void print(){
        ostringstream s;
        s << "Port " << port_id << " with attribute " << attribute << endl;
        fprintf(stderr, s.str().c_str());
    }

private:

    int port_id;
    P attribute;
};

//template <class V, class P>
//    struct cmpVertex{
//        bool operator () (const int& a,const int& b) const{
//            return a < b ;
//        }
//    };

template <class V, class P>
class Vertex {
public:
    Vertex () {}
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

    // new version
    void incOutgoingEdges(){
        n_outgoing_edges++;
    }

    void incIngoingEdges(){
        n_ingoing_edges++;
    }

    void decOutgoingEdges(){
        n_outgoing_edges--;
    }

    void decIngoingEdges(){
        n_ingoing_edges--;
    }

    int outGoingEdges(){ return n_outgoing_edges;}

    int inGoingEdges(){ return n_ingoing_edges;}

    int vertexId() { return vertex_id; }

    V getAttribute() { return attribute; }

    PortMap getPorts() { return ports; }

    int getPortsNum() { return ports.size(); }

    bool operator==(Vertex v) { return v.vertexId() == vertex_id; }

    bool operator<(Vertex v){
        return  v.vertexId() > vertex_id;
    }

    Port<P> getPort(const int port_id){
        //check -
        return ports[port_id];
    }

    void setVertexId(int id){
        vertex_id = id;
    }

    void setAttribute(V attr){
        attribute = attr;
    }

    void setPortAttr(int port_id,P attr){
        Port<P> p = Port<P>(port_id, attr);
        ports[port_id] = p;
    }

    void print() {
        ostringstream s;
        s << "Vertrex " << vertex_id << " with attribute " << attribute << ", with " << ports.size() << " ports." << endl;
        s << "Vertrex " << vertex_id << " with Ports  : " << endl ;
        fprintf(stderr, s.str().c_str());
        for (auto p : ports)
            p.second.print();
    }

//    Vertex<V,P>& operator=(Vertex<V,P> const &a){
//        this->vertex_id = a.vertexId();
//        this->attribute = a.getAttribute();
//        this->ports = a.getPorts();
//   }
    PortMap ports;
private:
    int vertex_id;
    V attribute;
    int n_outgoing_edges;
    int n_ingoing_edges;

};

template <class V, class P, class E>
class Edge {
public:
    Edge() = default ;
    Edge(vport src , vport dst , E attr = E()): source(src), dest(dst), attribute(attr) {
        vport_id p1 = std::make_pair(src.first.vertexId(),src.second.getPortId());
        vport_id p2 = std::make_pair(dst.first.vertexId(),dst.second.getPortId());
        id = std::make_pair(p1, p2) ;
    }

    vport getSource() { return source; }

    vport getDest() { return dest; }

    edge_id EdgeId() const { return id; }

    bool operator==(Edge e) { return e.id == id; }

    E getAttribute() { return attribute; }

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

    void print() {
        ostringstream s = toString();
        fprintf(stderr, s.str().c_str());
    }

    void printIds() {
        ostringstream s = toStringIds();
        fprintf(stderr, s.str().c_str());
    }

private:
    vport source;
    vport dest;
    edge_id id; // if ther is a problrm with the COMP make the feiled public
    E attribute;
};

struct cmpVport {
    bool operator()(const vport_id& a,const vport_id& b) const  {
        return (a.first < b.first) || (a.first == b.first && a.second < b.second);
    }
};

template <class V = int, class P = int, class E = int>
struct cmpEdge {
    // check-
    bool operator()(const Edge<V, P, E>& a,const Edge<V, P, E>& b) const {
        return a.EdgeId() < b.EdgeId() ;
    }
};

template <class V , class P , class E >
class DFSIterator;

template <class V , class P , class E >
class BFSIterator;

template <class V , class P , class E >
class DFSVertexIterator;

template <class V , class P , class E >
class BFSVertexIterator;

template <class V = int, class P = int, class E = int>
class PortGraph {
private:

    // maps vport vp (vertix and port) to all the edges that has an vp as a source
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adjacency_list;
    // maps vport vp (vertix and port) to all the transpose edges that has an vp as a source
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> transpose_adjacency_list;
    // maps vport vp (vertix and port) to all the vports that has an vp as a source or dest
    map<vport_id, set<vport_id>> adjacency_list_undirected;


    // flag for the transpose Graph
    bool is_transpose ;

    // maps vport_id (vertix and port) to its vport
    // vertex id -> vertex
    // vmap por_id -> attr
    map<int, Vertex<V,P>> vport_map;

    //maps vertex_id -> neighbors :- vertices
    map<int,map<int,bool>> vertex_neighbors;
    //maps vertex_id -> neighbors :- vertices for the transpose graph
    map<int,map<int,bool>> transpose_vertex_neighbors;

    // params for SCC
    const int UNSEEN = -1;
    const int SEEN = 1;

    // for Algo's in case of transpose
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& AdjacencyList(){
        return is_transpose ? transpose_adjacency_list : adjacency_list ;
    }

    // for Algo's in case of transpose
    map<int,map<int,bool>>& VertexNeighbors(){
        return is_transpose ? transpose_vertex_neighbors : vertex_neighbors ;
    }

    // for shortest paths
    typedef vector<edge_id> Path;
    map<pair<vport_id, vport_id>, double> shortest_paths_weights;
    map<pair<vport_id, vport_id>, Path> shortest_paths;

public:
    // Build PortGraph with
    PortGraph() = default;

    PortGraph(int n_vertices, vector<int> ports_num, vector<edge_id> edges_list,
              VerticesAttributes verticesAttributes =  VerticesAttributes(),
              vector<PortsAttributes> portsAttributes = vector<PortsAttributes>(),
              EdgesAttributes edgesAttributes = EdgesAttributes())
    {
        is_transpose = false;
        for (int i = 0; i < n_vertices; ++i) {
            if(verticesAttributes.empty() && portsAttributes .empty())
                addVertix(i, ports_num[i]);
            else if(portsAttributes .empty())
                addVertix(i, ports_num[i],verticesAttributes[i]);
            else if(verticesAttributes.empty())
                addVertix(i, ports_num[i],V(),portsAttributes[i]);
            else {
                addVertix(i, ports_num[i],verticesAttributes[i],portsAttributes[i]);
            }
        }

        //add transpose edge
        for (int i = 0; i < edges_list.size(); ++i) {
            if(edgesAttributes.empty())
                addEdge(edges_list[i]);
            else {
                addEdge(edges_list[i], edgesAttributes[i]);
            }
        }

    }

    void transposeGraph(){
        //check
        is_transpose = !is_transpose;
    }

    bool isTranspose(){
        return is_transpose;
    }

    void setTranspose(bool _is_transpose){
        is_transpose = _is_transpose;
    }

    void addVertix(int vertex_id, int ports_num, V attr = V(), PortsAttributes ports_attr = PortsAttributes())
    {
        for (int i = 0; i < ports_num; ++i) {
            adjacency_list[vport_id(vertex_id, i)] = set<Edge<V, P, E>, cmpEdge<V,P,E>>();
            // transpose Graph
            transpose_adjacency_list[vport_id(vertex_id, i)] = set<Edge<V, P, E>, cmpEdge<V,P,E>>();
            adjacency_list_undirected[vport_id(vertex_id, i)] = set<vport_id>();
            P portAttr =  ports_attr.empty() ? P() : ports_attr[i] ;
        }
        vport_map[vertex_id] = Vertex<V, P>(vertex_id, ports_num, attr, ports_attr);
        if(vertex_neighbors.find(vertex_id) == vertex_neighbors.end())
            vertex_neighbors[vertex_id] = map<int,bool>();
        transpose_vertex_neighbors[vertex_id] = map<int,bool>();
    }

    void addEdge(edge_id id, E attr = E())
    {
        // add edge
        Vertex<V,P> v1 = vport_map[id.first.first];
        Vertex<V,P> v2 = vport_map[id.second.first];
        vport vp1 = vport(v1,v1.getPort(id.first.second));
        vport vp2 = vport(v2,v2.getPort(id.second.second));
        adjacency_list[id.first].insert(Edge<V, P, E>(vp1, vp2, attr) );
        // add transpose Edge
        transpose_adjacency_list[id.second].insert(Edge<V, P, E>(vp2, vp1, attr));
        adjacency_list_undirected[id.first].insert(id.second);
        adjacency_list_undirected[id.second].insert(id.first);
        //check
        if(vertex_neighbors[id.first.first].find(id.second.first) == vertex_neighbors[id.first.first].end())
            vertex_neighbors[id.first.first].insert(pair<int,bool>(id.second.first,true));
        if(transpose_vertex_neighbors[id.second.first].find(id.first.first) == transpose_vertex_neighbors[id.second.first].end())
            transpose_vertex_neighbors[id.second.first].insert(pair<int,bool>(id.first.first, true));
        //assume not exist
        vport_map[id.first.first].incOutgoingEdges();
        vport_map[id.second.first].incIngoingEdges();
    }

    //Print the vertcies and ports per vertex
    void Print()
    {
        for (auto x : vport_map) {
            x.second.print();
        }
    }

    //Print the edges
    void PrintEdges() {

        for (auto it = adjacency_list.begin(); it != adjacency_list.end(); ++it)
            for(auto e : (*it).second)
                e.print();
    }

    // get outgoing edges from vport id
    vector<edge_id> getOutgoingEdges(vport_id id) {
        auto adj_list = AdjacencyList();
        vector<edge_id> to_return;
        for(auto e : adj_list[id]) {
            to_return.push_back(e.EdgeId());
        }
        return to_return;
    }

    vector<edge_id> getEdges() {
        vector<edge_id> to_return;
        auto adj_list = AdjacencyList();
        for (auto pr : adj_list) {
            for (auto edge : pr.second)
                to_return.push_back(edge.EdgeId());
        }
        return to_return;
    }

/********** Topological Sort **********/

// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: is g a DAG (return value), a topological ordering of g (order).
// comment: order is valid only if g is a DAG.
    bool topological_sort(){
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
        for(auto id : vport_order){
            for (auto edge : adj_list[id]) {
                vport_id neighbor = edge.EdgeId().second;
                vport_indegree[neighbor] -- ;
                if (vport_indegree[neighbor] == 0)
                    vport_order.push_back(neighbor);
            }
        }
        // check -
        return vport_order.size() == adj_list.size();
    }

/********** Strongly Connected Components **********/

    void KosarajuDFS (const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list
            ,vport_id id, map<vport_id,vport_id>& S, map<vport_id,int>& colorMap, int color){

        colorMap[id] = color;
        for(auto edge : adj_list[id]){
            vport_id neighbor = edge.EdgeId().second;
            if(colorMap[neighbor] == UNSEEN)
                KosarajuDFS(neighbor, S,colorMap,color);
        }
        S.insert({id,id});
    }

// Compute the number of SCCs and maps nodes to their corresponding SCCs.
// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: the number of SCCs (return value), a mapping from node to SCC color (components).
    int findSCC(map<vport_id ,int>& components){
        // first pass: record the `post-order' of original graph
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adj_list = AdjacencyList();
        map<vport_id ,vport_id> postOrder , dummy ;
        map<vport_id ,int>  seen ;
        for(auto it = adj_list.begin(); it != adj_list.end(); ++it){
            vport_id id = (*it).first;
            seen[id] = UNSEEN;
            components[id] = UNSEEN;
        }
        for(auto it = adj_list.begin(); it != adj_list.end(); ++it){
            vport_id id = (*it).first;
            if(seen[id] == UNSEEN)
                KosarajuDFS(adj_list , id,postOrder,seen,SEEN);
        }
        // reverse iterator
        // second pass: explore the SCCs based on first pass result
        int numSCC = 0;
        for(auto rit = adj_list.rbegin(); rit != adj_list.rend(); ++rit){
            vport_id id = (*rit).first;
            transposeGraph();
            map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> rev_adj_list = AdjacencyList();
            transposeGraph();
            if(components[postOrder[id]] == UNSEEN)
                KosarajuDFS(rev_adj_list, postOrder[id], dummy, components, numSCC++);
        }
        return numSCC;
    }

// Computes the SCC graph of a given digraph.
// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: strongly connected components graph of g (sccg).
    void findSCCgraph(vector<set<int>>& sccg){
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adj_list = AdjacencyList();
        map<vport_id, int> component ;
        int n = findSCC(component);
    }

/********** Min Spanning Tree **********/
private:
    struct unionfindPG  {
        map<vport_id ,int> rank;
        map<vport_id ,vport_id >parent;

        explicit unionfindPG (map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list ) {
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

// input: edges v1->v2 of the form (weight,(v1,v2)),
//        number of nodes (n), all nodes are between 0 and n-1.
// output: weight of a minimum spanning tree.
public:
    double Kruskal(double(*weightFunc)(const edge_id ,const E attr)){
        int n = AdjacencyList().size();
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
                mst_cost += e.first;
                components.unite(e.second.first, e.second.second);
            }
        }
        return mst_cost;
    }

/********** BFS/DFS **********/

    PGVportIterator vportEnd(){
        return PGVportIterator(vport_id(-1,-1));
    }

    PGVertexIterator vertexEnd(){
        return PGVertexIterator(-1);
    }

    vector<vport_id> getVports(){
        vector<vport_id> res;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for(auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id p = (*it).first;
            res.push_back(p);
        }
        return res;
    }

    vector<int> getVertices(){
        vector<int> res;
        for(auto it : vport_map)
            res.push_back(it.first);
        return res;
    }

    set<Edge<V,P,E>,cmpEdge<V,P,E>> getAdjList(vport_id id){
        return AdjacencyList()[id];
    }

    vector<int> getVertexAdjList(int src){
        vector<int> res;
        for(auto it : VertexNeighbors()[src])
            res.push_back(it.first);
        return res;
    }

/********** Bipartite **********/

    bool isNeighbers(const vport_id& src, const vport_id& dst){
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for(const Edge<V,P,E>& e : adj_list[src]){
            if(e.EdgeId().second == dst)
                return true;
        }
        return false;
    }

    bool isNeighbers(const int src, const int dst) {
        const PortMap mp_src = vport_map[src].getPorts();
        for(auto it1 = mp_src.begin(); it1 != mp_src.end();it1++){
            //Port<P> p1 = (*it1).second;
            vport_id id1 = vport_id(src,(*it1).second.getPortId());
            const PortMap mp_dst = vport_map[dst].getPorts();
            for(auto it2 = mp_dst.begin(); it2 != mp_dst.end();it2++){
                //Port<P> p2 = (*it2).second;
                vport_id id2 = vport_id(dst,(*it2).second.getPortId());
                if(isNeighbers(id1, id2))
                    return true;
            }
        }
        return false;
    }

// Function returns true if graph is Bipartite, else false
    bool isBipartite(vport_id src){
        // if src not in graph ...
        map<vport_id,int> colorArr;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            colorArr.insert(pair<vport_id, int>((*it).first,-1));
        }
        // Assign first color to source
        colorArr[src] = 1;
        // Create a queue (FIFO) of vertex
        // numbers and enqueue source vertex
        // for BFS traversal
        vector<vport_id> queue;
        queue.push_back(src);
        // Run while there are vertices
        // in queue (Similar to BFS)
        while(! queue.empty()){
            // Dequeue a vertex from queue
            vport_id u = queue.front();
            queue.erase(queue.begin());
            // Return false if there is a self-loop
            if(isNeighbers(u,u))
                return false;
            // Find all non-colored adjacent vertices
            for (auto it = adj_list.begin(); it !=adj_list.end(); ++it){
                vport_id v = (*it).first;
                // An edge from u to v exists and
                // destination v is not colored
                if(isNeighbers(u,v) && colorArr[v] == -1 ){
                    // Assign alternate color to this adjacent v of u
                    colorArr[v] = 1- colorArr[u];
                    queue.push_back(v);
                }
                    // An edge from u to v exists and destination
                    // v is colored with same color as u
                else if(isNeighbers(u,v) && colorArr[v] == colorArr[u])
                    return false;
            }
        }
        // If we reach here, then all adjacent
        // vertices can be colored with alternate color
        return true;
    }

/********** Induced Graph **********/

    Edge<V,E,P> getEdge(const edge_id id){
        assert(isNeighbers(id.first,id.second));
        auto adj_list = this->AdjacencyList();
        for(const Edge<V,P,E>& e : adj_list[id.first]){
            if(e.EdgeId() == id)
                return e;
        }
        //won't get here
        return Edge<V,P,E>();
    }

    V& getVertexAttr(const int i){
        return vport_map[i].getAttribute();
    }

    P& getPortAttr(const int vertex_id , const int port_id){
        return vport_map[vertex_id].getPort(port_id).getAttribute();
    }

    void addVport(vport_id id,V vertex_attr,P port_attr){
        //check -
        Vertex<V,P>  v = Vertex<V,P>();
        v.setVertexId(id.first);
        v.setAttribute(vertex_attr);
        v.setPortAttr(id.second,port_attr);
        vport_map[id.first] = v;
    }

    void addPort(vport_id id,P attr = P()){
        //check -
        vport_map[id.first].setPortAttr(id.second,attr);
    }

    //returns all the edges form vertex 'src' to vertex 'dst' ,else empty sub set
    vector<edge_id> getNeighbers(const int src, const int dst){
        vector<edge_id> res ;
        const PortMap mp_src = vport_map[src].getPorts();
        for(auto it1 = mp_src.begin(); it1 != mp_src.end();it1++){
            //Port<P> p1 = (*it1).second;
            vport_id id1 = vport_id(src,(*it1).second.getPortId());
            const PortMap mp_dst = vport_map[dst].getPorts();
            for(auto it2 = mp_dst.begin(); it2 != mp_dst.end();it2++){
                //Port<P> p2 = (*it2).second;
                vport_id id2 = vport_id(dst,(*it2).second.getPortId());
                if(isNeighbers(id1, id2))
                    res.push_back(edge_id(id1, id2));
            }
        }
        return res;
    }

    // adds vertex to pg
    // assume "new" vertex
    void addVertex(Vertex<V,P> v){
        vport_map[v.vertexId()] = v;
    }

    ///  Induced Graph by vport sub set
    PortGraph inducedGraph(bool(*pred)(const vport_id)){
        map<vport_id,int> sub_vport_set;
        // new induced graph
        PortGraph<V,P,E> pg = PortGraph<V,P,E>();
        map<int,int> hash;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        // get all vports that satisfies Pred
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if(pred(id)){
                sub_vport_set.insert(pair<vport_id,int>(id,1));
            }
        }
        // get all the edges for the induced graph
        for(auto it1 = sub_vport_set.begin();it1 != sub_vport_set.end();it1++){
            auto tmp = it1;
            tmp++;
            //check
            for(auto it2 = tmp; it2 != sub_vport_set.end();it2++){
                vport_id src = (*it1).first;
                vport_id dst = (*it2).first;
                // add new vport
                if(isNeighbers(src,dst)){
                    // new
                    if(dst.first == 6){
                        int x = 0 ;
                    }
                    int y = 0;
                    if(hash.find(src.first) == hash.end()){
                        hash[src.first] = 1;
                        pg.addVport(src,vport_map[src.first].getAttribute(),vport_map[src.first].getPort(src.second).getAttribute());
                        if(hash.find(dst.first) == hash.end()){
                            //new
                            hash[dst.first] = 1;
                            pg.addVport(dst,vport_map[dst.first].getAttribute(),vport_map[dst.first].getPort(dst.second).getAttribute());
                        }else{
                            pg.addPort(dst,vport_map[dst.first].getPort(dst.second).getAttribute());
                        }
                    }else {
                        pg.addPort(src,vport_map[src.first].getPort(src.second).getAttribute());
                        if (hash.find(dst.first) == hash.end()) {
                            //new
                            hash[dst.first] = 1;
                            pg.addVport(dst,vport_map[dst.first].getAttribute(),vport_map[dst.first].getPort(dst.second).getAttribute());
                        }else{
                            pg.addPort(dst,vport_map[dst.first].getPort(dst.second).getAttribute());
                        }
                    }
                    edge_id e_id = edge_id(src,dst);
                    P e_attr = getEdge(e_id).getAttribute();
                    //check -
                    // add rev
                    pg.addEdge(e_id,e_attr);
                }
            }
        }

        // new induced port graph
        pg.setTranspose(this->isTranspose());
        return pg;
    }

    ///  Induced Graph by vertex sub set
    PortGraph inducedGraph(bool(*pred)(const int)){
        map<int,int> sub_vertex_set;
        // new induced graph
        PortGraph<V,P,E> pg = PortGraph<V,P,E>();
        map<int,int> hash;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        // get all vports that satisfies Pred
        for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if(pred(id.first)){
                sub_vertex_set.insert(pair<int,int>(id.first,1));
            }
        }
        // get all the edges for the induced graph
        for(auto it1 = sub_vertex_set.begin();it1 !=sub_vertex_set.end();it1++){
            auto tmp = it1;
            tmp++;
            int src = (*it1).first;
            // add first vertex
            if(hash.find(src) == hash.end()){
                Vertex<V,P> v1 = Vertex<V,P>();
                v1.setVertexId(src);
                v1.setAttribute(vport_map[src].getAttribute());
                pg.addVertex(v1);
                hash[src] = 1;
            }
            for(auto it2 = tmp; it2 != sub_vertex_set.end();it2++){
                int dst = (*it2).first;
                vector<edge_id> vec = getNeighbers(src, dst);

                // add second vertex
                if(hash.find(dst) == hash.end()){
                    Vertex<V, P> v2 = Vertex<V, P>();
                    v2.setVertexId(dst);
                    v2.setAttribute(vport_map[dst].getAttribute());
                    pg.addVertex(v2);
                    hash[dst] = 1;
                }
                if(! vec.empty()){
                    //new edge
                    for(auto e_id : vec){
                        pg.addPort(e_id.first, vport_map[src].getPort(e_id.first.second).getAttribute());
                        pg.addPort(e_id.second, vport_map[dst].getPort(e_id.second.second).getAttribute());
                        Edge<V,P,E> edge = this->getEdge(e_id);
                        pg.addEdge(edge.EdgeId(), edge.getAttribute());
                    }
                }
            }
        }

        // new induced port graph
        pg.setTranspose(this->isTranspose());
        return pg;
    }

    ///  Induced Graph by edge sub set
    PortGraph inducedGraph(bool(*pred)(const vport_id, const vport_id)){
        PortGraph<V,P,E> pg = PortGraph<V,P,E>();
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        map<int,int> hash;
        for(auto it = adj_list.begin(); it != adj_list.end(); it++){
            for(auto e : (*it).second){
                // edge satisfies Pred
                vport_id src = e.EdgeId().first;
                vport_id dst = e.EdgeId().second;
                if(pred(src, dst)){
                    //add edge
                    if(hash.find(src.first) == hash.end()){
                        pg.addVport(src, vport_map[src.first].getAttribute(),vport_map[src.first].getPorts()[src.second].getAttribute());
                        if(hash.find(dst.first) == hash.end()){
                            pg.addVport(dst, vport_map[dst.first].getAttribute(),vport_map[dst.first].getPorts()[dst.second].getAttribute());
                        }else{
                            pg.addPort(src, vport_map[dst.first].getPorts()[dst.second].getAttribute());
                        }
                    }else{
                        pg.addPort(src,vport_map[src.first].getPorts()[src.second].getAttribute());
                        if(hash.find(dst.first) == hash.end()){
                            pg.addVport(src, vport_map[dst.first].getAttribute(),vport_map[dst.first].getPorts()[dst.second].getAttribute());
                        }else{
                            pg.addPort(dst, vport_map[dst.first].getPorts()[dst.second].getAttribute());
                        }
                    }
                    pg.addEdge(e.EdgeId(),e.getAttribute());
                }
            }
        }
        // new induced port graph
        pg.setTranspose(this->isTranspose());
        return pg;
    }


    /* Return true if dest is reachable from source
       Else return false
    */
    bool is_reachable(vertex_id source, vertex_id dest) {
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
    bool is_reachable(vport_id source, vport_id dest) {
        BFSIterator<V,P,E> itr = BFSIterator<V,P,E>(this, source);
        for (; itr != vportEnd(); itr = itr.next()) {
            if (*itr == dest)
                return true;
        }
        return false;
    }

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

    bool Empty(){
        return vport_map.empty();
    }

    /// vport version of finding Clique
    vector<vport_id> findVportClique(int k){
        assert(k>=2);
        vector<vport_id> res;
        vector<vport_id> candidates;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for(auto it : adj_list){
            vport_id id = it.first;
            if(adj_list[id].size() >= k-1)
                candidates.push_back(id);
        }
        findClique(candidates,res,k);
        return res;
    }

    /// vertex version of finding Clique
    vector<int> findVertexClique(int k){
        assert(k>=2);
        vector<int> res;
        vector<int> candidates;
        for(auto it : vport_map){
            if(it.second.outGoingEdges() >= k-1)
                candidates.push_back(it.first);
        }
        findClique(candidates,res,k);
        return res;
    }

    template <class T>
    bool findClique(vector<T>& candidates,vector<T>& current_set,int k){
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

/****************** SubGraph ******************/

    bool isSubGraph(PortGraph<V,P,E>& sub_graph,bool vertex_attr_check,bool ports_attr_check,bool edge_attr_check){
        const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& sub_adj_list = sub_graph.getAdjList();

        /// vertex - check
        for(int vertex_id : sub_graph.getVertices()){
            // vertex id
            if(vport_map[vertex_id] == vport_map.end())
                return false;
            // vertex attribute
            if(vertex_attr_check && sub_graph.getVertexAttr(vertex_id) != vport_map[vertex_id].getAttribute())
                return false;
        }

        /// vports - check
        for(auto id : sub_graph.getVports()){
            // vport id
            if(adj_list.find(id) == adj_list.end())
                return false;
            // port attribute
            if(ports_attr_check && sub_graph.getPortAttr(id.first,id.second) != getPortAttr(id.first,id.second))
                return false;
        }

        /// edge - check
        for(auto it = sub_adj_list.begin(); it != sub_adj_list.end(); it++){
            for(auto e : it->second){
                edge_id e_id = e.EdgeId();
                // edge id
                if(adj_list[e_id.first].find(e_id.second) == adj_list[e_id.first].end())
                    return false;
                // edge attribute
                if(edge_attr_check && sub_graph.getEdge(e_id).getAttribute() != getEdge(e_id).getAttribute())
                    return false;
            }
        }

        // is sub graph
        return true;
    }

};

/********** Iterator Implementation **********/

/// DFS BY VPORT
template <class V , class P , class E >
class DFSIterator : public PGVportIterator {
private:
    vector<vport_id> path;
    PortGraph<V,P,E>* pg;
    map<vport_id,bool> visited;
    map<vport_id,bool> not_visited;
    DFSIterator(vector<vport_id> _path,PortGraph<V,P,E>* _pg,map<vport_id,bool>& _visited,map<vport_id,bool>& _not_visited,vport_id _current){
        path = _path;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    DFSIterator(vport_id id){
        current = id;
    }

    DFSIterator(PortGraph<V,P,E>* _pg,vport_id src){
        this->pg = _pg;
        // if src not found in PG
        path.push_back(src);
        visited.insert(pair<vport_id,bool>(src,true));
        current = src;
        for(auto p : pg->getVports()) {
            if(p != src)
                not_visited.insert(pair<vport_id, bool>(p, true));
        }
    }

    ~DFSIterator(){}

    DFSIterator operator++() {
        while(!not_visited.empty()){
            vport_id current_src = current;
            for(Edge<V,P,E> e : pg->getAdjList(current_src)){
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
            if(!path.empty())
                current = path.back();
            else{// get random new 'not visited' vport
                auto it = not_visited.begin();
                current = (*it).first;
                path.push_back(current);
                visited.insert(pair<vport_id,bool>(current,true));
                not_visited.erase(current);
                return (*this);
            }
        }
        current = END_VPORT;
        return *this;
    }

    DFSIterator operator++(int) const {
        DFSIterator old = DFSIterator(path,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    vport_id operator*()const {
        return current;
    }

//  vport_id const& operator*(){
//        return *current;
//    }

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
    DFSVertexIterator(vector<int> _path,PortGraph<V,P,E>* _pg,map<int,bool>& _visited,map<int,bool>& _not_visited,int _current){
        path = _path;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    DFSVertexIterator(int id){
        current = id;
    }

    DFSVertexIterator(PortGraph<V,P,E>* _pg,int src){
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

    ~DFSVertexIterator(){}

    DFSVertexIterator operator++() {
        while(!not_visited.empty()){
            int current_src = current;
            //done
            for(int dst : pg->getVertexAdjList(current_src)){
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
            if(!path.empty())
                current = path.back();
            else{// get random new 'not visited' vertex
                auto it = not_visited.begin();
                current = (*it).first;
                path.push_back(current);
                visited.insert(pair<int,bool>(current,true));
                not_visited.erase(current);
                return (*this);
            }
        }
        current = END_VERTEX;
        return *this;
    }

    DFSVertexIterator operator++(int) const {
        DFSVertexIterator old = DFSVertexIterator(path,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    int operator*()const {
        return current;
    }

//  vport_id const& operator*(){
//        return *current;
//    }

    bool operator== (const PGVertexIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGVertexIterator& that)const {return !this->operator == (that) ;}
};

/// BFS BY VPORT
template <class V , class P , class E >
class BFSIterator : public PGVportIterator{
private:
    vector<vport_id> queue;
    PortGraph<V,P,E>* pg;
    map<vport_id,bool> visited;
    map<vport_id,bool> not_visited;
    //vport_id current;
    BFSIterator(vector<vport_id> _queue,PortGraph<V,P,E>* _pg,map<vport_id,bool>& _visited,map<vport_id,bool>& _not_visited,vport_id _current){
        queue = _queue;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    BFSIterator(vport_id id){
        current = id;
    }

    ~BFSIterator(){}

    BFSIterator(PortGraph<V,P,E>* _pg,vport_id src){
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
            current = vport_id(END_VPORT);
            return *this;
        }
        vport_id current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for(Edge<V,P,E> e : pg->getAdjList(current_src)){
            vport_id dst = e.EdgeId().second;
            if(visited.find(dst) == visited.end()){
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<vport_id,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if(!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else{// get new random 'not visited' vport
            if(not_visited.empty()){
                current = vport_id(END_VPORT);
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
            current = vport_id(END_VPORT);
            return *this;
        }
        vport_id current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for(Edge<V,P,E> e : pg->getAdjList(current_src)){
            vport_id dst = e.EdgeId().second;
            if(visited.find(dst) == visited.end()){
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<vport_id,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if(!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else{ // return end
            current = vport_id(END_VPORT);
            return *this;
        }
        return (*this);
    }

    BFSIterator operator++(int) const {
        BFSIterator old = BFSIterator(queue,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    vport_id operator*()const {
        return current;
    }

//  vport_id const& operator*(){
//        return *current;
//    }

    bool operator== (const PGVportIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGVportIterator& that)const {return !this->operator==(that);}
};

/// BFS BY VERTEX
template <class V , class P , class E >
class BFSVertexIterator : public PGVertexIterator{
private:
    vector<int> queue;
    PortGraph<V,P,E>* pg;
    map<int,bool> visited;
    map<int,bool> not_visited;
    //int current;
    BFSVertexIterator(vector<int> _queue,PortGraph<V,P,E>* _pg,map<int,bool>& _visited,map<int,bool>& _not_visited,int _current){
        queue = _queue;
        pg = _pg;
        visited = _visited;
        not_visited = _not_visited;
        current = _current;
    }
public:

    BFSVertexIterator(int id){
        current = id;
    }

    ~BFSVertexIterator(){}

    BFSVertexIterator(PortGraph<V,P,E>* _pg,int src){
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
            current = END_VERTEX;
            return *this;
        }
        int current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for(int dst : pg->getVertexAdjList(current_src)){
            if(visited.find(dst) == visited.end()){
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<int,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if(!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else{// get new random 'not visited' vertex
            if(not_visited.empty()){
                current = END_VERTEX;
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
            current = END_VERTEX;
            return *this;
        }
        int current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for(int  dst : pg->getVertexAdjList(current_src)){
            if(visited.find(dst) == visited.end()){
                assert(not_visited.find(dst) != not_visited.end());
                visited.insert(pair<int,bool>(dst,true));
                not_visited.erase(dst);
                queue.push_back(dst);
            }
        }
        if(!queue.empty()) {
            current = queue.front();
            return (*this);
        }
        else{ // return end
            current = END_VERTEX;
            return *this;
        }
        return (*this);
    }

    BFSVertexIterator operator++(int) const {
        BFSVertexIterator old = BFSVertexIterator(queue,pg,visited,not_visited,current);
        this->operator++();
        return old;
    }

    int operator*()const {
        return current;
    }

//  vport_id const& operator*(){
//        return *current;
//    }

    bool operator== (const PGVertexIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGVertexIterator& that)const {return !this->operator == (that);}
};

#endif //PROJECT1_PORTGRAPH_H

/* TODO:
 *  topological_sort -- DONE
 *  strongly_connected_components -- DONE
 *  transpose_graph -- Done
 *  min_spanning_tree -- Done
 *  {DFS | BFS} iterator (vports/vertices) -- DONE
 *  is_bipartite -- DONE
 *  induced_graph by (vports/vertices/edges) -- DONE
 *  shortestpath(weight_function) -- DONE
 *  findPathCost(vport, vport, cost_function) -- DONE
 *  is_reachable(vport source, vport dest) -- DONE
 *  findClique (vports/vertices) -- DONE
 *  isSubGraph -- DONE
 * is_reachable(vertex source, vertex dest) -- DONE
 * max_flow -- DONE

FERAS
min_cut() --

MARIO
testing -- PIW
error handling -- PIW
check validations -- PIW

4 later
make_connected // maybe with min edges --
*/