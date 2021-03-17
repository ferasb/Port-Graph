#ifndef PROJECT1_PORTGRAPH_H
#define PROJECT1_PORTGRAPH_H
#include "Utilities.h"
using namespace std;

template <class P>
#define PortSet set<Port<P>, cmpPort<P>>
struct cmpPort{
    bool operator () (const Port<P>& a,const Port<P>& b){
        return a.getPortId() < b.getPortId() ;
    }
};

template <class P>
class Port{
public:

    explicit Port(int portId = -1, P attr = P()) : port_id(portId) , attribute(attr) {}

    int getPortId() const { return port_id; }

    P getAttribute()  { return attribute; }

    bool operator==(const Port& p) { return p.getPortId() == port_id; }

    bool operator<(const Port& p){
        return this->port_id < p.getPortId();
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

template <class V, class P>
class Vertex {
public:
    Vertex () {}
    Vertex(int vertex_id , int ports_num , V vertex_attr = V(), PortsAttributes ports_attr = PortsAttributes())
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

    bool operator==(Vertex v) { return v.vertexId() == vertex_id; }

    bool operator<(Vertex v){
        return  v.vertexId() > vertex_id;
    }

//    Vertex<V,P>& operator=(Vertex<V,P> const &a){
//        this->vertex_id = a.vertexId();
//        this->attribute = a.getAttribute();
//        this->ports = a.getPorts();

//    }

    void print() {
        ostringstream s;
        s << "Vertrex " << vertex_id << " with attribute " << attribute << ", with " << ports.size() << " ports." << endl;
        s << "Vertrex " << vertex_id << " with Ports attribute : " << endl ;
        fprintf(stderr, s.str().c_str());
        for (auto p : ports)
            p.print();
    }

private:
    int vertex_id;
    V attribute;
    PortSet ports;
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

    edge_id getEdgeId() const { return id; }

    bool operator==(Edge e) { return e.id == id; }

    E getAttribute() { return attribute; }

     void print() {
         ostringstream s;
         s << "Edge (" << id.first.first << ", " << id.first.second << ") " << "-- (" << id.second.first << ", " << id.second.second << ") " ;
         s << "with attribute " << attribute << endl ;
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
        return a.getEdgeId() < b.getEdgeId() ;
    }
};

template <class V , class P , class E >
class DFSIterator;

template <class V , class P , class E >
class BFSIterator;

template <class V = int, class P = int, class E = int>
class PortGraph {
private:
    // maps vport vp (vertix and port) to all the edges that has an vp as a source
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adjacency_list;
    // maps vport vp (vertix and port) to all the transpose edges that has an vp as a source
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> transpose_adjacency_list;
    // flag for the transpose Graph
    bool is_transpose ;
    // maps vport_id (vertix and port) to its vport
    map<vport_id, vport, cmpVport> vport_map;
    // params for SCC
    const int UNSEEN = -1;
    const int SEEN = 1;
    // for Algo's in case of transpose
    map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& AdjacencyList(){
        return is_transpose ? transpose_adjacency_list : adjacency_list ;
    }

public:
    // Build PortGraph with
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
        is_transpose = !is_transpose;
    }
    void addVertix(int vertex_id, int ports_num, V attr = V(), PortsAttributes ports_attr = PortsAttributes())
    {
        for (int i = 0; i < ports_num; ++i) {
            adjacency_list[vport_id(vertex_id, i)] = set<Edge<V, P, E>, cmpEdge<V,P,E>>();
            // transpose Graph
            transpose_adjacency_list[vport_id(vertex_id, i)] = set<Edge<V, P, E>, cmpEdge<V,P,E>>();
            P portAttr =  ports_attr.empty() ? P() : ports_attr[i] ;
            vport_map[vport_id(vertex_id, i)]
                    = vport(Vertex<V, P>(vertex_id, ports_num, attr, ports_attr), Port<P>(i,portAttr));
        }
    }

    void addEdge(edge_id id, E attr = E())
    {
        /// debug
        Edge<V, P, E> e(Edge<V, P, E>(vport_map[id.first], vport_map[id.second], attr)) ;
        set<Edge<V, P, E>, cmpEdge<V,P,E>> sety = adjacency_list[id.first] ;
        sety.insert(Edge<V, P, E>(vport_map[id.first], vport_map[id.second], attr)) ;
        // add edge
        adjacency_list[id.first].insert(Edge<V, P, E>(vport_map[id.first], vport_map[id.second], attr) );
        // add transpose Edge
        transpose_adjacency_list[id.second].insert(Edge<V, P, E>(vport_map[id.second], vport_map[id.first], attr));
        ///  debug
        set<Edge<V, P, E>, cmpEdge<V,P,E>> setx = adjacency_list[id.first] ;
    }

    //Print the vertcies and ports per vertex
    void Print()
    {
        for (auto x : vport_map) {
            Vertex<V, P>& v = x.second.first;
            v.print();
        }
    }

    //Print the edges
    void PrintEdges() {

        for (typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adjacency_list.begin(); it != adjacency_list.end(); ++it)
            for(auto e : (*it).second)
                e.print();
    }

/********** Topological Sort **********/

// input: directed graph (g[u] contains the neighbors of u, nodes are named 0,1,...,|V|-1).
// output: is g a DAG (return value), a topological ordering of g (order).
// comment: order is valid only if g is a DAG.
    bool topological_sort(){
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport> adj_list = AdjacencyList();
        // compute indegree of all nodes
        map<vport_id,int> vport_indegree = map<vport_id,int>();
        for(typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it) {
            set<Edge<V, P, E>> edge_set = (*it).second;
            vport_id id = (*it).first;
            for(auto edge : edge_set) {
                vport_indegree[edge.getEdgeId().second]++ ;
            }
        }

        // order sources first
        vector<vport_id> vport_order = vector<vport_id>();
        for(typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it) {
            vport_id id = (*it).first;
            if (vport_indegree[id] == 0)
                vport_order.push_back(id);
        }
        // go over the ordered nodes and remove outgoing edges,
        // add new sources to the ordering
        for(auto id : vport_order){
            for (auto edge : adj_list[id]) {
                vport_id neighbor = edge.getEdgeId().second;
                vport_indegree[neighbor] -- ;
                if (vport_indegree[neighbor] == 0)
                    vport_order.push_back(neighbor);
            }
        }
        // check -
        return vport_order.size() == vport_map.size();
    }

/********** Strongly Connected Components **********/

    void KosarajuDFS (const map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list ,vport_id id, map<vport_id,vport_id>& S, map<vport_id,int>& colorMap, int color){

        colorMap[id] = color;
        for(auto edge : adj_list[id]){
            vport_id neighbor = edge.getEdgeId().second;
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
        for(typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it){
            vport_id id = (*it).first;
            seen[id] = UNSEEN;
            components[id] = UNSEEN;
        }
        for(typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it){
            vport_id id = (*it).first;
            if(seen[id] == UNSEEN)
                KosarajuDFS(adj_list , id,postOrder,seen,SEEN);
        }
        // reverse iterator
        // second pass: explore the SCCs based on first pass result
        int numSCC = 0;
        for(typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::reverse_iterator rit = adj_list.rbegin(); rit != adj_list.rend(); ++rit){
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

    struct unionfindPG  {
        map<vport_id ,int> rank;
        map<vport_id ,vport_id >parent;

        explicit unionfindPG (map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list ) {
            for (typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it){
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
    double Kruskal(double(*weightFunc)(const edge_id ,const E attr)){
        int n = AdjacencyList().size();
        // (weight , edge_id)
        vector<pair<double,edge_id>> edges;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for (typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it){
            set<Edge<V, P, E>,cmpEdge<V,P,E>> s_edge = (*it).second;
            for(auto e : s_edge)
                edges.push_back(pair<double,edge_id>(weightFunc(e.getEdgeId(),e.getAttribute()) ,e.getEdgeId()));
        }
        // with weightFunc
        sort(edges.begin(), edges.end());
        double mst_cost = 0;
        unionfindPG components(adj_list);
        for (pair<double,edge_id> e : edges) {
            if (components.find(e.second.first) != components.find(e.second.second)) {
                // W function on edges
                mst_cost += e.first;
                components.unite(e.second.first, e.second.second);
            }
        }
        return mst_cost;
    }

/********** BFS/DFS **********/
    PGIterator end(){
        return PGIterator(vport_id(-1,-1));
    }

    vector<vport_id> getVports(){
        vector<vport_id> res;
        for(typename map<vport_id, vport,cmpVport>::iterator it = vport_map.begin(); it != vport_map.end(); ++it) {
            vport_id p = (*it).first;
            res.push_back(p);
        }
        return res;
    }

    set<Edge<V,P,E>,cmpEdge<V,P,E>> getAdjList(vport_id id){
        return AdjacencyList()[id];
    }

/********** Bipartite **********/
// This function returns true if graph is Bipartite, else false
    bool isNeighbers(const vport_id& src, const vport_id& dst){
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for(const Edge<V,P,E>& e : adj_list[src]){
            if(e.getEdgeId().second == dst)
                return true;
        }
        return false;
    }

    bool isBipartite(vport_id src){
        // if src not in graph ...
        map<vport_id,int> colorArr;
        map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>& adj_list = AdjacencyList();
        for (typename map<vport_id, set<Edge<V, P, E>, cmpEdge<V,P,E>>, cmpVport>::iterator it = adj_list.begin(); it != adj_list.end(); ++it) {
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
            for (typename map<vport_id, vport, cmpVport>::iterator it = vport_map.begin(); it != vport_map.end(); ++it){
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

    /* TODO:
     *  topological_sort -- DONE
     *  strongly_connected_components -- DONE
     *  transpose_graph -- Done
     *  min_spanning_tree -- Done
     *  is_reachable(vport source, vport dest) -- Done
     *  {DFS | BFS} we should implement an {DFS | BFS} iterator -- DONE
     *  is_bipartite() -- DONE

    shortestpath(weight_function) --
    max_flow --
    min_cut() --
    colouring --
    induced_graph(vports/vertices/edges) // or filter --
    is_subgraph(subgraph) --
    netlistToPortGraph // maybe constructor --
    findPathCost(vport, vport, cost_function) --
    make_connected // maybe with min edges --
    */
};


template <class V , class P , class E >
class DFSIterator : public PGIterator {
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
                vport_id dst = e.getEdgeId().second;
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
        current = END;
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

    bool operator== (const PGIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGIterator& that)const {return !this->operator==(that);}
};

template <class V , class P , class E >
class BFSIterator : public PGIterator{
private:
    vector<vport_id> queue;
    PortGraph<V,P,E>* pg;
    map<vport_id,bool> visited;
    map<vport_id,bool> not_visited;
    vport_id current;
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
            current = vport_id(END);
            return *this;
        }
        vport_id current_src = queue.front();
        queue.erase(queue.begin());
        // add the next layer
        for(Edge<V,P,E> e : pg->getAdjList(current_src)){
            vport_id dst = e.getEdgeId().second;
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
                current = vport_id(END);
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

    bool operator== (const PGIterator& that)const {return (this->current) == (*that) ;}
    bool operator!= (const PGIterator& that)const {return !this->operator==(that);}
};

#endif //PROJECT1_PORTGRAPH_H
