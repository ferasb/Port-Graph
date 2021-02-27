#include <utility>
#include <set>
#include <vector>

using namespace std;

template<class V = int, class P = int>
class Vertex;

template<class P = int>
class Port;

template<class V = int, class P = int, class E = int>
class Edge;

#define vport pair<Vertex<V, P>, Port<P>>; 

typedef pair<int, int> vport_id;

typedef pair<vport_id, vport_id> edge_id;

#define PortSet set<Port<P>>;

#define PortsAttributes vector<P>;