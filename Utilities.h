#include <utility>
#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <sstream>
using namespace std;

// Forward declarations
template <class V = int, class P = int>
class Vertex;

template <class P = int>
class Port;

template <class V = int, class P = int, class E = int>
class Edge;

#define vport pair<Vertex<V, P>, Port<P> >

typedef pair<int, int> vport_id;

typedef pair<vport_id, vport_id> edge_id;

#define PortSet set<Port<P> >

#define VerticesAttributes vector<V>
#define PortsAttributes vector<P>
#define EdgesAttributes vector<E>