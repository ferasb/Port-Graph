#ifndef PORJECT_OFF_UTILITIES_H
#define PORJECT_OFF_UTILITIES_H


#include <utility>
#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t
#include <assert.h> // For assert
#include <algorithm> // For sort
#include <float.h>  // For DBL_MAX
#include <queue>   // For priority_queue
#include <limits>  // For numeric_limits<int>::max()
using namespace std;

// Forward declarations
template <class V = int, class P = int>
class Vertex;

template <class P = int>
class Port;

template <class V = int, class P = int, class E = int>
class Edge;

#define vport pair<Vertex<V, P>, Port<P>>

typedef pair<int, int> vport_id;

typedef pair<vport_id, vport_id> edge_id;

typedef pair<vport_id, vport_id> vport_pair_id;

typedef int vertex_id;

typedef double (*WeightFunction)(edge_id);

typedef int (*CapacityFunction)(edge_id);

//For Vport DFS/BFS Iterators
class PGVportIterator{
    const vport_id END_VPORT = vport_id(-1,-1);
public:
    vport_id current;
    PGVportIterator()= default;
    PGVportIterator(vport_id src){
        current = src;
    }
    vport_id operator*()const {
        return current;
    }

};

//For Vertex DFS/BFS Iterators
class PGVertexIterator{
    const int END_VERTEX = -1;
public:
    int current;
    PGVertexIterator()= default;
    PGVertexIterator(int src){
        current = src;
    }
    int operator*()const {
        return current;
    }

};

vector<string> split (string& s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

#define VerticesAttributes vector<V>
#define PortsAttributes vector<P>
#define EdgesAttributes vector<E>

#endif //PORJECT_OFF_UTILITIES_H
