//
// Created by mario barbara on 08/03/2021.
//

#ifndef PROJECT1_UTILITIES_H
#define PROJECT1_UTILITIES_H
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t
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

#define inf  std::numeric_limits<int>::infinity()

#define END vport_id(-1,-1)

typedef pair<vport_id, vport_id> edge_id;

//For DFS/BFS Iterators
class PGIterator{
public:
    vport_id current;
    PGIterator()= default;
    PGIterator(vport_id src){
        current = src;
    }
    vport_id operator*()const {
        return current;
    }

};

#define VerticesAttributes vector<V>
#define PortsAttributes vector<P>
#define EdgesAttributes vector<E>
#endif //PROJECT1_UTILITIES_H
