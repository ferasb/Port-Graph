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

#define vport pair<Vertex<V, P>, Port<P> >

typedef pair<int, int> vport_id;

#define END_VPORT vport_id(-1,-1)

#define END_VERTEX -1

typedef pair<vport_id, vport_id> edge_id;

typedef pair<vport_id, vport_id> vport_pair_id;


typedef double (*WeightFunction)(edge_id);

typedef int (*CapacityFunction)(edge_id);

//For Vport DFS/BFS Iterators
class PGVportIterator{
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

//For Vport DFS/BFS Iterators
class PGVertexIterator{
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


#define VerticesAttributes vector<V>
#define PortsAttributes vector<P>
#define EdgesAttributes vector<E>
#endif //PROJECT1_UTILITIES_H
