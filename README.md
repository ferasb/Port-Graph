# Port Graph
 A library for Port Graph.
 A Port Graph is an extension of a regular graph, in which each vertex have a set of ports, and each edge connects between a pair of a vertex and a port (denoted vport) to another vport.
 This is a basic library for Directed Port Graphs, and includes basic contruction operation (see code for more info):
 1. Building an empty PortGraph
 2. Building a PortGraph from a list of vertices, ports, edges and attributes.
 3. Adding a vertex with ports.
 4. Adding a port to a vertex.
 4. Adding an edge.

 It also includes some graph algorithims that were extented to Port Graph:
 1. Reverse Graph: Reverses each edge
 2. Topological Sort: If the Port Graph is DAG then this function outputs a toplogical sort of the Port Graph.
 3. findMST: Finds a minimum spanning tree and its weight.
 4. isBipartite: Checks if the Port Graph is bipartite.
 5. inducedGraph: Outputs a sub-PortGraph of the main PortGraph according to a predicator on the vertices/vports/edges.
 6. isReachable: Checks if two vertices/vports are reachable.
 7. ShortestPath: Calculates the shortest path between two vports according to a weight function.
 7. ShortestPathWeight: Calculates the shortest path weight between two vports according to a weight function.
 8. maxFlow: Calculates the max flow in a network between two vports according to a capacity function.
 9. findVportClique/findVertexClique: Tries to find a clique of a given size in the Port Graph.
 10. isSubGraph: checks if the PortGraph is a sub-portGraph of another PortGraph.
 11. BFS/DFS iterators.
