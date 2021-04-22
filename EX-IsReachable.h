//
// Created by mario barbara
//
#ifndef PORJECT_OFF_EX_BIPARTITE_H
#define PORJECT_OFF_EX_BIPARTITE_H

/**** IsReachable Example ****/

#include <iostream>
#include "PortGraph.h"

void PG_ISREACHABLE(){
    ostringstream s1;
    s1 << "Attributes : <int ,int ,int>" << endl ;
    s1 << "IsReachable Example " << endl;
    fprintf(stdout, s1.str().c_str());
    int n_vertex = 5;
    vport_id id00 = vport_id(0, 0);
    vport_id id01 = vport_id(0, 1);
    vport_id id10 = vport_id(1, 0);
    vport_id id11 = vport_id(1, 1);
    vport_id id20 = vport_id(2, 0);
    vport_id id21 = vport_id(2, 1);
    vport_id id30 = vport_id(3, 0);
    vport_id id31 = vport_id(3, 1);
    vport_id id32 = vport_id(3, 2);
    vport_id id40 = vport_id(4, 0);
    vector<int> ports_num({2,2,2,3,1});
    vector<edge_id> edges_list = {edge_id(id00, id10),
                                  edge_id(id00, id21),
                                  edge_id(id00, id31),
                                  edge_id(id01, id20),
                                  edge_id(id01, id11),
                                  edge_id(id20, id32),
                                  edge_id(id11, id20),
                                  edge_id(id11, id30),
                                  edge_id(id40, id21),
                                  edge_id(id40, id32),
    };

    /// Output : Is Reachable
    PortGraph<int, int, int> pg = PortGraph<int, int, int>(n_vertex, ports_num, edges_list);
    s1.str("");
    // Should be true
    s1 << "vport (2,1) is reachable from vport (0, 0): " << (pg.isReachable(id00, id21)? "True": "false") << endl;
    // Should be true
    s1 << "vport (3,2) is reachable from vport (0, 1): " << (pg.isReachable(id01, id32)? "True": "false") << endl;
    // Should be false
    s1 << "vport (3,0) is reachable from vport (0, 0): " << (pg.isReachable(id00, id30)? "True": "false") << endl;
    // Should be true after adding edge (1,0)->(3,0)
    pg.addEdge(edge_id(id10, id30));
    s1 << "adding edge (1,0)->(3,0)" << endl;
    s1 << "vport (3,0) is reachable from vport (0, 0): " << (pg.isReachable(id00, id30)? "True": "false") << endl;
    // Should be false
    s1 << "vport (4,0) is reachable from vport (2, 1): " << (pg.isReachable(id21, id40)? "True": "false") << endl;
    // Should be true
    s1 << "vertex 2 is reachable from vertex 0: " << (pg.isReachable(0, 2)? "True": "false") << endl;
    // Should be true
    s1 << "vertex 3 is reachable from vertex 4: " << (pg.isReachable(4, 3)? "True": "false") << endl;
    // Should be false
    s1 << "vertex 1 is reachable from vertex 2: " << (pg.isReachable(2, 1)? "True": "false") << endl;
    fprintf(stdout, s1.str().c_str());
}

#endif

