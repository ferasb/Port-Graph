//
// Created by mario barbara
//
#ifndef PORJECT_OFF_EX_BIPARTITE_H
#define PORJECT_OFF_EX_BIPARTITE_H

/**** MaxFlow Example ****/

#include <iostream>
#include "PortGraph.h"

int g(edge_id e) {
    vector<vport_id> ids;
    for(int i = 0; i < 6;i++)
        ids.push_back(vport_id(i,0));

    if (e == edge_id(ids[0], ids[2]))
        return 2;
    else if (e == edge_id(ids[0], ids[1]))
        return 2;
    else if (e == edge_id(ids[1], ids[3]))
        return 2;
    else if (e == edge_id(ids[2], ids[4]))
        return 2;
    else if (e == edge_id(ids[3], ids[5]))
        return 1;
    else if (e == edge_id(ids[4], ids[5]))
        return 1;
}

void PG_SHORTESTPATHS() {
    ostringstream s;
    s << "Max Flow Example" << endl;
    s << "Attributes : <int ,int ,int>" << endl ;
    fprintf(stdout, s.str().c_str());
    vector<vport_id> ids;
    vector<int> ports_num;
    for(int i = 0; i < 6;i++){
        ids.push_back(vport_id(i,0));
        ports_num.push_back(1);
    }
    vector<edge_id> edges_list({edge_id(ids[0], ids[2]),
                                edge_id(ids[0], ids[1]),
                                edge_id(ids[1], ids[3]),
                                edge_id(ids[2], ids[4]),
                                edge_id(ids[3], ids[5]),
                                edge_id(ids[4], ids[5])});

    PortGraph<int, int, int> pg = PortGraph<int, int, int>(6, ports_num, edges_list);

    vport_id src = ids[0];
    vport_id dst = ids[5];
    s.str("");
    int flow = pg.maxFlow(g, src, dst);

    s << "Flow from vport (0,0) to vport (5,0) is " << flow << "." << endl;
    fprintf(stdout, s.str().c_str());
}

#endif

