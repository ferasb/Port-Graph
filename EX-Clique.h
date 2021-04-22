//
// Created by mario barbara
//
#ifndef PORJECT_OFF_EX_BIPARTITE_H
#define PORJECT_OFF_EX_BIPARTITE_H

/**** Clique Example ****/

#include <iostream>
#include "PortGraph.h"

void PG_CLIQUE(){
    ostringstream s;
    s << "Attributes : <int ,int ,int>" << endl ;
    s << "Clique Example" << endl << endl;
    fprintf(stderr, s.str().c_str());
    vector<vport_id> ids;
    vector<int> ports_num;
    for(int i = 0; i < 4;i++){
        ids.push_back(vport_id(i,0));
        ports_num.push_back(2);
    }
    vector<edge_id> edges_list({edge_id(ids[0], ids[1]),
                                edge_id(ids[0], ids[2]),
                                edge_id(ids[0], ids[3]),
                                edge_id(ids[1], ids[0]),
                                edge_id(ids[1], ids[2]),
                                edge_id(ids[1], ids[3]),
                                edge_id(ids[2], ids[0]),
                                edge_id(ids[2], ids[1]),
                                edge_id(ids[2], ids[3]),
                                edge_id(ids[3], ids[0]),
                                edge_id(ids[3], ids[1]),
                                edge_id(ids[3], ids[2]),
                                edge_id(ids[4], ids[0]),
                                edge_id(ids[4], ids[2]),
                                edge_id(ids[4], ids[3])});

    PortGraph<int, int, int> pg1 = PortGraph<int, int, int>(5, ports_num, edges_list);
    auto res1 = pg1.findVportClique(3);
    s.str("");
    s << "Clique of size 3." << endl;
    fprintf(stdout, s.str().c_str());
    res1.PrintEdgesIds();
    auto res2 = pg1.findVportClique(4);
    s.str("");
    s << endl << "Clique of size 4." << endl;
    fprintf(stdout, s.str().c_str());
    res2.PrintEdgesIds();                           
}

#endif

