//
// Created by mario barbara
//
#ifndef PORJECT_OFF_EX_BIPARTITE_H
#define PORJECT_OFF_EX_BIPARTITE_H

/**** Transpose Example ****/

#include <iostream>
#include "PortGraph.h"

void PG_TRANSPOSE(){
    ostringstream s;
    s << "Attributes : <int ,int ,int>" << endl ;
    s << "Transpose Example" << endl << endl;
    fprintf(stdout, s.str().c_str());
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
    edge_id e00_10 = edge_id(id00, id10);
    edge_id e00_21 = edge_id(id00, id21);
    edge_id e00_31 = edge_id(id00, id31);
    edge_id e01_20 = edge_id(id01, id20);
    edge_id e01_11 = edge_id(id01, id11);
    edge_id e20_32 = edge_id(id20, id32);
    edge_id e11_20 = edge_id(id11, id20);
    edge_id e40_21 = edge_id(id40, id21);
    edge_id e40_32 = edge_id(id40, id32);
    edge_id e21_32 = edge_id(id21, id32);
    edge_id e10_32 = edge_id(id01, id32);
    edge_id e21_31 = edge_id(id21, id31);


    vector<edge_id> edges_list = {e00_10,
                                  e00_21,
                                  e00_31,
                                  e01_20,
                                  e01_11,
                                  e20_32,
                                  e11_20,
                                  e40_21,
                                  e40_32,
                                  e21_32,
                                  e10_32,
                                  e21_31
    };

    PortGraph<int, int, int> pg = PortGraph<int, int, int>(n_vertex, ports_num, edges_list);
    s.str("");
    s << "Print graph before transpose" << endl;
    fprintf(stdout, s.str().c_str());
    pg.PrintEdgesIds();
    pg.transposeGraph();
    s.str("");
    s << endl << "Print graph after transpose" << endl;
    fprintf(stdout, s.str().c_str());
    pg.PrintEdgesIds();
}

#endif

