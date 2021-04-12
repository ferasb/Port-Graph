//// Created by mario barbara//#ifndef PORJECT_OFF_EX_BIPARTITE_H#define PORJECT_OFF_EX_BIPARTITE_H/**** BIPARTITE Example ****/#include <iostream>#include "PortGraph.h"/// EX1 : isBIPARTITE/// EX2 : isNotBIPARTITEvoid PG_BIPARTITE(){    ostringstream s1;    s1 << "Attributes : < Color ,Color ,String >" << endl ;    s1 << "BIPARTITE Example " << endl;    fprintf(stdout, s1.str().c_str());    typedef enum{blue,red,green,white,yellow} Color;    int n_vertex = 5;    vport_id id00 = vport_id(0, 0);    vport_id id01 = vport_id(0, 1);    vport_id id10 = vport_id(1, 0);    vport_id id20 = vport_id(2, 0);    vport_id id21 = vport_id(2, 1);    vport_id id22 = vport_id(2, 2);    vport_id id30 = vport_id(3, 0);    vport_id id31 = vport_id(3, 1);    vport_id id40 = vport_id(4, 0);    vector<int> ports_num({2,1,3,2,1});    vector<edge_id> edges_list = {edge_id(id00,id10),                                  edge_id(id01,id10),                                  edge_id(id20,id30),                                  edge_id(id21,id31),                                  edge_id(id31,id20),                                  edge_id(id22,id40),                                  };    vector<Color> verAttr = {white,white,white,white,white};    vector<vector<Color>> portAttr = {{blue,blue},{green},{green,green,green},{blue,blue},{blue}};    vector<std::string> edgeAttr = {            "V0:P0 -> V1:P0 : default",            "V0:P1 -> V1:P0 : default",            "V2:P0 -> V3:P0 : default",            "V2:P1 -> V3:P1 : default",            "V3:P1 -> V2:P0 : default",            "V2:P2 -> V4:P0 : default",    };    /// Output : is BIPARTITE    PortGraph<Color, Color, std::string> pg = PortGraph<Color, Color, std::string>(n_vertex, ports_num, edges_list,verAttr,portAttr,edgeAttr);    std::string res1 = pg.isBipartite() ? "is " : "isn't ";    ostringstream s_inner1;    s_inner1 << "The following Port Graph " << res1 << "BIPARTITE" << endl ;    fprintf(stdout, s_inner1.str().c_str());    pg.Print();    pg.PrintEdges();    // add new edges    ostringstream s2;    s2 << "After adding the following edges : " << endl ;    s2 << "V2:P1 -> V4:P0 with \"added\" attribute" << endl;    s2 << "V3:P1 -> V4:P0 with \"added\" attribute" << endl;    fprintf(stdout, s2.str().c_str());    pg.addEdge(edge_id(id21,id40),"V2:P1 -> V4:P0 : Added");    pg.addEdge(edge_id(id31,id40),"V3:P1 -> V4:P0 : Added");    std::string res2 = pg.isBipartite() ? "is " : "isn't ";    ostringstream s_inner2;    s_inner2 << "The following Port Graph " << res2 << "BIPARTITE" << endl ;    fprintf(stdout, s_inner2.str().c_str());    pg.Print();    pg.PrintEdges();}#endif