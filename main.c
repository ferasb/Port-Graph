#include <iostream>
#include "PortGraph.h"

// Basic test - without Attr
void test1(){
    ostringstream s;
    s << "TEST 1" << endl;
    s << "testing - int , int ,int" << endl ;
    fprintf(stdout, s.str().c_str());
    vport_id id1 = vport_id(0, 0);
    vport_id id2 = vport_id(0, 1);
    vport_id id3 = vport_id(1, 0);
    vport_id id4 = vport_id(1, 1);
    vport_id id5 = vport_id(2, 0);
    vport_id id6 = vport_id(2, 1);
    vport_id id7 = vport_id(3, 0);
    vport_id id8 = vport_id(3, 1);
    vector<int> ports_num({2,2,2,2,2});
    vector<edge_id> edges_list({edge_id(id1, id3),
                                edge_id(id1, id5),
                                edge_id(id4, id8),
                                edge_id(id6, id7)});
    PortGraph<int, int, int> pg = PortGraph<int, int, int>(4, ports_num, edges_list);
    // ----
    vport_id id9 = vport_id(4, 0);
    vport_id id10 = vport_id(4, 1);
    pg.addVertix(4,2);
    pg.Print();
    pg.addEdge({edge_id(id1,id9)});
    pg.addEdge({edge_id(id3,id9)});
    pg.addEdge({edge_id(id5,id9)});
    pg.addEdge({edge_id(id7,id9)});
    pg.PrintEdges();
    //---
}

// Basic test - with Attr
void test2(){
    ostringstream s;
    s << "TEST 2" << endl;
    s << "testing - int , int ,int : Attr" << endl ;
    fprintf(stdout, s.str().c_str());
    vport_id id1 = vport_id(0, 0);
    vport_id id2 = vport_id(0, 1);
    vport_id id3 = vport_id(1, 0);
    vport_id id4 = vport_id(1, 1);
    vport_id id5 = vport_id(2, 0);
    vport_id id6 = vport_id(2, 1);
    vport_id id7 = vport_id(3, 0);
    vport_id id8 = vport_id(3, 1);
    vector<int> ports_num({2,2,2,2,2});
    vector<edge_id> edges_list({edge_id(id1, id3),
                                edge_id(id1, id5),
                                edge_id(id4, id8),
                                edge_id(id6, id7)});


    vector<int> verAttr = {0,10,20,30};
    vector<vector<int>> portAttr = {{00,10},{100,110},{200,210},{300,310}};
    vector<std::string> edgeAttr = {"0->1","0->2","1->3","2->3"};
    PortGraph<int, int, std::string> pg = PortGraph<int, int, std::string>(4, ports_num, edges_list,verAttr,portAttr,edgeAttr);
    // ----
    vport_id id9 = vport_id(4, 0);
    vport_id id10 = vport_id(4, 1);
    pg.addVertix(4,1,40,{400});
    pg.addEdge(edge_id(id1,id9),"1->4");
    pg.Print();
    pg.PrintEdges();

}


int main()
{
    //test1();
    //test2();
    return 0;
}