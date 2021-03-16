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
    // ---- Insert Test
    vport_id id9 = vport_id(4, 0);
    vport_id id10 = vport_id(4, 1);
    pg.addVertix(4,2);
    pg.Print();
    pg.addEdge({edge_id(id1,id9)});
    pg.addEdge({edge_id(id3,id9)});
    pg.addEdge({edge_id(id5,id9)});
    pg.addEdge({edge_id(id7,id9)});
    pg.PrintEdges();
    pg.PrintEdges();
    //--- DFS Test
    DFSIterator<int,int,int> it1(&pg,id1);
    for(;it1 != pg.end() ;++it1){
        vport_id curr = (*it1);
        ostringstream s1;
        s1<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s1.str().c_str());
    }
    //---BFS Test
    BFSIterator<int,int,int> it2(&pg,id1);
    for(;it2 != pg.end() ;++it2){
        vport_id curr = (*it2);
        ostringstream s2;
        s2<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s2.str().c_str());
    }

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

// Tree Test - without Attr
// Link to graph : http://graphonline.ru/en/?graph=TmJTylolZRuHGGkP
void test3(){
    ostringstream s;
    s << "TEST 3" << endl;
    s << "testing TREE - int , int ,int" << endl ;
    fprintf(stdout, s.str().c_str());
    vector<vport_id> ids;
    vector<int> ports_num;
    for(int i = 0; i < 12;i++){
        ids.push_back(vport_id(i,0));
        ports_num.push_back(1);
    }
    vector<edge_id> edges_list({edge_id(ids[0], ids[1]),
                                edge_id(ids[0], ids[2]),
                                edge_id(ids[1], ids[3]),
                                edge_id(ids[1], ids[4]),
                                edge_id(ids[2], ids[5]),
                                edge_id(ids[2], ids[6]),
                                edge_id(ids[3], ids[7]),
                                edge_id(ids[3], ids[8]),
                                edge_id(ids[4], ids[9]),
                                edge_id(ids[5], ids[10]),
                                edge_id(ids[6], ids[11])});

    PortGraph<int, int, int> pg = PortGraph<int, int, int>(12, ports_num, edges_list);
    // ----
    pg.Print();
    pg.PrintEdges();
    // DFS TEST
    for(DFSIterator<int,int,int> it1(&pg,ids[0]);it1 != pg.end() ;++it1){
        vport_id curr = (*it1);
        ostringstream s1;
        s1<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s1.str().c_str());
    }
    s<<endl;
    fprintf(stdout,s.str().c_str());
    //---BFS Test
    for(BFSIterator<int,int,int> it2(&pg,ids[0]);it2 != pg.end() ;++it2){
        vport_id curr = (*it2);
        ostringstream s2;
        s2<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s2.str().c_str());
    }
}

int main()
{
    //test1();
    //test2();
    test3();
    return 0;
}