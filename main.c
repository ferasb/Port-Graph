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
    fprintf(stderr, s.str().c_str());
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
    // DFS TEST
    for(DFSIterator<int,int,int> it1(&pg,ids[0]);it1 != pg.end() ;++it1){
        vport_id curr = (*it1);
        ostringstream s1;
        s1<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s1.str().c_str());
    }
    s<<endl;
    fprintf(stderr,s.str().c_str());
    //---BFS Test
    for(BFSIterator<int,int,int> it2(&pg,ids[0]);it2 != pg.end() ;++it2){
        vport_id curr = (*it2);
        ostringstream s2;
        s2<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s2.str().c_str());
    }
}

// MST Test
// Link to graph : http://graphonline.ru/en/?graph=uLcpFMUCtmmKCgGl
double w(const edge_id id, const double e){
    return e;
}
void test4(){
    ostringstream s;
    s << "TEST 4" << endl;
    s << "testing MST - int , int ,int" << endl ;
    fprintf(stderr, s.str().c_str());
    vector<vport_id> ids;
    vector<int> ports_num;
    for(int i = 0; i < 12;i++){
        ids.push_back(vport_id(i,0));
        ports_num.push_back(1);
    }
    vector<edge_id> edges_list({edge_id(ids[0], ids[2]),
                                edge_id(ids[0], ids[4]),
                                edge_id(ids[0], ids[3]),
                                edge_id(ids[2], ids[5]),
                                edge_id(ids[2], ids[1]),
                                edge_id(ids[4], ids[7]),
                                edge_id(ids[4], ids[11]),
                                edge_id(ids[1], ids[9]),
                                edge_id(ids[1], ids[10]),
                                edge_id(ids[7], ids[10]),
                                edge_id(ids[9], ids[6]),
                                edge_id(ids[10], ids[6]),
                                edge_id(ids[11], ids[8])});

    vector<double> edgeAttr({0,0,0,5,1.5,1.5,1,2,1,0.5,1,1,2.5,0});
    PortGraph<int, int, double> pg = PortGraph<int, int, double >(12, ports_num, edges_list,vector<int>(),vector<vector<int>>(),edgeAttr);
    // ----
    // DFS TEST
    //pg.PrintEdges();
    for(DFSIterator<int,int,double > it1(&pg,ids[0]);it1 != pg.end() ;++it1){
        vport_id curr = (*it1);
        ostringstream s1;
        s1<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s1.str().c_str());
    }
    s<<endl;
    fprintf(stderr,s.str().c_str());
    //---BFS Test
    for(BFSIterator<int,int,double> it2(&pg,ids[0]);it2 != pg.end() ;++it2){
        vport_id curr = (*it2);
        ostringstream s2;
        s2<< curr.first << " " << curr.second<< endl;
        fprintf(stderr,s2.str().c_str());
    }
    //---MST Test
    assert(pg.Kruskal(w) == 13.5);
    //--- isBipartite Test
    assert(pg.isBipartite(ids[0])==true);
    
    //---
    //ADD edge mst=12 ,not 2 colored
    pg.addEdge(edge_id(ids[6], ids[0]),0);
    assert(pg.Kruskal(w) == 12);
    //--- isBipartite Test
    assert(pg.isBipartite(ids[0])==false);
}


void test5(){
    ostringstream s;
    s << "TEST 5" << endl;
    s << "testing isReachable - int , int ,int" << endl ;
    fprintf(stderr, s.str().c_str());
    vector<vport_id> ids;
    vector<int> ports_num;
    for(int i = 0; i < 12;i++){
        ids.push_back(vport_id(i,0));
        ports_num.push_back(1);
    }
    vector<edge_id> edges_list({edge_id(ids[0], ids[2]),
                                edge_id(ids[0], ids[4]),
                                edge_id(ids[0], ids[3]),
                                edge_id(ids[2], ids[5]),
                                edge_id(ids[2], ids[1]),
                                edge_id(ids[4], ids[7]),
                                edge_id(ids[4], ids[11]),
                                edge_id(ids[1], ids[9]),
                                edge_id(ids[1], ids[10]),
                                edge_id(ids[7], ids[10]),
                                edge_id(ids[9], ids[6]),
                                edge_id(ids[10], ids[6]),
                                edge_id(ids[11], ids[8])});

    vector<double> edgeAttr({0,0,0,5,1.5,1.5,1,2,1,0.5,1,1,2.5,0});
    PortGraph<int, int, double> pg = PortGraph<int, int, double >(12, ports_num, edges_list,vector<int>(),vector<vector<int>>(),edgeAttr);

    s.str("");
    bool reach00_50 = false; // true
    bool reach00_10 = false; // true
    bool reach80_00 = false; // false
    reach00_50 = pg.is_reachable(vport_id(0,0), vport_id(5,0));
    reach00_10 = pg.is_reachable(vport_id(0,0), vport_id(1,0));
    reach80_00 = pg.is_reachable(vport_id(8,0), vport_id(0,0));

    s << "reach00_50 is " << (reach00_50? "true" : "false") << endl;
    s << "reach00_10 is " << (reach00_10? "true" : "false") << endl;
    s << "reach80_00 is " << (reach80_00? "true" : "false") << endl;
    fprintf(stderr, s.str().c_str());

}

double f(edge_id id) { return 1;}

void test6(){
    ostringstream s;
    s << "TEST 6" << endl;
    s << "testing shortest paths - int , int ,int" << endl ;
    fprintf(stderr, s.str().c_str());
    vector<vport_id> ids;
    vector<int> ports_num;
    for(int i = 0; i < 12;i++){
        ids.push_back(vport_id(i,0));
        ports_num.push_back(1);
    }
    vector<edge_id> edges_list({edge_id(ids[0], ids[2]),
                                edge_id(ids[0], ids[4]),
                                edge_id(ids[0], ids[3]),
                                edge_id(ids[2], ids[5]),
                                edge_id(ids[2], ids[1]),
                                edge_id(ids[4], ids[7]),
                                edge_id(ids[4], ids[11]),
                                edge_id(ids[1], ids[9]),
                                edge_id(ids[1], ids[10]),
                                edge_id(ids[7], ids[10]),
                                edge_id(ids[9], ids[6]),
                                edge_id(ids[10], ids[6]),
                                edge_id(ids[11], ids[8]),
                                edge_id(ids[0], ids[8])});

    vector<double> edgeAttr({0,0,0,5,1.5,1.5,1,2,1,0.5,1,1,2.5,0});
    PortGraph<int, int, double> pg = PortGraph<int, int, double >(12, ports_num, edges_list,vector<int>(),vector<vector<int>>(),edgeAttr);

    WeightFunction wf = f;
    vport_id src = ids[0];
    vport_id dst = ids[6];
    s.str("");
    auto pth = pg.shortestPath(wf, src, dst);
    double weight = pg.shortestPathWeight(wf, src, dst);

    s << "shorest path weight from vport 00 to vport 60 " << weight << "." << endl;
    s << "shortest path is: ";
    for (auto& id: pth) {
        s << "(" << id.first.first << ", " << id.first.second << ") " << "-- (" << id.second.first << ", " << id.second.second << "), ";
    } 
    fprintf(stderr, s.str().c_str());

}

int main()
{
    //test1();
    //test2();
    //test3();
    //test4();
    //test5();
    test6();

    return 0;

}