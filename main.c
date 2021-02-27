#include "PortGraph.h"

int main()
{
    vport_id id1 = vport_id(0, 0);
    vport_id id2 = vport_id(1, 1);
    vector<int> ports_num({1, 2, 3, 4 ,5});
    vector<edge_id> edges_list({edge_id(id1, id2)});
    PortGraph<int, int, int> pg = PortGraph<int, int, int>(5, ports_num, edges_list);
    pg.Print();
    return 0;
}