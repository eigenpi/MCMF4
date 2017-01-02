#include <assert.h>
#include <vector>
#include <deque>

#include "mcmf4.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//
// MCMF_EDMONDS
//
////////////////////////////////////////////////////////////////////////////////


// Min cost max flow * (Edmonds-Karp relabelling + fast heap Dijkstra)
// Takes a directed graph where each edge has a capacity ('cap') and a 
// cost per unit of flow ('cost') and returns a maximum flow network
// of minimal cost ('fcost') from s to t. USE mcmf3.cpp FOR DENSE GRAPHS!
//
// PARAMETERS:
//      - cap (global): adjacency matrix where cap[u][v] is the capacity
//          of the edge u->v. cap[u][v] is 0 for non-existent edges.
//      - cost (global): a matrix where cost[u][v] is the cost per unit
//          of flow along the edge u->v. If cap[u][v] == 0, cost[u][v] is
//          ignored. ALL COSTS MUST BE NON-NEGATIVE!
//      - n: the number of vertices ([0, n-1] are considered as vertices).
//      - s: source vertex.
//      - t: sink.
// RETURNS:
//      - the flow
//      - the total cost through 'fcost'
//      - fnet contains the flow network. Careful: both fnet[u][v] and
//          fnet[v][u] could be positive. Take the difference.
// COMPLEXITY:
//      - Worst case: O(m*log(m)*flow  <?  n*m*log(m)*fcost)
// FIELD TESTING:
//      - Valladolid 10594: Data Flow
// REFERENCE:
//      Edmonds, J., Karp, R.  "Theoretical Improvements in Algorithmic
//          Efficieincy for Network Flow Problems".
//      This is a slight improvement of Frank Chu's implementation.


bool MCMF_EDMONDS::dijkstra( int n, int s, int t)
{
    CLEAR_MEM_1D( _d, HALF_INFINITY, _NN); // Note: HALF_INFINITY was 0x3F;
    CLEAR_MEM_1D( _par, -1, _NN);
    CLEAR_MEM_1D( _inq, -1, _NN);
    //for ( int i = 0; i < n; i++ ) _d[i] = HALF_INFINITY, _par[i] = -1;
    _d[s] = _qs = 0;
    _inq[_q[_qs++] = s] = 0;
    _par[s] = n;

    while ( _qs) {
        // get the minimum from _q and bubble down
        int u = _q[0];
        _inq[u] = -1;
        _q[0] = _q[--_qs];
        if ( _qs) { _inq[_q[0]] = 0; }

        for ( int i = 0, j = 2*i + 1, t; j < _qs; i = j, j = 2*i + 1) {
            if ( j + 1 < _qs && _d[_q[j + 1]] < _d[_q[j]]) j++;
            if ( _d[_q[j]] >= _d[_q[i]]) break;
            BUBL;
        }

        // relax edge (u,i) or (i,u) for all i;
        for ( int k = 0, v = _adj[u][k]; k < _deg[u]; v = _adj[u][++k]) {
            // try undoing edge v->u      
            if ( _fnet[v][u] && _d[v] > POTENTIAL(u,v) - _cost[v][u] ) 
                _d[v] = POTENTIAL(u,v) - _cost[v][_par[v] = u];

            // try using edge u->v
            if( _fnet[u][v] < _cap[u][v] && _d[v] > POTENTIAL(u,v) + _cost[u][v]) 
                _d[v] = POTENTIAL(u,v) + _cost[_par[v] = u][v];

            if ( _par[v] == u) {
                // bubble up or decrease key
                if ( _inq[v] < 0) { 
                    _inq[_q[_qs] = v] = _qs; _qs++;
                }
                for ( int i = _inq[v], j = ( i - 1 )/2, t; _d[_q[i]] < _d[_q[j]]; 
                      i = j, j = ( i - 1 )/2 ) {
                    BUBL;
                }
            }
        }
    }

    for ( int i = 0; i < n; i++) { 
        if ( _pi[i] < HALF_INFINITY ) _pi[i] += _d[i];
    }

    return ( _par[t] >= 0);
}


int MCMF_EDMONDS::run_edmonds( int n, int s, int t, int &flow_cost) 
{
    // build the adjacency list
    CLEAR_MEM_1D( _deg, 0, _NN);
    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            if ( _cap[i][j] || _cap[j][i] ) _adj[i][_deg[i]++] = j;
        }
    }

    CLEAR_MEM_2D( _fnet, 0, _NN);
    CLEAR_MEM_1D( _pi, 0, _NN);
    int flow = flow_cost = 0;

    // repeatedly, find a cheapest path from s to t
    while ( dijkstra( n, s, t) ) {
        // get the bottleneck capacity
        int bot = INT_MAX;
        int this_bot = INT_MAX;
        for ( int v = t, u = _par[v]; v != s; u = _par[v = u] ) {
            //bot <?= (_fnet[v][u] ? _fnet[v][u] : ( _cap[u][v] - _fnet[u][v]));
            this_bot = (_fnet[v][u] ? _fnet[v][u] : ( _cap[u][v] - _fnet[u][v]));
            bot = (bot < this_bot) ? bot : this_bot;
        }
        
        // update the flow network
        for ( int v = t, u = _par[v]; v != s; u = _par[v = u] ) {
            if ( _fnet[v][u] ) {
                _fnet[v][u] -= bot; 
                flow_cost -= bot * _cost[v][u]; }
            else { 
                _fnet[u][v] += bot; 
                flow_cost += bot * _cost[u][v]; 
            }
        }

        flow += bot;
    }

    return flow;
}


////////////////////////////////////////////////////////////////////////////////
//
// main
//
////////////////////////////////////////////////////////////////////////////////


int main( int argc, char *argv[])
{
    // use on a small cooked example; debugging purposes only;
    int num_vertices = 18;
    int num_edges = 27;
    MCMF_EDMONDS my_mcmf_problem( num_vertices);

    int bound_nn = 2*3 + 15 + 2 ; // 2 * _feeders_count + _ties_count + 2;
    my_mcmf_problem.set_NN_and_allocate_arrays( bound_nn);

    my_mcmf_problem.set_edge(0, 1, 1, 94);
    my_mcmf_problem.set_edge(1, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 3, 1, 66);
    my_mcmf_problem.set_edge(3, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 4, 1, 35);
    my_mcmf_problem.set_edge(4, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 5, 1, 1);
    my_mcmf_problem.set_edge(5, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 6, 1, 26);
    my_mcmf_problem.set_edge(6, 2, 1, 0);

    my_mcmf_problem.set_edge(7, 8, 1, 78);
    my_mcmf_problem.set_edge(8, 2, 1, 0);
    my_mcmf_problem.set_edge(7, 9, 1, 80);
    my_mcmf_problem.set_edge(9, 2, 1, 0);

    my_mcmf_problem.set_edge(10, 11, 1, 87);
    my_mcmf_problem.set_edge(11, 12, 1, 0);
    my_mcmf_problem.set_edge(10, 13, 1, 41);
    my_mcmf_problem.set_edge(13, 12, 1, 0);
    my_mcmf_problem.set_edge(10, 14, 1, 68);
    my_mcmf_problem.set_edge(14, 12, 1, 0);
    my_mcmf_problem.set_edge(10, 15, 1, 59);
    my_mcmf_problem.set_edge(15, 12, 1, 0);

    my_mcmf_problem.set_edge(16, 0, 1, 0);
    my_mcmf_problem.set_edge(16, 7, 1, 0);
    my_mcmf_problem.set_edge(16, 10, 1, 0);

    my_mcmf_problem.set_edge(2, 17, 1, 0);
    my_mcmf_problem.set_edge(12, 17, 1, 0);


    int flow_cost = 0;
    int max_flow = my_mcmf_problem.run_edmonds( num_vertices, 16, 17, flow_cost); // S,T;
    printf("\n flow: %d", max_flow);
    printf("\n cost: %d\n", flow_cost);
    for ( int i = 0; i < num_vertices; i++) {
        for ( int j = 0; j < num_vertices; j++) {
            int this_flow = my_mcmf_problem.get_fnet(i,j);          
            if ( this_flow)
                printf("\n %d -> %d: %d cost: %d", i, j, 
                       my_mcmf_problem.get_fnet(i,j), my_mcmf_problem.get_cost(i,j));
        }
    }

    return 0;
}
