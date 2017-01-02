// Cristinel Ababei, December 2009, Fargo ND
// E-mail: cristinel.ababei@ndsu.edu
//
// This is an adapted (i.e., ported to C++) version of the Edmonds-Karp 
// algo from: http://shygypsy.com/tools/mcmf4.cpp
// The original authors of that implementation: Frank Chu and Igor Naverniouk.

#ifndef _MCMF4_H_
#define _MCMF4_H_

#include <stdlib.h>
#include <stdio.h>
#include <climits>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//
// MCMF_EDMONDS
//
////////////////////////////////////////////////////////////////////////////////


#define CLEAR_MEM_1D(p,v,n) { for(int i=0; i<n; i++) p[i]=v; }
#define CLEAR_MEM_2D(p,v,n) { for(int i=0; i<n; i++) for(int j=0; j<n; j++) p[i][j]=v; }
#define HALF_INFINITY (INT_MAX/2)


class MCMF_EDMONDS
{
 private:

#define BUBL {                                                          \
        t = _q[i]; _q[i] = _q[j]; _q[j] = t;                            \
        t = _inq[_q[i]]; _inq[_q[i]] = _inq[_q[j]]; _inq[_q[j]] = t; }

 // Dijkstra's using non-negative edge weights (cost + potential)
#define POTENTIAL(u,v) ( _d[u] + _pi[u] - _pi[v])

 private:
     int _num_vertices;
 int _qs;
 int _NN;
 // adjacency matrix (fill this up)
 int **_cap;
 // cost per unit of flow matrix (fill this up)
 int **_cost;
 // flow network and adjacency list
 int **_fnet, **_adj;
 int *_deg;
 // labelling function
 int *_pi;
 // Dijkstra's predecessor, depth and priority queue
 int *_par, *_d, *_q, *_inq;

 public:
 MCMF_EDMONDS() { _NN = 0; }
 MCMF_EDMONDS(int num_vertices) {
     _NN = 0;
     _num_vertices = num_vertices;
 }
 ~MCMF_EDMONDS() {}

 void clear_all() {
     if ( _NN > 0) {
         CLEAR_MEM_2D( _cap, 0, _NN);
         CLEAR_MEM_2D( _cost, 0, _NN);
         CLEAR_MEM_2D( _fnet, 0, _NN);
         CLEAR_MEM_2D( _adj, 0, _NN);
         CLEAR_MEM_1D( _deg, 0, _NN);
         CLEAR_MEM_1D( _pi, 0, _NN);
         CLEAR_MEM_1D( _par, 0, _NN);
         CLEAR_MEM_1D( _d, 0, _NN);
         CLEAR_MEM_1D( _q, 0, _NN);
         CLEAR_MEM_1D( _inq, 0, _NN);
     }
 }
 void set_NN_and_allocate_arrays( int nn) {
     _NN = nn;
     _cap = (int**) malloc(_NN * sizeof(int *));
     _cost = (int**) malloc(_NN * sizeof(int *));
     _fnet = (int**) malloc(_NN * sizeof(int *));
     _adj = (int**) malloc(_NN * sizeof(int *));
     if ( _cap == NULL) { printf("\nError: malloc 1\n"); exit(1); }
     if ( _cost == NULL) { printf("\nError: malloc 2\n"); exit(1); }
     if ( _fnet == NULL) { printf("\nError: malloc 3\n"); exit(1); }
     if ( _adj == NULL) { printf("\nError: malloc 4\n"); exit(1); }
     for ( int i = 0; i < _NN; i++) {
         _cap[i] = (int *) malloc(_NN * sizeof(int));
         _cost[i] = (int *) malloc(_NN * sizeof(int));
         _fnet[i] = (int *) malloc(_NN * sizeof(int));
         _adj[i] = (int *) malloc(_NN * sizeof(int));
         if ( _cap[i] == NULL) { printf("\nError: malloc 5\n"); exit(1); }
         if ( _cost[i] == NULL) { printf("\nError: malloc 6\n"); exit(1); }
         if ( _fnet[i] == NULL) { printf("\nError: malloc 7\n"); exit(1); }
         if ( _adj[i] == NULL) { printf("\nError: malloc 8\n"); exit(1); }
     }
     _deg = (int*) malloc(_NN * sizeof(int));
     _pi = (int*) malloc(_NN * sizeof(int));
     _par = (int*) malloc(_NN * sizeof(int));
     _d = (int*) malloc(_NN * sizeof(int));
     _q = (int*) malloc(_NN * sizeof(int));
     _inq = (int*) malloc(_NN * sizeof(int));
     if ( _deg == NULL) { printf("\nError: malloc 9\n"); exit(1); }
     if ( _pi == NULL) { printf("\nError: malloc 10\n"); exit(1); }
     if ( _par == NULL) { printf("\nError: malloc 11\n"); exit(1); }
     if ( _d == NULL) { printf("\nError: malloc 12\n"); exit(1); }
     if ( _q == NULL) { printf("\nError: malloc 13\n"); exit(1); }
     if ( _inq == NULL) { printf("\nError: malloc 14\n"); exit(1); }
     //printf ("\nAllocated arrays of MCMF solver.");

     // reset all allocated mem;
     clear_all();
 }
 // Note: the CS2 engine has deallocate_arrays() counterpart;
 void free_arrays() {
     for ( int i = 0; i < _NN; i++) {
         free( _cap[i]); free( _cost[i]); free( _fnet[i]); free( _adj[i]);
     }
     free( _cap); free( _cost); free( _fnet); free( _adj);
     free( _deg); free( _pi); free( _par); free( _d); free( _q); free( _inq);
 }
 void set_num_vertices( int num_vertices) { _num_vertices = num_vertices; }
 void set_edge(int i, int j, int edge_capacity, int edge_cost) {
     assert( i >= 0 && i < _num_vertices);
     assert( j >= 0 && j < _num_vertices);
     _cap[i][j] = edge_capacity; 
     _cost[i][j] = edge_cost; 
 }
 bool dijkstra( int n, int s, int t);
 int run_edmonds( int n, int s, int t, int &flow_cost);
 inline int potential(int u, int v) {
     return ( _d[u] + _pi[u] - _pi[v]);
 }
 void buble(int i, int j, int &t) {
     t = _q[i]; _q[i] = _q[j]; _q[j] = t;
     t = _inq[_q[i]]; _inq[_q[i]] = _inq[_q[j]]; _inq[_q[j]] = t;
 }
 int get_fnet(int u, int v) const { return _fnet[u][v]; }
 int get_cost(int u, int v) const { return _cost[u][v]; }
};

#endif
