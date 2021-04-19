#include <iostream>
#include <vector>
#include "linear_eq.h"
#include "or_lab.h"
#include "debug.h"

using namespace std;
typedef double R;

void test_le()
{
    int ne, nv;
    cin >> ne >> nv;
    linear_eq le(ne, nv);
    le.input();
    cout << le;

    OR_lab orl;
    vector<vector<R>> bfs;
    orl.basic_feasible_soln(bfs, le);
    le.row_echilon_form(1);
    cout << le;
    debug(bfs);
}

void test_simplex()
{
    int bs, nbs;
    cout << "Number of basic variables: ";
    cin >> bs;
    cout << "Number of non-basic variables: ";
    cin >> nbs;
    simplex si(bs, nbs);
    si.input();
    cout << si;
    si.compute_table();
    cout << si;
    si.print_solution();
}

// void test_big_m()
// {
//     int eq,bs,nbs,art;
//     cout<<"Number of equations: ";
//     cin>>eq;
//     cout<<"Number of basic variables: ";
//     cin>>bs;
//     cout<<"Number of non-basic variables: ";
//     cin>>nbs;
//     cout<<"Number of artificial variables: ";
//     cin>>art;
//     big_m bm(eq,bs,nbs,art);
//     bm.input();
//     cout<<bm;
//     bm.compute_table();
//     cout<<bm;
//     bm.print_solution();
// }

// objective values for different problem
vector<R> zcoeff1 = {2, 5};
vector<R> zcoeff2 = {4, 3, 6};
vector<R> zcoeff3 = {12, 15, 14};
vector<R> zcoeff4 = {1, -3, 3};
vector<R> zcoeff5 = {3, 2, 2};

R objective(vector<R> &xs, int prob_id)
{
    R val = 0;
    vector<R> zcoeff;
    if(prob_id==1)
    {
        zcoeff=zcoeff1;
    }else if(prob_id==2)
    {
        zcoeff=zcoeff2;
    }else if(prob_id==3)
    {
        zcoeff=zcoeff3;
    }else if(prob_id==4)
    {
        zcoeff=zcoeff4;
    }else if(prob_id==5)
    {
        zcoeff=zcoeff5;
    }
    for (int i = 0; i < zcoeff.size(); i++)
    {
        val += zcoeff[i] * xs[i];
    }
    return val;
}

void print_solutions(vector<vector<R>> &bfs, int prob_id)
{
    int oid = 0;
    int id = 0;
    R mxv = -1.0;

    for (auto s : bfs)
    {
        bool feasible = true;
        cout << "[ ";
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                continue;
            }
            if (s[i] < 0)
                feasible = false;
            cout << "x" << i + 1 << "=" << s[i] << ", ";
        }

        if (feasible)
        {
            cout << "] Feasible and value of Z=";
            R z_val = objective(s,prob_id);
            cout << z_val;
            if (z_val > mxv)
            {
                oid = id;
                mxv = z_val;
            }
        }
        else
        {
            cout << "] Not Feasible and value of Z=";
            R z_val = objective(s,prob_id);
            cout << z_val;
        }
        cout << "\n";
        id++;
    }

    cout << "\n Among them the optimal solution is ";
    bool feasible = true;
    cout << "[ ";
    for (int i = 0; i < bfs[oid].size(); i++)
    {
        if (i == bfs[oid].size() - 1)
        {
            continue;
        }
        if (bfs[oid][i] < 0)
            feasible = false;
        cout << "x" << i + 1 << "=" << bfs[oid][i] << ", ";
    }

    cout << "] Feasible and value of Z=" << mxv << "\n";
}

void print_min_solutions(vector<vector<R>> &bfs, int prob_id)
{
    int oid = 0;
    int id = 0;
    R mxv = 1000000.0;

    for (auto s : bfs)
    {
        bool feasible = true;
        cout << "[ ";
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                continue;
            }
            if (s[i] < 0)
                feasible = false;
            cout << "x" << i + 1 << "=" << s[i] << ", ";
        }

        if (feasible)
        {
            cout << "] Feasible and value of Z=";
            R z_val = objective(s,prob_id);
            cout << z_val;
            if (z_val < mxv)
            {
                oid = id;
                mxv = z_val;
            }
        }
        else
        {
            cout << "] Not Feasible and value of Z=";
            R z_val = objective(s,prob_id);
            cout << z_val;
        }
        cout << "\n";
        id++;
    }

    cout << "\n Among them the optimal solution is ";
    bool feasible = true;
    cout << "[ ";
    for (int i = 0; i < bfs[oid].size(); i++)
    {
        if (i == bfs[oid].size() - 1)
        {
            continue;
        }
        if (bfs[oid][i] < 0)
            feasible = false;
        cout << "x" << i + 1 << "=" << bfs[oid][i] << ", ";
    }

    cout << "] Feasible and value of Z=" << mxv << "\n";
}

void lab_2()
{
    int ne, nv;
    cin >> ne >> nv;
    linear_eq le(ne, nv);
    le.input();
    // cout << le;

    OR_lab orl;
    vector<vector<R>> bfs;
    orl.basic_feasible_soln(bfs, le);
    // le.row_echilon_form(1);
    // cout<<le;
    // debug(bfs);
    cout<<"BFS's are:\n";
    // for maximization problems
    // in 2nd input just give the problem number
    // print_solutions(bfs,5);

    // for minimization problems
    print_min_solutions(bfs,4);
    
}

/*

*/

int main()
{
    // test_simplex();
    lab_2();

    return 0;
}

/*
Number of basic variables: 4
Number of non-basic variables: 2
for each next 5 rows input 8 space separtaed numbers in the format :
z x[i][1] ... x[i][n_nbs] s[i][1] ... s[i][n_slc] Sol[i]
1 -5 -4 0 0 0 0 0
0 6 4 1 0 0 0 24
0 1 2 0 1 0 0 6
0 -1 1 0 0 1 0 1
0 0 1 0 0 0 1 2
0 [ 1 -5 -4 0 0 0 0  | 0 ]
3 [ 0 6 4 1 0 0 0  | 24 ]
4 [ 0 1 2 0 1 0 0  | 6 ]
5 [ 0 -1 1 0 0 1 0  | 1 ]
6 [ 0 0 1 0 0 0 1  | 2 ]
0 [ 1 0 0 0.75 0.5 0 0  | 21 ]
1 [ 0 1 0 0.25 -0.5 0 0  | 3 ]
2 [ 0 0 1 -0.125 0.75 0 0  | 1.5 ]
5 [ 0 0 0 0.375 -1.25 1 0  | 2.5 ]
6 [ 0 0 0 0.125 -0.75 0 1  | 0.5 ]
z=21
x1=3
x2=1.5
s3=2.5
s4=0.5

*/