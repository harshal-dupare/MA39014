#pragma once

#include <iostream>
#include <vector>
#include "debug.h"
#include "linear_eq.h"

using namespace std;

typedef double R;

/*
Baseline model, linear programming problem, convex sets, convex functions and
their properties, basic feasible solution, optimal solution, related theorems.

Graphical method for solving two and three variable problems, 

simplex method,
Big M method, degenerate LP problem, product form of inverse of a matrix,
revised simplex method, duality theorems, complementary slackness principle,
primal-dual simplex algorithm, sensitivity analysis, parametric programming,
linear integer programming problem, Gomory cutting plane method, branch and
bound algorithm, 0-1 implicit enumeration, transportation problem, assignment
problem with their solution methodologies. 

Theory of games, two-person zero-sum
games with and without saddle-points, pure and mixed strategies, graphical
method of solution of a 2ï´n game, solution of an mï´n game by simplex method.

*/

struct simplex
{
    vector<vector<R>> A;
    vector<R> sol;
    vector<R> z;
    vector<int> xids;
    int n_bsc;
    int n_nbsc;
    int n_tot;
    R EPS = 1e-6;
    R _pEPS = 1e-9;
    R ninf = -1e20;
    R inf = 1e20;

    simplex(int _n_bsc, int _n_nbsc)
    {
        this->n_bsc = _n_bsc;
        this->n_nbsc = _n_nbsc;
        this->n_tot = _n_nbsc + _n_bsc;
        this->A.assign(_n_bsc + 1, vector<R>(_n_nbsc + _n_bsc + 1));
        this->sol.assign(_n_bsc + 1, 0);
        this->xids.assign(_n_bsc + 1, 0);
        for (int i = 1; i < _n_bsc + 1; i++)
        {
            this->xids[i] = i + _n_nbsc;
        }
    }

    void rswap(int i, int j)
    {
        vector<R> tp = this->A[j];
        this->A[j] = this->A[i];
        this->A[i] = tp;

        R btp = this->sol[j];
        this->sol[j] = this->sol[i];
        this->sol[i] = btp;
    }

    // R[i]<-R[i]-C*R[j]
    void r_reduce(int i, int j, R C)
    {
        for (int k = 0; k < this->n_tot; k++)
        {
            this->A[i][k] -= this->A[j][k] * C;
        }
        this->sol[i] -= this->sol[j] * C;
    }

    void scale(int j, R C)
    {
        for (int i = 0; i < n_tot; i++)
        {
            this->A[j][i] /= C;
        }
        this->sol[j] /= C;
    }

    // optamilaity condition
    int new_entering_variable(int min_or_max = 1)
    {
        R mxnv = 0.0;
        int mxid = -1;
        for (int i = 0; i < n_tot + 1; i++)
        {
            if ((this->A[0][i] * ((R)min_or_max)) < 0.0 && abs(this->A[0][i]) > mxnv)
            {
                mxnv = abs(this->A[0][i]);
                mxid = i;
            }
        }

        return mxid;
    }

    // feasibility condiion
    int new_leaving_variable(int k)
    {
        int mnid = -1;
        R mnv = inf;
        for (int i = 1; i < n_bsc + 1; i++)
        {
            if (this->A[i][k] > 0 && ((this->sol[i] / this->A[i][k]) < mnv))
            {
                mnv = (this->sol[i] / this->A[i][k]);
                mnid = i;
            }
        }

        return mnid;
    }

    // 1 for max problem -1 for min problem
    void compute_table(int min_or_max = 1)
    {
        while (true)
        {
            int new_etr = this->new_entering_variable(min_or_max);
            if (new_etr == -1)
            {
                break;
            }

            int new_lev = this->new_leaving_variable(new_etr);
            if (new_lev == -1)
            {
                break;
            }

            this->scale(new_lev, this->A[new_lev][new_etr]);

            for (int i = 0; i < n_bsc + 1; i++)
            {
                if (i != new_lev)
                {
                    this->r_reduce(i, new_lev, this->A[i][new_etr]);
                }
            }

            this->xids[new_lev] = new_etr;
        }
    }

    void input()
    {

        cout << "for each next " << this->n_bsc + 1 << " rows input " << this->n_tot + 2 << " space separtaed numbers in the format :\nz x[i][1] ... x[i][n_nbs] s[i][1] ... s[i][n_slc] Sol[i]\n";
        for (int i = 0; i < n_bsc + 1; i++)
        {
            for (int j = 0; j < n_tot + 1; j++)
            {
                cin >> this->A[i][j];
            }
            cin >> this->sol[i];
        }
    }

    friend ostream &operator<<(ostream &os, simplex &si)
    {
        for (int i = 0; i < si.n_bsc + 1; i++)
        {
            os << si.xids[i] << " [ ";
            for (int j = 0; j < si.n_tot + 1; j++)
            {
                os << si.A[i][j] << " ";
            }
            os << " | " << si.sol[i];
            os << " ]\n";
        }

        return os;
    }

    void print_solution()
    {
        cout << "z=" << sol[0] << "\n";
        for (int i = 1; i < n_bsc + 1; i++)
        {
            if (xids[i] < n_nbsc + 1)
                cout << "x" << xids[i] << "=" << sol[i] << "\n";
            else
                cout << "s" << xids[i] - n_nbsc << "=" << sol[i] << "\n";
        }
    }
};

struct OR_lab
{
    // goes through all subset of size m and applies GS/GE and saves soln
    void compute_bfs(linear_eq &le, vector<int> &subset, int lst, int k, vector<vector<R>> &bfs)
    {
        if (subset.size() < le.n_eq)
        {
            for (int i = lst; i < le.n_var; i++)
            {
                subset.push_back(i);
                compute_bfs(le, subset, i + 1, k, bfs);
                subset.pop_back();
            }
            return;
        }
        linear_eq le_new(le.n_eq, le.n_eq);

        for (int i = 0; i < le.n_eq; i++)
        {
            for (int j = 0; j < subset.size(); j++)
            {
                le_new.A[i][j] = le.A[i][subset[j]];
            }
            le_new.B[i] = le.B[i];
        }
        vector<R> x;

        // // if using gauss sidel
        // bool converge = le_new.gauss_seidel(k,x);

        // if using gauss elimination
        bool converge = le_new.gauss_elemination(true);
        x.assign(le_new.n_eq, 0);
        for (int i = 0; i < le.n_eq; i++)
        {
            if (le_new.A[i][i] < le._pEPS)
            {
                if (le_new.B[i] >= le._pEPS)
                {
                    x[i] = -1;
                }
            }
            else
            {
                x[i] = le_new.B[i];
            }
        }

        vector<R> xx(le.n_var, 0);
        for (int i = 0; i < subset.size(); i++)
        {
            xx[subset[i]] = x[i];
        }

        // save soln
        if (converge)
            xx.push_back(1);
        else
            xx.push_back(0);

        bfs.push_back(xx);
    }

    void basic_feasible_soln(vector<vector<R>> &bfs, linear_eq &le, int k = 100)
    {
        vector<int> subset;
        compute_bfs(le, subset, 0, 1, bfs);
    }

    void find_optimal_fs(vector<vector<R>> &bfs, R (*f)(vector<R> &))
    {
    }

    // big m

    // dual simplex

    // revise simplex

    // 2 phase simplex

    // branch and bound

    // cutting plain

    // transport
};

// // need to complete
// struct big_m
// {
//     vector<vector<R>> A;
//     vector<R> sol;
//     vector<R> z;
//     vector<int> xids;
//     R M = 300.0;
//     int n_eq;
//     int n_gt;
//     int n_lt;
//     int n_art;
//     int n_slc;
//     int n_tot;
//     R EPS = 1e-6;
//     R _pEPS = 1e-9;
//     R ninf = -1e20;
//     R inf = 1e20;

//     big_m(int _n_lt, int _n_eq, int _n_gt)
//     {
//         this->n_eq = _n_eq;
//         this->n_art = _n_gt+_n_eq;
//         this->n_tot = _n_nbsc + _n_bsc + _n_art;
//         this->A.assign(_n_eq + 1, vector<R>(this->n_tot+1));
//         this->sol.assign(_n_eq + 1, 0);
//         this->xids.assign(_n_eq + 1, 0);
//         for (int i = 1; i < _n_eq + 1; i++)
//         {
//             this->xids[i] = i + _n_nbsc;
//         }
//     }

//     void rswap(int i, int j)
//     {
//         vector<R> tp = this->A[j];
//         this->A[j] = this->A[i];
//         this->A[i] = tp;

//         R btp = this->sol[j];
//         this->sol[j] = this->sol[i];
//         this->sol[i] = btp;
//     }

//     // R[i]<-R[i]-C*R[j]
//     void r_reduce(int i, int j, R C)
//     {
//         for (int k = 0; k < this->n_tot; k++)
//         {
//             this->A[i][k] -= this->A[j][k] * C;
//         }
//         this->sol[i] -= this->sol[j] * C;
//     }

//     void scale(int j, R C)
//     {
//         for (int i = 0; i < n_tot; i++)
//         {
//             this->A[j][i] /= C;
//         }
//         this->sol[j] /= C;
//     }

//     // optamilaity condition
//     int new_entering_variable(int min_or_max = 1)
//     {
//         R mxnv = 0.0;
//         int mxid = -1;
//         for (int i = 0; i < n_bsc + 1; i++)
//         {
//             if ((this->A[0][i] * ((R)min_or_max)) < 0.0 && abs(this->A[0][i]) > mxnv)
//             {
//                 mxnv = abs(this->A[0][i]);
//                 mxid = i;
//             }
//         }

//         return mxid;
//     }

//     // feasibility condiion
//     int new_leaving_variable(int k)
//     {
//         int mnid = -1;
//         R mnv = inf;
//         for (int i = 1; i < n_bsc + 1; i++)
//         {
//             if (this->A[i][k] > 0 && ((this->sol[i] / this->A[i][k]) < mnv))
//             {
//                 mnv = (this->sol[i] / this->A[i][k]);
//                 mnid = i;
//             }
//         }

//         return mnid;
//     }

//     // 1 for max problem -1 for min problem
//     void compute_table(int min_or_max = 1)
//     {
//         // preprocess the z row
//         vector<R> art_coeff(this->n_art);
//         for(int i=n_bsc+n_nbsc+1;i<n_tot+1;i++)
//         {
//             art_coeff[i-n_bsc-n_nbsc-1] = A[0][i];
//         }

//         for(int i=n_bsc+n_nbsc+1;i<n_tot+1;i++)
//         {
//             this->r_reduce(0,)
//         }


//         while (true)
//         {
//             int new_etr = this->new_entering_variable(min_or_max);
//             if (new_etr == -1)
//             {
//                 break;
//             }

//             int new_lev = this->new_leaving_variable(new_etr);
//             if (new_lev == -1)
//             {
//                 break;
//             }

//             this->scale(new_lev, this->A[new_lev][new_etr]);

//             for (int i = 0; i < n_bsc + 1; i++)
//             {
//                 if (i != new_lev)
//                 {
//                     this->r_reduce(i, new_lev, this->A[i][new_etr]);
//                 }
//             }

//             this->xids[new_lev] = new_etr;
//         }
//     }

//     void input()
//     {

//         cout << "for each next " << this->n_bsc + 1 << " rows input " << this->n_tot + 2 << " space separtaed numbers in the format :\nz x[i][1]..x[i][n_nbs] s[i][1]..s[i][n_slc] art[i][1]..art[i][n_art] Sol[i]\n";
//         for (int i = 0; i < n_eq + 1; i++)
//         {
//             for (int j = 0; j < n_tot + 1; j++)
//             {
//                 cin >> this->A[i][j];
//             }
//             cin >> this->sol[i];
//         }
//     }

//     friend ostream &operator<<(ostream &os, big_m &si)
//     {
//         for (int i = 0; i < si.n_eq + 1; i++)
//         {
//             os << si.xids[i] << " [ ";
//             for (int j = 0; j < si.n_tot + 1; j++)
//             {
//                 os << si.A[i][j] << " ";
//             }
//             os << " | " << si.sol[i];
//             os << " ]\n";
//         }

//         return os;
//     }

//     void print_solution()
//     {
//         cout << "z=" << sol[0] << "\n";
//         for (int i = 1; i < n_bsc + 1; i++)
//         {
//             if (xids[i] < n_nbsc + 1)
//                 cout << "x" << xids[i] << "=" << sol[i] << "\n";
//             else
//                 cout << "s" << xids[i] - n_nbsc << "=" << sol[i] << "\n";
//         }
//     }
// };
