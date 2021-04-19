#pragma once

#include <iostream>
#include <vector>
#include "debug.h"

using namespace std;

typedef double R;

struct linear_eq
{
    vector<vector<R>> A;
    vector<R> B;
    vector<int> xids;
    int n_eq;
    int n_var;
    R EPS = 1e-6;
    R _pEPS = 1e-9;
    R ninf = -1e20;

    linear_eq(int _n_eq, int _n_var)
    {
        this->n_eq = _n_eq;
        this->n_var = _n_var;
        this->A.assign(_n_eq, vector<R>(_n_var));
        this->B.assign(_n_var, 0);
        this->xids.assign(_n_var, 0);
    }

    void rswap(int i, int j)
    {
        vector<R> tp = this->A[j];
        this->A[j] = this->A[i];
        this->A[i] = tp;

        R btp = this->B[j];
        this->B[j] = this->B[i];
        this->B[i] = btp;
    }

    // R[i]<-R[i]-C*R[j]
    void r_reduce(int i, int j, R C)
    {
        for (int k = 0; k < this->n_var; k++)
        {
            this->A[i][k] -= this->A[j][k] * C;
        }
        this->B[i] -= this->B[j] * C;
    }

    void scale(int j, R C)
    {
        for (int i = 0; i < n_var; i++)
        {
            this->A[j][i] /= C;
        }
        this->B[j] /= C;
    }

    void input()
    {
        cout << "for each next " << this->n_eq << " rows input " << this->n_var + 1 << " space separtaed numbers in the format :\na[i][1] a[i][2] .. a[i][n] b[i]\n";
        for (int i = 0; i < n_eq; i++)
        {
            for (int j = 0; j < n_var; j++)
            {
                cin >> this->A[i][j];
            }
            cin >> this->B[i];
        }
    }

    void gauss_seidel_itr(int n, vector<vector<R>> &A, vector<R> &B, vector<R> &x)
    {
        for (int i = 0; i < n; i++)
        {
            x[i] = B[i];
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    x[i] -= A[i][j] * x[j];
                }
            }
            x[i] /= A[i][i];
        }
    }

    bool gauss_seidel(int k, vector<R> &x)
    {
        x.resize(this->n_eq, 0);
        for (int i = 0; i < this->n_eq; i++)
            x[i] = 0;

        double err = 1;
        vector<R> xl = x;

        while (k > 0 && err > EPS)
        {
            gauss_seidel_itr(this->n_eq, A, B, x);
            err = 0;
            for (int j = 0; j < this->n_eq; j++)
            {
                err += abs(xl[j] - x[j]);
            }
            xl = x;
            k--;
        }

        double errlim = (double)this->n_var * 10.0;
        if (err > errlim)
        {
            return false;
        }
        return true;
    }

    bool row_echilon_form(bool reduced = false)
    {
        bool full_rank = true;
        for (int i = 0; i < min(n_eq, n_var); i++)
        {
            int mxid = -1;
            R mxv = 0;
            for (int j = i; j < n_eq; j++)
            {
                if (abs(this->A[i][j]) >= mxv)
                {
                    mxid = j;
                    mxv = abs(this->A[i][j]);
                }
            }

            if (mxv < EPS)
            {
                full_rank = false;
                continue;
            }

            if (mxid != i)
                this->rswap(mxid, i);

            for (int j = i + 1; j < n_eq; j++)
            {
                this->r_reduce(j, i, this->A[j][i] / this->A[i][i]);
            }

            if (reduced)
            {
                this->scale(i, this->A[i][i]);
            }
        }

        return full_rank;
    }

    bool gauss_elemination(bool reduce = false)
    {
        bool full_rank = true;
        for (int i = 0; i < min(n_eq, n_var); i++)
        {
            int mxid = -1;
            R mxv = 0;
            for (int j = i; j < n_eq; j++)
            {
                if (abs(this->A[i][j]) >= mxv)
                {
                    mxid = j;
                    mxv = abs(this->A[i][j]);
                }
            }

            if (mxv < _pEPS)
            {
                full_rank = false;
                continue;
            }

            if (mxid != i)
                this->rswap(mxid, i);

            for (int j = 0; j < n_eq; j++)
            {
                if (j != i)
                {
                    this->r_reduce(j, i, this->A[j][i] / this->A[i][i]);
                }
            }

            if (reduce)
            {
                this->scale(i, this->A[i][i]);
            }
        }

        return full_rank;
    }

    friend ostream &operator<<(ostream &os, linear_eq &le)
    {
        for (int i = 0; i < le.n_eq; i++)
        {
            os << "[ ";
            for (int j = 0; j < le.n_var; j++)
            {
                os << le.A[i][j] << " ";
            }
            os << " | " << le.B[i];
            os << " ]\n";
        }

        return os;
    }

};