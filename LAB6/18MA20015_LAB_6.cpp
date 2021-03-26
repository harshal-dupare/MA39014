/*
18MA20015 | Harshal Dupare | Lab 6 | 17-2-2021 |
*/

#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;
typedef double R;
typedef long long ll;
typedef vector<vector<R>> vvr;
typedef vector<R> vr;
ll itrcount = 0;
ll wd = 14;
R Z_sum;

struct linear_eq
{
    vector<vector<R>> A;
    vector<R> B;
    vector<ll> xids;
    ll n_eq;
    ll n_var;
    R EPS = 1e-6;
    R _pEPS = 1e-9;
    R ninf = -1e20;

    linear_eq(ll _n_eq, ll _n_var)
    {
        this->n_eq = _n_eq;
        this->n_var = _n_var;
        this->A.assign(_n_eq, vector<R>(_n_var, (R)0));
        this->B.assign(_n_var, 0);
        this->xids.assign(_n_var, 0);
    }

    void rswap(ll i, ll j)
    {
        vector<R> tp = this->A[j];
        this->A[j] = this->A[i];
        this->A[i] = tp;

        R btp = this->B[j];
        this->B[j] = this->B[i];
        this->B[i] = btp;
    }

    // R[i]<-R[i]-C*R[j]
    void r_reduce(ll i, ll j, R C)
    {
        for (ll k = 0; k < this->n_var; k++)
        {
            this->A[i][k] -= this->A[j][k] * C;
        }
        this->B[i] -= this->B[j] * C;
    }

    void scale(ll j, R C)
    {
        for (ll i = 0; i < n_var; i++)
        {
            this->A[j][i] /= C;
        }
        this->B[j] /= C;
    }

    void input()
    {
        cout << "for each next " << this->n_eq << " rows input " << this->n_var + 1 << " space separtaed numbers in the format :\na[i][1] a[i][2] .. a[i][n] b[i]\n";
        for (ll i = 0; i < n_eq; i++)
        {
            for (ll j = 0; j < n_var; j++)
            {
                cin >> this->A[i][j];
            }
            cin >> this->B[i];
        }
    }

    void gauss_seidel_itr(ll n, vector<vector<R>> &A, vector<R> &B, vector<R> &x)
    {
        for (ll i = 0; i < n; i++)
        {
            x[i] = B[i];
            for (ll j = 0; j < n; j++)
            {
                if (i != j)
                {
                    x[i] -= A[i][j] * x[j];
                }
            }
            x[i] /= A[i][i];
        }
    }

    bool gauss_seidel(ll k, vector<R> &x)
    {
        x.resize(this->n_eq, 0);
        for (ll i = 0; i < this->n_eq; i++)
            x[i] = 0;

        R err = 1;
        vector<R> xl = x;

        while (k > 0 && err > EPS)
        {
            gauss_seidel_itr(this->n_eq, A, B, x);
            err = 0;
            for (ll j = 0; j < this->n_eq; j++)
            {
                err += abs(xl[j] - x[j]);
            }
            xl = x;
            k--;
        }

        R errlim = (R)this->n_var * 10.0;
        if (err > errlim)
        {
            return false;
        }
        return true;
    }

    bool row_echilon_form(bool reduced = false)
    {
        bool full_rank = true;
        for (ll i = 0; i < min(n_eq, n_var); i++)
        {
            ll mxid = -1;
            R mxv = 0;
            for (ll j = i; j < n_eq; j++)
            {
                if (abs(this->A[j][i]) >= mxv)
                {
                    mxid = j;
                    mxv = abs(this->A[j][i]);
                }
            }

            if (mxv < EPS)
            {
                full_rank = false;
                continue;
            }

            if (mxid != i)
                this->rswap(mxid, i);

            for (ll j = i + 1; j < n_eq; j++)
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
        for (ll i = 0; i < min(n_eq, n_var); i++)
        {
            ll mxid = -1;
            R mxv = 0;
            for (ll j = i; j < n_eq; j++)
            {
                if (abs(this->A[j][i]) >= mxv)
                {
                    mxid = j;
                    mxv = abs(this->A[j][i]);
                }
            }

            if (mxv < EPS)
            {
                full_rank = false;
                continue;
            }
            // debug(mxv);
            // debug(mxid);
            // debug(i);
            if (mxid != i)
                this->rswap(mxid, i);

            for (ll j = 0; j < n_eq; j++)
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
        for (ll i = 0; i < le.n_eq; i++)
        {
            os << "[ ";
            for (ll j = 0; j < le.n_var; j++)
            {
                os << le.A[i][j] << " ";
            }
            os << " | " << le.B[i];
            os << " ]\n";
        }

        return os;
    }
};

struct DS
{
    vvr A;
    vr B, C;
    int n, m;
    int K_col, K_row;
    R MX_mz = 0;
    R MN_sol = 1e6;
    vector<int> Basic_vars;
    vr Z;
    vr C_mn_Z;
    vr sol;
    vr ratio;
    vr K_row_val;
    vr K_col_val;
    vr coeff_Basic_vars;
    R Z_sum;
    int k;
    vvr SIMPLEX_TABLE;
    int maxitr = 10;

    DS(int _n, int _m)
    {
        n = _n;
        m = _m;
        A = vvr(m, vr(n));
        B = vr(m);
        C = vr(n);
        Basic_vars = vector<int>(m);
        Z = vr(n + m);
        C_mn_Z = vr(n + m);
        sol = vr(m);
        ratio = vr(m);
        K_row_val = vr(n + m);
        K_col_val = vr(m);
        coeff_Basic_vars = vr(n + m);
        k = n + m;
        SIMPLEX_TABLE = vvr(m, vr(k));
    }

    void input()
    {
        cout << "Enter the coefficients of variables in objective\nif cx[i] is coefficient of x_i then give input in format\ncx[1] cx[2] ... cx[n]\n";
        for (int i = 0; i < n; i++)
        {
            cin >> this->C[i];
        }

        cout << "For each condition give input of the coefficient of the variables in < condition, in the format below\n "
             << "coeff[i][0] .... coeff[i][number_variable] val[i]\n ";

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cin >> this->A[i][j];
            }
            cin >> this->B[i];
        }
    }

    void init()
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                SIMPLEX_TABLE[i][j] = A[i][j];
            }
            for (int j = n; j < (n + m); j++)
            {
                SIMPLEX_TABLE[i][j] = 0.0;
            }
            SIMPLEX_TABLE[i][n + i] = 1.0;
        }
        vr coeff_Basic_vars(n + m);
        for (int i = 0; i < n; i++)
        {
            coeff_Basic_vars[i] = C[i];
        }
        for (int j = n; j < (n + m); j++)
        {
            coeff_Basic_vars[j] = 0.0;
        }

        for (int i = 0; i < m; i++)
        {
            Basic_vars[i] = i + n;
            sol[i] = B[i];
        }
    }
};

void dual_Simplex(vvr &A, vr &B, vr &C, int n, int m)
{
    int i, j, k = n + m;
    vvr SIMPLEX_TABLE(m, vr(k));
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            SIMPLEX_TABLE[i][j] = A[i][j];
        }
        for (j = n; j < (n + m); j++)
        {
            SIMPLEX_TABLE[i][j] = 0.0;
        }
        SIMPLEX_TABLE[i][n + i] = 1.0;
    }
    vr coeff_Basic_vars(n + m);
    for (i = 0; i < n; i++)
    {
        coeff_Basic_vars[i] = C[i];
    }
    for (j = n; j < (n + m); j++)
    {
        coeff_Basic_vars[j] = 0.0;
    }
    int K_col, K_row;
    R MX_mz = 0.0;
    R MN_sol = 999999.0;
    vector<int> Basic_vars(m);
    vr Z(n + m);
    vr C_mn_Z(n + m);
    vr sol(m);
    vr ratio(m);
    vr K_row_val(n + m);
    vr K_col_val(m);
    for (i = 0; i < m; i++)
    {
        Basic_vars[i] = i + n;
        sol[i] = B[i];
    }
    int check = 1;
    int iter = 0;
    R sol_K;
    R Z_sol;

    while (check)
    {
        if (iter >= 11)
            break;
        if (iter != 0)
        {
            Basic_vars[K_row] = K_col;
            for (i = 0; i < m; i++)
            {
                if (i == K_row)
                {
                    sol[i] = sol[i] / K_row_val[K_col];
                }
                else
                {
                    sol[i] = sol[i] - (K_col_val[i] * sol_K) / K_row_val[K_col];
                }

                for (j = 0; j < (m + n); j++)
                {
                    if (i == K_row)
                    {
                        SIMPLEX_TABLE[i][j] = SIMPLEX_TABLE[i][j] / K_row_val[K_col];
                    }
                    else
                    {
                        SIMPLEX_TABLE[i][j] = SIMPLEX_TABLE[i][j] - (K_col_val[i] * K_row_val[j]) / K_row_val[K_col];
                    }
                }
            }
        }
        check = 0;

        cout << "\n\t Iteration : " << iter << endl;
        cout << "\n\t CB_i \t C_j ";
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", coeff_Basic_vars[i]);
        }
        cout << "\n \t \t BV. ";
        for (i = 0; i < (m + n); i++)
        {
            cout << "\t     x_" << i + 1;
        }
        cout << "\t Solution\n";
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            cout << "_";
        }
        cout << endl;
        for (i = 0; i < m; i++)
        {
            printf("\t %0.2lf    x_%d ", coeff_Basic_vars[Basic_vars[i]], Basic_vars[i] + 1);
            for (j = 0; j < (m + n); j++)
            {
                cout << "\t" << SIMPLEX_TABLE[i][j];
            }
            cout << "\t" << sol[i] << endl;
        }
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            cout << "_";
        }
        for (i = 0; i < (m + n); i++)
        {
            Z_sum = 0;
            for (j = 0; j < m; j++)
            {
                Z_sum += coeff_Basic_vars[Basic_vars[j]] * SIMPLEX_TABLE[j][i];
            }
            Z[i] = Z_sum;
            C_mn_Z[i] = coeff_Basic_vars[i] - Z[i];
        }
        printf("\n\t Z_j \t ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", Z[i]);
        }
        Z_sum = 0;
        for (i = 0; i < m; i++)
        {
            Z_sum += coeff_Basic_vars[Basic_vars[i]] * sol[i];
        }
        Z_sol = Z_sum;
        printf("\t %lf", Z_sol);
        printf("\n \t C_j - Z-j ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", C_mn_Z[i]);
        }
        printf("\n");
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            cout << "_";
        }
        printf("\n");

        MN_sol = 0;
        for (i = 0; i < m; i++)
        {
            if (sol[i] < MN_sol)
            {
                MN_sol = sol[i];
                K_row = i;
            }
            if (sol[i] < 0)
            {
                check = 1;
            }
        }
        if (check == 0)
        {
            break;
        }
        printf("The leaving variable is : x_%d \n", Basic_vars[K_row] + 1);
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("_");
        }
        printf("\n");
        printf("Variables ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t   x_%d       ", i + 1);
        }
        printf("\n");
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("_");
        }
        printf("\n -(C_j - Z-j)");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", -1 * C_mn_Z[i]);
        }
        printf("\n x_%d \t", K_row + 1);
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", SIMPLEX_TABLE[K_row][i]);
        }
        R minratio = 100000;
        R ratio2;
        printf("\n Ratio \t");
        for (i = 0; i < (m + n); i++)
        {
            if (SIMPLEX_TABLE[K_row][i] < 0)
            {
                ratio2 = (-1 * C_mn_Z[i]) / SIMPLEX_TABLE[K_row][i];
                printf("\t %lf", ratio2);
                if (ratio2 < minratio)
                {
                    minratio = ratio2;
                    K_col = i;
                }
            }
            else
            {
                printf("\t  --      ");
            }
        }
        printf("\n");
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("_");
        }
        cout << "\nMinimum ratio is " << minratio << "So the entering variable is x_" << K_col + 1 << endl;
        for (i = 0; i < (m + n); i++)
        {
            K_row_val[i] = SIMPLEX_TABLE[K_row][i];
        }
        sol_K = sol[K_row];
        for (i = 0; i < m; i++)
        {
            K_col_val[i] = SIMPLEX_TABLE[i][K_col];
        }
        iter++;
    }

    if (iter < 10)
    {
        cout << "\nThe Optimal solutions is : ";
        for (i = 0; i < m; i++)
        {
            cout << " x_" << Basic_vars[i] + 1 << "=" << sol[i] << ",";
        }
        cout << "\nAnd the optimal value of Z is : " << Z_sol << endl;
    }
}

void dl()
{
    int n, m, i, j, k;
    cout << "Enter number of variable : ";
    cin >> n;
    cout << "Enter number of equation : ";
    cin >> m;
    vvr A(m, vr(n));
    vr B(m);
    vr C(n);
    cout << "Enter the coefficients of variables in objective\nif cx[i] is coefficient of x_i then give input in format\ncx[1] cx[2] ... cx[n]\n";
    for (int i = 0; i < n; i++)
    {
        cin >> C[i];
    }

    cout << "For each condition give input of the coefficient of the variables in < condition, in the format below\n "
         << "coeff[i][0] .... coeff[i][number_variable] val[i]\n";

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> A[i][j];
        }
        cin >> B[i];
    }

    dual_Simplex(A, B, C, n, m);
}

int main()
{
    dl();
    return 0;
}