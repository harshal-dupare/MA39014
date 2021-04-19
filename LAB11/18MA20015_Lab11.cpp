/*
18MA20015 | Harshal Dupare | Lab 11 | 31-3-2021 |
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
typedef double R;
typedef vector<vector<R>> vvr;
typedef vector<R> vr;
typedef long long ll;
ll itrcount = 0;
ll prtitr = 0;
ll wd = 10;
R _EPSILON = 1e-4;
R infinity_R = 1e15;
bool printtbl = false;
vector<R> zcoeff;

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

        double err = 1;
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

struct simplex
{
    vector<vector<R>> A;
    vector<R> sol;
    vector<R> z;
    vector<ll> xids;
    ll n_bsc;
    ll n_nbsc;
    ll n_tot;
    R EPS = 1e-6;
    R _pEPS = 1e-9;
    R ninf = -1e20;
    R inf = 1e20;

    simplex()
    {
    }

    simplex(ll _n_bsc, ll _n_nbsc)
    {
        this->n_bsc = _n_bsc;
        this->n_nbsc = _n_nbsc;
        this->n_tot = _n_nbsc + _n_bsc;
        this->A.assign(_n_bsc + 1, vector<R>(_n_nbsc + _n_bsc + 1));
        this->sol.assign(_n_bsc + 1, 0);
        this->xids.assign(_n_bsc + 1, 0);
        for (ll i = 1; i < _n_bsc + 1; i++)
        {
            this->xids[i] = i + _n_nbsc;
        }
    }

    void init(ll _n_bsc, ll _n_nbsc)
    {
        this->n_bsc = _n_bsc;
        this->n_nbsc = _n_nbsc;
        this->n_tot = _n_nbsc + _n_bsc;
        this->A.assign(_n_bsc + 1, vector<R>(_n_nbsc + _n_bsc + 1));
        this->sol.assign(_n_bsc + 1, 0);
        this->xids.assign(_n_bsc + 1, 0);
        for (ll i = 1; i < _n_bsc + 1; i++)
        {
            this->xids[i] = i + _n_nbsc;
        }
    }

    void rswap(ll i, ll j)
    {
        vector<R> tp = this->A[j];
        this->A[j] = this->A[i];
        this->A[i] = tp;

        R btp = this->sol[j];
        this->sol[j] = this->sol[i];
        this->sol[i] = btp;
    }

    // R[i]<-R[i]-C*R[j]
    void r_reduce(ll i, ll j, R C)
    {
        for (ll k = 0; k < this->n_tot + 1; k++)
        {
            this->A[i][k] -= this->A[j][k] * C;
        }
        this->sol[i] -= this->sol[j] * C;
    }

    void scale(ll j, R C)
    {
        for (ll i = 0; i < n_tot + 1; i++)
        {
            this->A[j][i] /= C;
        }
        this->sol[j] /= C;
    }

    // optamilaity condition
    ll new_entering_variable(ll min_or_max = 1)
    {
        R mxnv = 0.0;
        ll mxid = -1;
        for (ll i = 1; i < n_tot + 1; i++)
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
    ll new_leaving_variable(ll k)
    {
        ll mnid = -1;
        R mnv = inf;
        for (ll i = 1; i < n_bsc + 1; i++)
        {
            if (this->A[i][k] > 0 && ((this->sol[i] / this->A[i][k]) < mnv))
            {
                mnv = (this->sol[i] / this->A[i][k]);
                mnid = i;
            }
        }

        return mnid;
    }

    void printvar_category()
    {
        vector<ll> bsc(this->n_tot + 1, 0);
        vector<ll> bvr, nbvr;
        for (auto j : this->xids)
        {
            bsc[j] = 1;
        }

        for (ll i = 1; i < this->n_tot + 1; i++)
        {
            if (bsc[i] == 1)
            {
                bvr.push_back(i);
            }
            else
            {
                nbvr.push_back(i);
            }
        }

        string ss;
        cout << "Non Basic Variables Are :";
        for (auto j : bvr)
        {
            if (j > this->n_nbsc)
            {
                ss = "s_" + to_string(j - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(j);
            }
            cout << ss << ", ";
        }
        cout << "\n";

        cout << "Basic Variables Are :";
        for (auto j : nbvr)
        {
            if (j > this->n_nbsc)
            {
                ss = "s_" + to_string(j - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(j);
            }
            cout << ss << ", ";
        }

        cout << "\n\n";
    }

    // 1 for max problem -1 for min problem
    void compute_table(ll min_or_max = 1)
    {
        ll itrrcount = 1;
        while (true)
        {
            ll new_etr = this->new_entering_variable(min_or_max);
            if (new_etr == -1)
            {
                cout << "Coudn't find entering variable, Simplex iteration completed\n\n";
                break;
            }

            string ss;
            if (new_etr > this->n_nbsc)
            {
                ss = "s_" + to_string(new_etr - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(new_etr);
            }
            cout << "Entering Variable at " << itrrcount << " iteration is " << ss << "\n";

            ll new_lev = this->new_leaving_variable(new_etr);
            if (new_lev == -1)
            {
                cout << "Coudn't find leaving variable\n\n";
                break;
            }

            ll new_levid = this->xids[new_lev];

            if (new_levid > this->n_nbsc)
            {
                ss = "s_" + to_string(new_levid - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(new_levid);
            }

            cout << "Leaving Variable at " << itrrcount << " iteration is " << ss << "\n\n";
            this->scale(new_lev, this->A[new_lev][new_etr]);

            for (ll i = 0; i < n_bsc + 1; i++)
            {
                if (i != new_lev)
                {
                    this->r_reduce(i, new_lev, this->A[i][new_etr]);
                }
            }

            this->xids[new_lev] = new_etr;
            if (printtbl)
            {
                cout << "After " << itrrcount << "th interation the table is:\n";
                cout << (*this);
            }
            this->printvar_category();
            itrrcount++;
        }
    }

    void input()
    {

        cout << "for each next " << this->n_bsc + 1 << " rows input " << this->n_tot + 2 << " space separtaed numbers in the format :\nz x[i][1] ... x[i][n_nbs] s[i][1] ... s[i][n_slc] Sol[i]\n";
        for (ll i = 0; i < n_bsc + 1; i++)
        {
            for (ll j = 0; j < n_tot + 1; j++)
            {
                cin >> this->A[i][j];
            }
            cin >> this->sol[i];
        }
    }

    friend ostream &operator<<(ostream &os, simplex &si)
    {
        string s = "Basic";
        os.width(wd);
        os << s;

        s = "z";
        os.width(wd);
        os << s;

        for (ll i = 0; i < si.n_nbsc; i++)
        {
            s = "x_";
            s += to_string(i + 1);
            os.width(wd);
            os << s;
        }
        for (ll i = 0; i < si.n_bsc; i++)
        {
            s = "s_";
            s += to_string(i + 1);
            os.width(wd);
            os << s;
        }
        os.width(wd);
        os << "Solution\n";

        for (ll i = 0; i < si.n_bsc + 1; i++)
        {
            if (si.xids[i] == 0)
            {
                s = "z";
            }
            else if (si.xids[i] > si.n_nbsc)
            {
                s = "s_" + to_string(si.xids[i] - si.n_nbsc);
            }
            else
            {
                s = "x_" + to_string(si.xids[i]);
            }

            os.width(wd);
            os << s;

            for (ll j = 0; j < si.n_tot + 1; j++)
            {
                os.width(wd);
                os << si.A[i][j];
            }

            os.width(wd);
            os << si.sol[i];
            os << "\n";
        }

        return os;
    }

    void print_solution()
    {
        cout << "Final solution is:\n";
        cout << "z=" << sol[0] << "\n";
        vector<R> vvl(this->n_tot + 1, 0);
        for (ll i = 1; i < n_bsc + 1; i++)
        {
            vvl[xids[i]] = sol[i];
        }

        string ss;
        for (ll i = 1; i < this->n_tot + 1; i++)
        {
            if (i > this->n_nbsc)
            {
                ss = "s_" + to_string(i - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(i);
            }
            cout << ss << "=" << vvl[i] << "\n";
        }
        cout << "\n\n";
    }
};

// function to compute the objective value
R objective(vector<R> &xs)
{
    R val = 0;
    for (ll i = 0; i < zcoeff.size(); i++)
    {
        val += zcoeff[i] * xs[i];
    }
    return val;
}

void standardize(vvr &AA, vr &BB, vector<char> &etyp, linear_eq &LE)
{
    ll sct = 0;
    ll ndv = AA[0].size();

    for (ll i = 0; i < etyp.size(); i++)
    {
        if (etyp[i] == '<')
        {
            for (ll j = 0; j < ndv; j++)
            {
                LE.A[i][j] = AA[i][j];
            }
            LE.A[i][ndv + sct] = (R)1;
            LE.B[i] = BB[i];
            sct++;
        }
        else if (etyp[i] == '=')
        {
            for (ll j = 0; j < ndv; j++)
            {
                LE.A[i][j] = AA[i][j];
            }
            LE.B[i] = BB[i];
        }
        else if (etyp[i] == '>')
        {
            for (ll j = 0; j < ndv; j++)
            {
                LE.A[i][j] = AA[i][j];
            }
            LE.A[i][ndv + sct] = (R)(-1);
            LE.B[i] = BB[i];
            sct++;
        }
    }
}

void init_simplex_table(linear_eq &LE, simplex &si)
{
    si.A[0][0] = 1;
    for (ll i = 1; i < si.n_nbsc + 1; i++)
    {
        si.A[0][i] = (-1) * zcoeff[i - 1];
    }

    for (ll i = 1; i < si.n_bsc + 1; i++)
    {
        for (ll j = 1; j < si.n_tot + 1; j++)
        {
            si.A[i][j] = LE.A[i - 1][j - 1];
        }

        si.sol[i] = LE.B[i - 1];
    }
}

struct game
{
    ll n, m;
    vvr payoff;
    vvr orgpayoff;
    ll which_player;
    simplex simplex_table;
    ll max_or_min;
    vr row_sol;
    vr col_sol;
    ll pure_row = -1;
    ll pure_col = -1;

    game(ll _n, ll _m)
    {
        this->n = _n;
        this->m = _m;
        this->which_player = -1;
    }

    void input()
    {
        cout << "Enter number of strategies for row player (player A) :";
        cin >> this->n;
        cout << "Enter number of strategies for column player (player B) :";
        cin >> this->m;

        this->row_sol.assign(this->n, 0);
        this->col_sol.assign(this->m, 0);

        this->payoff.assign(this->n, vr(this->m));
        this->orgpayoff.assign(this->n, vr(this->m));

        cout << "For next " << this->n << " rows enter " << this->m << " space saperated values that are payoffs in the format\n"
             << "p[i][0] ...p[i][m-1]\n\n";

        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                cin >> this->payoff[i][j];
            }
        }

        this->orgpayoff = this->payoff;
    }

    pair<R, R> pure_strategy_sols()
    {
        R row_val = -infinity_R, col_val = infinity_R;
        for (ll i = 0; i < this->n; i++)
        {
            R trow_val = infinity_R;
            for (ll j = 0; j < this->m; j++)
            {
                trow_val = min(trow_val, this->payoff[i][j]);
            }

            if (row_val <= trow_val)
            {
                row_val = trow_val;
                this->pure_row = i;
            }
        }

        for (ll i = 0; i < this->m; i++)
        {
            R tcol_val = -infinity_R;
            for (ll j = 0; j < this->n; j++)
            {
                tcol_val = max(tcol_val, this->payoff[j][i]);
            }

            if (tcol_val <= col_val)
            {
                col_val = tcol_val;
                this->pure_col = i;
            }
        }

        return {row_val, col_val};
    }

    // default convert wrt column player
    void convert_to_lpp(ll wrt_which_player = 1)
    {

        ll morm = 1, ne, nv, tnv, ndv, nc;
        zcoeff.clear();

        if (wrt_which_player == 1)
        {
            // column player
            morm = 1ll;
            this->max_or_min = 1ll;
            ne = this->n;
            nv = this->m;
            this->which_player = 1ll;
        }
        else
        {
            // row player
            morm = -1ll;
            this->max_or_min = -1ll;
            ne = this->m;
            nv = this->n;
            this->which_player = 0ll;
        }

        tnv = nv;
        zcoeff.assign(nv, (R)0);
        for (ll i = 0; i < nv; i++)
        {
            zcoeff[i] = 1;
        }

        nc = ne;
        vvr AA(nc, vr(nv));
        vr VAL(nc);
        vector<char> etyp(nc);

        for (ll i = 0; i < nc; i++)
        {
            for (ll j = 0; j < nv; j++)
            {
                if (wrt_which_player == 1)
                {
                    AA[i][j] = this->payoff[i][j];
                }
                else
                {
                    AA[i][j] = this->payoff[j][i];
                }
            }
            if (wrt_which_player == 1)
            {
                etyp[i] = '<';
                tnv++;
            }
            else
            {
                etyp[i] = '>';
                tnv++;
            }
            VAL[i] = 1;
        }

        linear_eq LE(ne, tnv);
        standardize(AA, VAL, etyp, LE); // standardize the problem
        cout << " ||   Standard Format of Problem is  ||\n";
        cout << LE;

        this->simplex_table.init(tnv - nv, nv);
        init_simplex_table(LE, this->simplex_table); // initilize simplex table

        return;
    }

    void solve_simplex()
    {
        cout << "Initial Simplex table:\n ";
        cout << this->simplex_table;

        printtbl = true;
        this->simplex_table.compute_table(this->max_or_min);
        cout << this->simplex_table;
        this->simplex_table.print_solution();
    }

    void get_game_sol()
    {

        R Z = this->simplex_table.sol[0];
        vector<R> vvl(this->simplex_table.n_tot + 1, 0);
        for (ll i = 1; i < this->simplex_table.n_bsc + 1; i++)
        {
            vvl[this->simplex_table.xids[i]] = this->simplex_table.sol[i];
            if (this->which_player == 1)
            {
                this->col_sol[this->simplex_table.xids[i] - 1] = this->simplex_table.sol[i] / Z;
            }
            else
            {
                this->row_sol[this->simplex_table.xids[i] - 1] = this->simplex_table.sol[i] / Z;
            }
        }

        if (this->which_player == 1)
        {
            for (ll j = this->simplex_table.n_nbsc + 1; j < this->simplex_table.n_tot + 1; j++)
            {
                // debug(j);
                this->row_sol[j - this->simplex_table.n_nbsc - 1] = this->simplex_table.A[0][j] / Z;
            }
        }
        else
        {
            for (ll j = this->simplex_table.n_nbsc + 1; j < this->simplex_table.n_tot + 1; j++)
            {
                this->col_sol[j - this->simplex_table.n_nbsc - 1] = this->simplex_table.A[0][j] / Z;
            }
        }

        cout << "\nOptimal mixed strategy for player A:\n";
        for (ll i = 0; i < this->row_sol.size(); i++)
        {
            cout << "p_" << i + 1 << " = " << this->row_sol[i] << ", ";
        }
        cout << "\n";

        cout << "Optimal mixed strategy for player B:\n";
        for (ll i = 0; i < this->col_sol.size(); i++)
        {
            cout << "p_" << i + 1 << " = " << this->col_sol[i] << ", ";
        }
        cout << "\n";

        R game_val = 0;

        for (ll j = 0; j < this->m; j++)
        {
            game_val += this->orgpayoff[0][j] * this->col_sol[j];
        }

        cout << "Value of game is: " << game_val << "\n";
    }

    void add_shift(R C)
    {
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                this->payoff[i][j] += C;
            }
        }
    }

    void solve(ll wrt_which_player = 1)
    {
        auto rc_sol = this->pure_strategy_sols();
        if (abs(rc_sol.first - rc_sol.second) <= _EPSILON)
        {
            cout << "\nGame is Stable\nHence, optimal solution is a pure Strategy with moves: \nRow player(A) move: " << this->pure_row + 1 << "\nColumn player(B) move: " << this->pure_col + 1 << "\nWith Value of the game being " << rc_sol.first << "\n\n";
            return;
        }

        cout << "\nGame is not Stable\noptimal solution is a mixed Strategy,\nand value of game lies between " << rc_sol.first << " and " << rc_sol.second << "\n";
        cout << "\n";

        R C = rc_sol.second + (R)2;
        cout << "Adding constant value " << C << " to all payoff values\n\n";
        // C = 0;

        this->add_shift(C);

        cout << "Converting game to LPP problem\n";
        this->convert_to_lpp(wrt_which_player);

        cout << "Solving game using simplex method\n";
        this->solve_simplex();

        cout << "Solution of the game is\n";
        this->get_game_sol();
    }
};

void lab_11()
{
    game g(0, 0);
    g.input();
    g.solve();
    return;
}

int main()
{
    lab_11();
    return 0;
}
