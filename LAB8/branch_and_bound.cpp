/*
18MA20015 | Harshal Dupare | Lab 5 | 10-2-2021 |
*/

#include <iostream>
#include <vector>
#include "library/debug.h"

using namespace std;
typedef double R;
typedef long long ll;
typedef vector<vector<R>> vvr;
typedef vector<R> vr;
ll itrcount = 0;
ll wd = 14;
R _EPSILON = 1e-4;
R infinity_R = 1e15;

// to store the coefficients of objective function
vector<R> zcoeff;

R fractional_part(R x)
{
    R x_ssc = (R)(((ll)(x)) + 3ll);
    R x_sc = (R)(((ll)(x)) + 2ll);
    R x_c = (R)(((ll)(x)) + 1ll);
    R x_f = (R)((ll)(x));
    R x_s = (R)(((ll)(x)) - 1ll);
    R x_sf = (R)(((ll)(x)) - 2ll);
    R x_ssf = (R)(((ll)(x)) - 3ll);

    for(ll i=-2;i<3;i++)
    {
        R x_i = (R)(i);
        if(x==x_i+_EPSILON||x==x_i-_EPSILON)
        {
            return 0;
        }
    }

    if (abs(x_c - x) <= _EPSILON || abs(x_f - x) <= _EPSILON || abs(x-x_s)<=_EPSILON || abs(x-x_sc)<=_EPSILON || abs(x-x_sf)<=_EPSILON || abs(x-x_ssf)<=_EPSILON || abs(x-x_ssc)<=_EPSILON)
    {
        return 0;
    }

    if (x < 0)
    {
        return 1 - x_f + x;
    }
    else
    {
        return x - x_f;
    }

    return 0;
}

R integral_part(R x)
{
    R x_sc = (R)(((ll)(x)) + 2ll);
    R x_c = (R)(((ll)(x)) + 1ll);
    R x_f = (R)((ll)(x));
    R x_s = (R)(((ll)(x)) - 1ll);
    R x_sf = (R)(((ll)(x)) - 2ll);

    for(ll i=-2;i<3;i++)
    {
        R x_i = (R)(i);
        if(x==x_i+_EPSILON||x==x_i-_EPSILON)
        {
            return x_i;
        }
    }

    // cout<<x_c<<":"<<x_f<<endl;
    if(abs(x)<=_EPSILON)
    {
        return 0;
    }

    if (abs(x_f - x) <= _EPSILON)
    {
        return x_f;
    }

    if (abs(x_c - x) <= _EPSILON)
    {
        return x_c;
    }

    if (abs(x_sc - x) <= _EPSILON)
    {
        return x_sc;
    }

    if (abs(x_sf - x) <= _EPSILON)
    {
        return x_sf;
    }

    if (abs(x_s - x) <= _EPSILON)
    {
        return x_s;
    }

    if (x > 0)
    {
        return x_f;
    }
    else if (x < 0)
    {
        return x_f - 1;
    }

    return 0;
}

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
};

struct dual_simplex
{
    vector<vector<R>> A;
    vector<R> sol;
    vector<R> z;
    vector<ll> xids;
    ll n_dv;
    ll n_slc;
    ll n_arf;
    ll n_srp;
    ll n_tot;
    ll n_nbsc;
    ll n_bsc;
    R M = 1e5;
    R EPS = 1e-6;
    R _pEPS = 1e-9;
    R ninf = -1e20;
    R inf = 1e20;
    bool printtbl = true;
    ll min_or_max = 1;

    void log_details()
    {
        cout << "n_dv : " << n_dv << "\n";
        cout << "n_slc : " << n_slc << "\n";
        cout << "n_arf : " << n_arf << "\n";
        cout << "n_srp : " << n_srp << "\n";
        cout << "n_tot : " << n_tot << "\n";
        cout << "n_nbsc : " << n_nbsc << "\n";
        cout << "n_bsc : " << n_bsc << "\n";
        cout << "A.size() : " << A.size() << "\n";
        cout << "A[0].size() : " << A[0].size() << "\n";
        cout << "sol.size() : " << sol.size() << "\n";
        cout << "xids.size() : " << xids.size() << "\n";
    }

    dual_simplex(ll _n_dv, ll _n_slc, ll _n_srp, ll _n_arf)
    {
        this->n_dv = _n_dv;
        this->n_slc = _n_slc;
        this->n_srp = _n_srp;
        this->n_arf = _n_arf;
        this->n_nbsc = _n_slc + _n_arf;
        this->n_bsc = _n_srp + _n_dv;
        this->n_tot = _n_dv + _n_slc + _n_srp + _n_arf;
        this->A.assign(_n_slc + _n_arf + 1, vector<R>(this->n_tot + 1));
        this->sol.assign(_n_slc + _n_arf + 1, 0);
        this->xids.assign(_n_slc + _n_arf + 1, 0);
        for (ll i = 1; i < _n_slc + _n_arf + 1; i++)
        {
            this->xids[i] = i + _n_srp + _n_dv;
        }
    }

    string id_to_var(ll id)
    {
        if (id == 0)
        {
            return "z";
        }

        if (id < this->n_dv + 1)
        {
            return "x_" + to_string(id);
        }
        else if (id < this->n_dv + this->n_srp + 1)
        {
            return "srp_" + to_string(id - this->n_dv);
        }
        else if (id < this->n_dv + this->n_arf + this->n_srp + 1)
        {
            return "art_" + to_string(id - (this->n_dv + this->n_srp));
        }
        else
        {
            return "s_" + to_string(id - (this->n_dv + this->n_srp + this->n_arf));
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

    // makes the table consistant
    void make_consistant()
    {
        for (ll i = this->n_bsc + 1; i < this->n_bsc + this->n_arf + 1; i++)
        {
            this->r_reduce(0ll, i - this->n_bsc, (this->min_or_max) * this->M);
        }
    }

    // optamilaity condition
    ll new_entering_variable(ll k)
    {
        R mnv = inf;
        ll mid = -1;
        for (ll i = 1; i < n_tot + 1; i++)
        {
            if ( this->A[k][i] < 0.0 &&  abs(this->A[0][i]/this->A[k][i]) < mnv)
            {
                mnv = abs(this->A[0][i]/this->A[k][i]) ;
                mid = i;
            }
        }

        return mid;
    }

    // feasibility condiion
    ll new_leaving_variable()
    {
        ll mnid = -1;
        R mnv = inf;
        for (ll i = 1; i < n_nbsc + 1; i++)
        {
            // most negative value for leaving - minimization
            if (this->sol[i] < 0 && this->sol[i] < mnv)
            {
                mnv = this->sol[i];
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
        cout << "Basic Variables Are :";
        for (auto j : bvr)
        {
            cout << this->id_to_var(j) << ", ";
        }
        cout << "\n";

        cout << "Non Basic Variables Are :";
        for (auto j : nbvr)
        {
            cout << this->id_to_var(j) << ", ";
        }

        cout << "\n";
    }

    // 1 for max problem -1 for min problem
    void compute_table()
    {
        ll itrrcount = 1;
        while (true)
        {
            ll new_lev = this->new_leaving_variable();
            if (new_lev == -1)
            {
                if (this->printtbl)
                    cout << "Coudn't find leaving variable, Simplex iteration completed\n\n";
                break;
            }

            if (this->printtbl)
                cout << "Leaving Variable at " << itrrcount << " iteration is " << this->id_to_var(this->xids[new_lev]) << "\n";

            ll new_etr = this->new_entering_variable(new_lev);
            if (new_etr == -1)
            {
                if (this->printtbl)
                    cout << "Coudn't find entering variable\n\n";
                break;
            }

            if (this->printtbl)
                cout << "Entering Variable at " << itrrcount << " iteration is " << new_etr << "\n\n";

            this->scale(new_lev, this->A[new_lev][new_etr]);

            for (ll i = 0; i < n_nbsc + 1; i++)
            {
                if (i != new_lev)
                {
                    this->r_reduce(i, new_lev, this->A[i][new_etr]);
                }
            }

            this->xids[new_lev] = new_etr;
            if (printtbl)
            {
                if (this->printtbl)
                    cout << "After " << itrrcount << "th interation the table is:\n";
                cout << (*this);
            }
            if (this->printtbl)
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

    friend ostream &operator<<(ostream &os, dual_simplex &si)
    {
        string s = "Basic";
        os.width(wd);
        os << s;

        for (ll i = 0; i < si.n_tot + 1; i++)
        {
            os.width(wd);
            os << si.id_to_var(i);
        }
        os.width(wd);
        os << "Solution\n";

        for (ll i = 0; i < si.n_nbsc + 1; i++)
        {
            os.width(wd);
            os << si.id_to_var(si.xids[i]);

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

        bool feasible = true;
        vector<R> vvl(this->n_tot + 1, 0);
        for (ll i = 1; i < n_nbsc + 1; i++)
        {
            vvl[xids[i]] = sol[i];
            if (xids[i] >= this->n_bsc + 1 && xids[i] < this->n_bsc + this->n_arf + 1)
            {
                feasible = false;
            }
        }

        if (!feasible)
        {
            cout << "Probelem has NO SOLUTION / INFEASIBLE SOLUTION as there is/are artificial variables is in the non basic variables\n\n";
            return;
        }

        cout << "z=" << sol[0] << "\n";
        string ss;
        for (ll i = 1; i < this->n_tot + 1; i++)
        {
            ss = this->id_to_var(i);
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

vr optimal_solutions(ll mxormn, vector<vector<R>> &bfs)
{
    int oid = 0;
    int id = 0;
    R ov=infinity_R;
    if(mxormn==1)
    {
        ov = (-1)*infinity_R;
    }

    vr ox;
    for (auto s : bfs)
    {
        bool converges = true;
        bool feasible = true;
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                if(s[i]==0)
                {
                    converges = false;
                }
            }
            else if (s[i] < 0)
            {
                feasible = false;
            }
        }

        if(feasible==false||converges==false)
        {
            continue;
        }

        R z_val = objective(s);

        if(mxormn==1 && z_val > ov)
        {
            ox = s;
            ov = z_val;
        }

        if(mxormn==-1 && z_val < ov)
        {
            ox = s;
            ov = z_val;
        }
        
    }

    for(ll i = 0;i<min(zcoeff.size(),ox.size());i++)
    {
        cout<<"x_"<<i+1<<" = "<<ox[i]<<", ";
    }
    if(ox.size()>0)
    {
        cout<<"Z="<<objective(ox)<<endl;
    }
    else
    {
        cout<<"Infeasible\n";
    }
    return ox;
}

// to standardize the table
void standardize(vvr &AA, vr &BB, vector<char> &etyp, linear_eq &LE, ll tsrp)
{
    ll sct = 0;
    ll srct = 0;
    ll ndv = AA[0].size();

    for (ll i = 0; i < etyp.size(); i++)
    {
        for (ll j = 0; j < ndv; j++)
        {
            LE.A[i][j] = AA[i][j];
        }
        LE.B[i] = BB[i];
        if (etyp[i] == '<')
        {
            LE.A[i][ndv + sct + tsrp] = (R)1;
            sct++;
        }
        else if (etyp[i] == '=')
        {
            LE.A[i][ndv + sct + tsrp] = (R)1;
            sct++;
        }
        else if (etyp[i] == '>')
        {
            LE.A[i][ndv + sct + tsrp] = (R)(1);
            LE.A[i][ndv + srct] = (R)(-1);
            srct++;
            sct++;
        }
    }
}

// to initilize the big-M table
void init_dual_simplex_table(linear_eq &LE, dual_simplex &bm)
{
    bm.A[0][0] = 1;

    for (ll i = 1; i < bm.n_dv + 1; i++)
    {
        bm.A[0][i] = (-1) * zcoeff[i - 1];
    }

    for (ll i = bm.n_bsc + 1; i < bm.n_bsc + bm.n_arf + 1; i++)
    {
        bm.A[0][i] = 0;
    }

    for (ll i = 1; i < bm.n_nbsc + 1; i++)
    {
        for (ll j = 1; j < bm.n_tot + 1; j++)
        {
            bm.A[i][j] = LE.A[i - 1][j - 1];
        }

        bm.sol[i] = LE.B[i - 1];
    }
}

void convert_to_dual(ll morm, vvr &A, vvr &DA, vr &B, vr &DB, vector<char> &et, vector<char> &Det, vr &zc, vr &dzc, vector<ll> &vsg, vector<ll> &Dvsg)
{
    for (ll i = 0; i < A.size(); i++)
    {
        for (ll j = 0; j < A[0].size(); j++)
        {
            DA[j][i] = A[i][j];
        }
    }

    for (ll i = 0; i < A[0].size(); i++)
    {
        DB[i] = zc[i];
        if (morm == 1)
        {
            if (vsg[i] == 1)
                Det[i] = '>';
            else if (vsg[i] == -1)
                Det[i] = '<';
            else
                Det[i] = '=';
        }
        else
        {
            if (vsg[i] == 1)
                Det[i] = '<';
            else if (vsg[i] == -1)
                Det[i] = '>';
            else
                Det[i] = '=';
        }
    }

    for (ll i = 0; i < dzc.size(); i++)
    {
        dzc[i] = B[i];
        if (morm == 1)
        {
            if (et[i] == '>')
                Dvsg[i] = -1;
            else if (et[i] == '<')
                Dvsg[i] = 1;
            else
                Dvsg[i] = 0;
        }
        else
        {
            if (et[i] == '>')
                Dvsg[i] = 1;
            else if (et[i] == '<')
                Dvsg[i] = -1;
            else
                Dvsg[i] = 0;
        }
    }

    for (ll i = 0; i < dzc.size(); i++)
    {
        if (Dvsg[i] == -1)
        {
            dzc[i] *= (-1);
            for (ll j = 0; j < DA.size(); j++)
            {
                DA[j][i] *= (-1);
            }
            Dvsg[i] = 1;
        }
    }
}

void print_LPP(ll mxormn, vvr &A, vr &B, vector<char> &et, vector<ll> &vsgn, vr &zc)
{
    ll n = zc.size();
    ll m = A.size();
    if (mxormn == 1)
    {
        cout << "Maximize \t ";
    }
    else
    {
        cout << "Minimize \t ";
    }
    cout << "Z = ";
    for (ll i = 0; i < n; i++)
    {
        if (zc[i] != 0)
        {
            cout << std::showpos << zc[i] << std::noshowpos;
            cout << "x_" << (i + 1) << " ";
        }
    }
    cout << endl
         << "Subject to:\n";

    for (ll i = 0; i < m; i++)
    {
        cout << "\t";
        for (ll j = 0; j < n; j++)
        {
            if (A[i][j] != 0)
            {
                cout << std::showpos << A[i][j] << std::noshowpos;
                cout << "x_" << (j + 1) << " ";
            }
        }
        cout << et[i] << "= " << B[i] << endl;
    }

    for (ll i = 0; i < n; i++)
    {
        if (vsgn[i] == 1)
        {
            cout << "x_" << (i + 1) << " >= 0,";
        }
        else if (vsgn[i] == 0)
        {
            cout << "x_" << (i + 1) << " is unristrected, ";
        }
        else
        {
            cout << "x_" << (i + 1) << " <= 0,";
        }
    }
    cout << endl
         << endl;
}

vr solve(ll mxormn, vvr &A, vr &B, vector<char> &et, vector<ll> &vsgn, vr &zc)
{
    ll ne = B.size(), nv = A[0].size(), tnv, ndv, nc= B.size();
    ll tsrp = 0, tarf = 0, tslc = 0;
    tnv = nv;
    ndv = nv;

    for (ll i = 0; i < nc; i++)
    {
        if (et[i] == '<')
        {
            tnv++;
            tslc++;
        }
        else if (et[i] == '=')
        {
            tnv++;
            tarf++;
        }
        else if (et[i] == '>')
        {
            tnv++;
            tsrp++;
            tarf++;
        }
    }

    linear_eq LE(ne, tnv + tsrp);

    standardize(A, B, et, LE, tsrp);

    linear_eq le = LE;
    // debug(LE);
    OR_lab orl;
    vector<vector<R>> bfs;
    orl.basic_feasible_soln(bfs, le);

    // debug(bfs);

    vr ox = optimal_solutions(mxormn,bfs);

    return ox;
}

vr branch_and_bound(ll mxormn, vvr &A, vr &B, vector<char> &et, vector<ll> &vsgn, vr &zc)
{
    print_LPP(mxormn, A, B, et, vsgn, zc);
    
    ll n = vsgn.size();
    vr x_i = solve(mxormn, A, B, et, vsgn, zc);
    // if(A.size()>3) return x_i;
    ll frc_id = -1;

    for(ll i=0;i<min((ll)x_i.size(),n);i++)
    {
        if(abs(fractional_part(x_i[i]))>_EPSILON)
        {
            frc_id = i;
            break;
        }
    }

    if(frc_id==-1)
    {
        return x_i;
    }

    A.push_back(vr(n));
    A[A.size()-1][frc_id]=(R)1.0;
    B.push_back(integral_part(x_i[frc_id]));
    et.push_back('<');

    vr lx_i = branch_and_bound(mxormn, A, B, et, vsgn, zc);

    et[et.size()-1]='<';
    B[B.size()-1]+=(R)1.0;
    B[B.size()-1]*=(R)(-1.0);
    A[A.size()-1][frc_id]=(R)(-1.0);
    vr gx_i = branch_and_bound(mxormn, A, B, et, vsgn, zc);

    if(gx_i.size()==0)
    {
        return lx_i;
    }

    if(lx_i.size()==0)
    {
        return gx_i;
    }

    R gZ = objective(gx_i); 
    R lZ = objective(lx_i); 

    if(mxormn==1)
    {
        if(gZ > lZ)
        {
            return gx_i;
        }

        return lx_i;
    }
    else
    {
        if(gZ < lZ)
        {
            return gx_i;
        }

        return lx_i;
    }
}

void lab()
{
    zcoeff.clear();

    ll morm = 1, ne, nv, tnv, ndv, nc;
    ll tsrp = 0, tarf = 0, tslc = 0;

    cout << "Enter 1 is its maximization problem if minimization problem then -1:";
    cin >> morm;
    cout << "number of variables:";
    cin >> nv;
    tnv = nv;
    ndv = nv;

    zcoeff.assign(nv, (R)0);

    cout << "Enter the coefficients of variables in objective\nif cx[i] is coefficient of x_i then give input in format\ncx[1] cx[2] ... cx[n]\n";
    for (ll i = 0; i < nv; i++)
    {
        cin >> zcoeff[i];
    }

    cout << "Enter number of conditions: ";
    cin >> nc;
    ne = nc;

    vector<ll> vsg(nv, 1);

    vvr AA(nc, vr(nv));
    vr VAL(nc);

    vector<char> etyp(nc);

    cout << "For each condition give input of the coefficient of the variables in condition\n "
         << "for greater than equal to condition: coeff[i][0] .... coeff[i][number_variable] > val[i]\n "
         << "for              equal to condition: coeff[i][0] .... coeff[i][number_variable] = val[i]\n "
         << "for    less than equal to condition: coeff[i][0] .... coeff[i][number_variable] < val[i]\n";

    for (ll i = 0; i < nc; i++)
    {
        for (ll j = 0; j < nv; j++)
        {
            cin >> AA[i][j];
        }
        cin >> etyp[i];
        cin >> VAL[i];
        if (etyp[i] == '<')
        {
            tnv++;
            tslc++;
        }
        else if (etyp[i] == '=')
        {
            tnv++;
            tarf++;
        }
        else if (etyp[i] == '>')
        {
            tnv++;
            tsrp++;
            tarf++;
        }
    }

    vr sln_x = branch_and_bound(morm, AA, VAL, etyp, vsg, zcoeff);

    for(ll i = 0;i<min(zcoeff.size(),sln_x.size());i++)
    {
        cout<<"x_"<<i+1<<" = "<<sln_x[i]<<", ";
    }
    if(sln_x.size()==0)
    {
        cout<<"NO solution\n";
    }

    return;

    print_LPP(morm, AA, VAL, etyp, vsg, zcoeff);

    linear_eq LE(ne, tnv + tsrp);

    standardize(AA, VAL, etyp, LE, tsrp); // standardize the problem

    cout << " ||   Standard Format of Problem is  ||\n";
    cout << LE;
    
    dual_simplex bm(ndv, tslc, tsrp, tarf);
    bm.M = 1;
    bm.min_or_max = morm;

    init_dual_simplex_table(LE, bm); // initilize 2-phase simplex table
    cout << bm;

    // bm.log_details();

    bm.compute_table();
    cout << bm;

    bm.print_solution();
    return;
}

int main()
{
    lab();
    return 0;
}

/*
1
2
5 4
2
1 1 < 5
10 6 < 45
*/
