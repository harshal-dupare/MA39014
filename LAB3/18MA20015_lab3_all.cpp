#include <iostream>
#include <vector>

using namespace std;
typedef double R;
typedef vector<vector<R>> vvr;
typedef vector<R> vr;
int itrcount = 0;
int prtitr = 0;
int wd = 10;
bool printtbl = false;

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
        this->A.assign(_n_eq, vector<R>(_n_var, (R)0));
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

    void printvar_category()
    {
        vector<int> bsc(this->n_tot + 1, 0);
        vector<int> bvr, nbvr;
        for (auto j : this->xids)
        {
            bsc[j] = 1;
        }

        for (int i = 1; i < this->n_tot + 1; i++)
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
            if (j > this->n_nbsc)
            {
                ss = "s_" + to_string(j - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(j);
            }
            cout<<ss<<", ";
        }
        cout << "\n";

        cout << "Non Basic Variables Are :";
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
            cout<<ss<<", ";
        }

        cout << "\n\n";
    }

    // 1 for max problem -1 for min problem
    void compute_table(int min_or_max = 1)
    {
        int itrrcount = 1;
        while (true)
        {
            int new_etr = this->new_entering_variable(min_or_max);
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

            int new_lev = this->new_leaving_variable(new_etr);
            if (new_lev == -1)
            {
                cout << "Coudn't find leaving variable\n\n";
                break;
            }

            int new_levid = this->xids[new_lev];

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

            for (int i = 0; i < n_bsc + 1; i++)
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
        string s = "Basic";
        os.width(wd);
        os << s;

        s = "z";
        os.width(wd);
        os << s;

        for (int i = 0; i < si.n_nbsc; i++)
        {
            s = "x_";
            s += to_string(i + 1);
            os.width(wd);
            os << s;
        }
        for (int i = 0; i < si.n_bsc; i++)
        {
            s = "s_";
            s += to_string(i + 1);
            os.width(wd);
            os << s;
        }
        os.width(wd);
        os << "Solution\n";

        for (int i = 0; i < si.n_bsc + 1; i++)
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

            for (int j = 0; j < si.n_tot + 1; j++)
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
        for (int i = 1; i < n_bsc + 1; i++)
        {
            vvl[xids[i]] = sol[i];
        }

        string ss;
        for (int i=1;i<this->n_tot+1;i++)
        {
            if (i > this->n_nbsc)
            {
                ss = "s_" + to_string(i - this->n_nbsc);
            }
            else
            {
                ss = "x_" + to_string(i);
            }
            cout<<ss<<"="<<vvl[i]<<"\n";
        }
        cout << "\n\n";
    }
};

vector<R> zcoeff;

// function to compute the objective value
R objective(vector<R> &xs)
{
    R val = 0;
    for (int i = 0; i < zcoeff.size(); i++)
    {
        val += zcoeff[i] * xs[i];
    }
    return val;
}

// printing solutions of max problems
void print_solutions(vector<vector<R>> &bfs)
{
    int oid = 0;
    int id = 0;
    R mxv = -1.0;

    for (int i = 0; i < bfs[0].size() - 1; i++)
    {
        string s = "x_";
        s += to_string(i + 1);
        cout.width(12);
        cout << s;
    }
    cout.width(12);
    cout << "Converges?";
    cout.width(12);
    cout << "Feasible?";
    cout.width(12);
    cout << "Z";
    cout << "\n";
    for (auto s : bfs)
    {
        bool feasible = true;
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                cout.width(12);
                cout << ((s[i]) ? "Yes" : "No");
                continue;
            }
            if (s[i] < 0)
                feasible = false;
            cout.width(12);
            cout << s[i];
        }
        cout.width(12);
        cout << ((feasible) ? "YES" : "NO");
        R z_val = objective(s);
        cout.width(12);
        cout << z_val;
        cout << "\n";

        if (feasible && z_val > mxv)
        {
            mxv = z_val;
            oid = id;
        }
        id++;
    }

    cout << "\n Among them the optimal solution is \n";
    for (int i = 0; i < bfs[0].size() - 1; i++)
    {
        string s = "x_";
        s += to_string(i + 1);
        cout.width(12);
        cout << s;
    }
    cout.width(12);
    cout << "Converges?";
    cout.width(12);
    cout << "Feasible?";
    cout.width(12);
    cout << "Z";
    cout << "\n";
    bool feasible = true;
    for (int i = 0; i < bfs[oid].size(); i++)
    {
        if (i == bfs[oid].size() - 1)
        {
            cout.width(12);
            cout << ((bfs[oid][i]) ? "Yes" : "No");
            continue;
        }
        if (bfs[oid][i] < 0)
            feasible = false;
        cout.width(12);
        cout << bfs[oid][i];
    }

    cout.width(12);
    cout << (feasible ? "Yes" : "No");
    cout.width(12);
    cout << mxv;
    cout << "\n";
}

// printing solutions of min problems
void print_min_solutions(vector<vector<R>> &bfs)
{
    int oid = 0;
    int id = 0;
    R mxv = -1.0;

    for (int i = 0; i < bfs[0].size() - 1; i++)
    {
        string s = "x_";
        s += to_string(i + 1);
        cout.width(12);
        cout << s;
    }
    cout.width(12);
    cout << "Converges?";
    cout.width(12);
    cout << "Feasible?";
    cout.width(12);
    cout << "Z";
    cout << "\n";
    for (auto s : bfs)
    {
        bool feasible = true;
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                cout.width(12);
                cout << ((s[i]) ? "Yes" : "No");
                continue;
            }
            if (s[i] < 0)
                feasible = false;
            cout.width(12);
            cout << s[i];
        }
        cout.width(12);
        cout << ((feasible) ? "YES" : "NO");
        R z_val = objective(s);
        cout.width(12);
        cout << z_val;
        cout << "\n";

        if (feasible && z_val < mxv)
        {
            mxv = z_val;
            oid = id;
        }
        id++;
    }

    cout << "\n Among them the optimal solution is \n";
    for (int i = 0; i < bfs[0].size() - 1; i++)
    {
        string s = "x_";
        s += to_string(i + 1);
        cout.width(12);
        cout << s;
    }
    cout.width(12);
    cout << "Converges?";
    cout.width(12);
    cout << "Feasible?";
    cout.width(12);
    cout << "Z";
    cout << "\n";
    bool feasible = true;
    for (int i = 0; i < bfs[oid].size(); i++)
    {
        if (i == bfs[oid].size() - 1)
        {
            cout.width(12);
            cout << ((bfs[oid][i]) ? "Yes" : "No");
            continue;
        }
        if (bfs[oid][i] < 0)
            feasible = false;
        cout.width(12);
        cout << bfs[oid][i];
    }

    cout.width(12);
    cout << (feasible ? "Yes" : "No");
    cout.width(12);
    cout << mxv;
    cout << "\n";
}

void standardize(vvr &AA, vr &BB, vector<char> &etyp, linear_eq &LE)
{
    int sct = 0;
    int ndv = AA[0].size();

    for (int i = 0; i < etyp.size(); i++)
    {
        if (etyp[i] == '<')
        {
            for (int j = 0; j < ndv; j++)
            {
                LE.A[i][j] = AA[i][j];
            }
            LE.A[i][ndv + sct] = (R)1;
            LE.B[i] = BB[i];
            sct++;
        }
        else if (etyp[i] == '=')
        {
            for (int j = 0; j < ndv; j++)
            {
                LE.A[i][j] = AA[i][j];
            }
            LE.B[i] = BB[i];
        }
        else if (etyp[i] == '>')
        {
            for (int j = 0; j < ndv; j++)
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
    for (int i = 1; i < si.n_nbsc + 1; i++)
    {
        si.A[0][i] = (-1) * zcoeff[i - 1];
    }

    for (int i = 1; i < si.n_bsc + 1; i++)
    {
        for (int j = 1; j < si.n_tot + 1; j++)
        {
            si.A[i][j] = LE.A[i - 1][j - 1];
        }

        si.sol[i] = LE.B[i - 1];
    }
}

void lab_3()
{
    zcoeff.clear();
    int morm = 1, ne, nv, tnv, ndv, nc;
    cout << "Enter 1 is its maximization problem if minimization problem then -1:";
    cin >> morm;
    cout << "number of variables:";
    cin >> nv;
    tnv = nv;
    zcoeff.assign(nv, (R)0);
    cout << "Enter the coefficients of variables in objective\nif cx[i] is coefficient of x_i then give input in format\ncx[1] cx[2] ... cx[n]\n";
    for (int i = 0; i < nv; i++)
    {
        cin >> zcoeff[i];
    }
    cout << "Enter number of conditions: ";
    cin >> nc;
    ne = nc;
    vvr AA(nc, vr(nv));
    vr VAL(nc);
    vector<char> etyp(nc);

    cout << "For each condition give input of the coefficient of the variables in condition\n "
         << "for greater than equal to condition: coeff[i][0] .... coeff[i][number_variable] > val[i]\n "
         << "for              equal to condition: coeff[i][0] .... coeff[i][number_variable] = val[i]\n "
         << "for    less than equal to condition: coeff[i][0] .... coeff[i][number_variable] < val[i]\n";
    for (int i = 0; i < nc; i++)
    {
        for (int j = 0; j < nv; j++)
        {
            cin >> AA[i][j];
        }
        cin >> etyp[i];
        cin >> VAL[i];
        if (etyp[i] == '<')
        {
            tnv++;
        }
        else if (etyp[i] == '=')
        {
        }
        else if (etyp[i] == '>')
        {
            tnv++;
        }
    }

    linear_eq LE(ne, tnv);
    standardize(AA, VAL, etyp, LE); // standardize the problem
    cout << " ||   Standard Format of Problem is  ||\n";
    cout << LE;

    simplex si(tnv - nv, nv);
    init_simplex_table(LE, si); // initilize simplex table

    cout << "Initial Simplex table:\n ";
    cout << si;

    printtbl = true;
    si.compute_table(morm);
    cout << si;
    si.print_solution();

    int TT = 1;
    while (true && TT > 0)
    {
        int casse;
        cout << "\n\n|| Menu ||\n"
             << " for Exit, enter : 0\n"
             << " for list the BFS, enter : 1 \n"
             << " for Number of Iterations to solve the problem, enter : 2\n"
             << " for List of all Non-basic variables along with net evaluations in i'th(user input) iteration, enter : 3\n"
             << " for List of Basic variables along with min ratios inithiteration, enter : 4\n"
             << " for optimal solution, enter : 5 \n"
             << "\nEnter Menu option: ";
        // cin >> casse;
        cout << "\n";
        casse = 1;
        if (casse == 0)
        {
            cout << "\n||     Thankyou     ||\n";
            break;
        }
        else if (casse == 1)
        {
            linear_eq le = LE;
            OR_lab orl;
            vector<vector<R>> bfs;
            orl.basic_feasible_soln(bfs, le);

            cout << "BFS's are:\n";
            if (morm == 1)
            {
                print_solutions(bfs);
            }
            else
            {
                print_min_solutions(bfs);
            }
            cout << "\n";
        }
        else if (casse == 2)
        {
            int ith;
            cout << "Enter the number of iteration to solve the problem:";
            cin >> ith;
            itrcount = ith;
        }
        else if (casse == 3)
        {
        }
        else if (casse == 4)
        {
        }
        else if (casse == 5)
        {
        }
        else
        {
            cout << "            WRONG INPUT : TRY AGAIN \n";
        }

        TT--;
    }

    return;
}

int main()
{
    lab_3();
    return 0;
}