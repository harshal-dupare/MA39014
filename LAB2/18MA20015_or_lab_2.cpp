#include <iostream>
#include <vector>

using namespace std;
typedef double R;

struct linear_eq
{
    // x1+x2=3
    // 2x1+3x2=5
    // A = [[1,1],[2,3]]
    // B = [3,5]
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

// objective values for different problem
vector<R> zcoeff1 = {2, 5};
vector<R> zcoeff2 = {4, 3, 6};
vector<R> zcoeff3 = {12, 15, 14};
vector<R> zcoeff4 = {1, -3, 3};
vector<R> zcoeff5 = {3, 2, 2};

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

    for(int i=0;i<bfs[0].size()-1;i++)
    {
        string s="x_";
        s+=to_string(i+1);
        cout.width(12);
        cout<<s;
    }
    cout.width(12);
    cout<<"Converges?";
    cout.width(12);
    cout<<"Feasible?";
    cout.width(12);
    cout<<"Z";
    cout<<"\n";
    for (auto s : bfs)
    {
        bool feasible = true;
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                cout.width(12);
                cout<<((s[i])?"Yes":"No");
                continue;
            }
            if (s[i] < 0)
                feasible = false;
            cout.width(12);
            cout << s[i];
        }
        cout.width(12);
        cout<<((feasible)?"YES":"NO");
        R z_val = objective(s);
        cout.width(12);
        cout << z_val;
        cout << "\n";

        if(feasible&&z_val>mxv)
        {
            mxv=z_val;
            oid=id;
        }
        id++;
    }

    cout << "\n Among them the optimal solution is \n";
    for(int i=0;i<bfs[0].size()-1;i++)
    {
        string s="x_";
        s+=to_string(i+1);
        cout.width(12);
        cout<<s;
    }
    cout.width(12);
    cout<<"Converges?";
    cout.width(12);
    cout<<"Feasible?";
    cout.width(12);
    cout<<"Z";
    cout<<"\n";
    bool feasible = true;
    for (int i = 0; i < bfs[oid].size(); i++)
    {
        if (i == bfs[oid].size() - 1)
        {
            cout.width(12);
            cout<<((bfs[oid][i])?"Yes":"No");
            continue;
        }
        if (bfs[oid][i] < 0)
            feasible = false;
        cout.width(12);
        cout << bfs[oid][i];
    }

    cout.width(12);
    cout<<(feasible?"Yes":"No");
    cout.width(12);
    cout<<mxv;
    cout<<"\n";

}

// printing solutions of min problems
void print_min_solutions(vector<vector<R>> &bfs)
{
    int oid = 0;
    int id = 0;
    R mxv = -1.0;

    for(int i=0;i<bfs[0].size()-1;i++)
    {
        string s="x_";
        s+=to_string(i+1);
        cout.width(12);
        cout<<s;
    }
    cout.width(12);
    cout<<"Converges?";
    cout.width(12);
    cout<<"Feasible?";
    cout.width(12);
    cout<<"Z";
    cout<<"\n";
    for (auto s : bfs)
    {
        bool feasible = true;
        for (int i = 0; i < s.size(); i++)
        {
            if (i == s.size() - 1)
            {
                cout.width(12);
                cout<<((s[i])?"Yes":"No");
                continue;
            }
            if (s[i] < 0)
                feasible = false;
            cout.width(12);
            cout << s[i];
        }
        cout.width(12);
        cout<<((feasible)?"YES":"NO");
        R z_val = objective(s);
        cout.width(12);
        cout << z_val;
        cout << "\n";

        if(feasible&&z_val<mxv)
        {
            mxv=z_val;
            oid=id;
        }
        id++;
    }

    cout << "\n Among them the optimal solution is \n";
    for(int i=0;i<bfs[0].size()-1;i++)
    {
        string s="x_";
        s+=to_string(i+1);
        cout.width(12);
        cout<<s;
    }
    cout.width(12);
    cout<<"Converges?";
    cout.width(12);
    cout<<"Feasible?";
    cout.width(12);
    cout<<"Z";
    cout<<"\n";
    bool feasible = true;
    for (int i = 0; i < bfs[oid].size(); i++)
    {
        if (i == bfs[oid].size() - 1)
        {
            cout.width(12);
            cout<<((bfs[oid][i])?"Yes":"No");
            continue;
        }
        if (bfs[oid][i] < 0)
            feasible = false;
        cout.width(12);
        cout << bfs[oid][i];
    }

    cout.width(12);
    cout<<(feasible?"Yes":"No");
    cout.width(12);
    cout<<mxv;
    cout<<"\n";

}


void lab_2()
{
    int N=1;
    for(int k=0;k<N;k++)
    {
        zcoeff.clear();
        int ne, nv,ndv;
        cout<<"Give the inputs for "<<k+1<<"'th problem\n";
        cout<<"number of equations:";
        cin >> ne;
        cout<<"number of variables:";
        cin >> nv;
        cout<<"number of decision variables:";
        cin >> ndv;
        zcoeff.assign(ndv,(R)0);
        cout<<"Enter the number of coefficients of decision variables in problem\n if cx[i] is coefficient of x_i then give input in format, n= number of decision variable:\ncx[1] cx[2] ... cxn[0]\n";
        for(int i=0;i<ndv;i++)
        {
            cin>>zcoeff[i];
        }

        int morm = 1;
        cout<<"Enter 1 is its maximization problem if minimization problem then 0:";
        cin>>morm;

        linear_eq le(ne, nv);
        le.input();

        OR_lab orl;
        vector<vector<R>> bfs;
        orl.basic_feasible_soln(bfs, le);

        cout<<"BFS's are:\n";

        // for maximization problems
        // in 2nd input just give the problem number
        if(morm==1) print_solutions(bfs);

        // for minimization problems
        if(morm==0) print_min_solutions(bfs);
    }
    
}

int main()
{
    lab_2();
    return 0;
}