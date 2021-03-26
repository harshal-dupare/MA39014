/*

Harshal Dupare
18MA20015
OR LAB 13 Jan

*/

#include <iostream>
#include <vector>

#define EPS 1e-4

using namespace std;

typedef long double ldouble;
typedef vector<ldouble> vf;
typedef vector<int> vi;
typedef vector<vector<ldouble>> vvf;

// one itr of GS
void GS_itr(int n, vvf &A, vf &B, vf &x)
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

// GS
bool GS(int n, int k, vvf &A, vf &B, vf &x)
{
    for (int i = 0; i < n; i++)
        x[i] = 0;
    vf xl = x;
    ldouble err = 1.0;
    while (k > 0)
    {
        GS_itr(n, A, B, x);
        err = 0;
        for (int j = 0; j < n; j++)
        {
            err += abs(xl[j] - x[j]);
        }
        xl = x;
        k--;
    }

    ldouble errlim = (ldouble)n * 10.0;
    if (err > errlim)
    {
        return false;
    }
    return true;
}

// goes through all subset of size m and applies GS and saves soln
void compute_bfs(vi &subset, int lst, int n, int m, int k, vvf &A, vf &B, vvf &bfs)
{
    int total = 0;
    for (int i = 0; i < n; i++)
    {
        if (subset[i] == 1)
            total++;
    }

    if (total < m)
    {
        for (int i = lst; i < n; i++)
        {
            if (subset[i] == 0)
            {
                subset[i] = 1;
                compute_bfs(subset, i + 1, n, m, k, A, B, bfs);
                subset[i] = 0;
            }
        }
        return;
    }

    vvf a(m, vf(m));
    vf b(m), x(m, 0);

    for (int i = 0; i < m; i++)
    {
        int jd = 0;
        for (int j = 0; j < n; j++)
        {
            if (subset[j] == 1)
            {
                a[i][jd] = A[i][j];
                jd++;
            }
        }
        b[i] = B[i];
    }

    bool converge = GS(m, k, a, b, x);
    int id = 0;
    vf xx(n, 0);
    for (int i = 0; i < n; i++)
    {
        if (subset[i] == 1)
        {
            xx[i] = x[id];
            id++;
        }
    }

    // save soln
    if (converge)
        xx.push_back(1);
    else
        xx.push_back(0);

    bfs.push_back(xx);
}

int main()
{
    vvf a(100, vf(100)), bfs;
    vf b(100), x(100);
    int n, m;
    cout << "Enter number of variables: ";
    cin >> n;
    cout << "Enter number of equation: ";
    cin >> m;
    int k = 15000; // number of iteration might need to change to get the right answer
    cout << "for each next " << m << " rows output " << n + 1 << " space ssepartaed numbers as : a[i][1] a[i][2] .. a[i][n] b[i]\n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    vi subset(n, 0);
    compute_bfs(subset, 0, n, m, k, a, b, bfs);

    for (int i = 0; i < bfs.size(); i++)
    {
        bool feas = true;
        for (int j = 0; j < n; j++)
        {
            cout << bfs[i][j] << " ";
            if (bfs[i][j] < 0.0)
            {
                feas = false;
            }
        }
        if (bfs[i][n])
        {
            cout << " CONVERGES";
        }
        else
        {
            cout << " DIVERGES";
        }
        if (feas)
        {
            cout << " FEASIBLE" << endl;
        }
        else
        {
            cout << " NOT FEASIBLE" << endl;
        }
    }

    return 0;
}
