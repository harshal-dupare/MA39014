#include <iostream>
#include <vector>

using namespace std;
#define EPS 1e-5

typedef vector<double> vf;
typedef vector<vector<double>> vvf;

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

bool GS(int n, int k, vvf &A, vf &B, vf &x)
{
    for (int i = 0; i < n; i++)
        x[i] = 0;

    double err = 1;
    vf xl = x;

    while (k > 0 && err > EPS)
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

    double errlim = (double)n * 10.0;
    if (err > errlim)
    {
        return false;
    }
    return true;
}

int main()
{
    vvf a(100, vf(100));
    vf b(100), x(100);
    int n;
    cout << "Enter number of variables: ";
    cin >> n;
    int k = 10; // number of iteration
    cout << "for each next " << n << " rows output " << n + 1 << " space ssepartaed numbers as : a[i][1] a[i][2] .. a[i][n] b[i]\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    bool converges = GS(n, 10, a, b, x);
    if(converges) cout<<"Converges\n";
    else cout<<"Diverges\n";

    for (int i = 0; i < n; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    return 0;
}
