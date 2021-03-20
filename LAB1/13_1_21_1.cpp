/*

Harshal Dupare
18MA20015
OR LAB 13 Jan

*/

#include <iostream>
#include <vector>

using namespace std;

typedef vector<double> vf;
typedef vector<vector<double>> vvf;

void GS_itr(int n, vvf &A, vf &b, vf &x)
{
    for (int i = 0; i < n; i++)
    {
        x[i] = b[i];
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

void GS(int n, int k, vvf &A, vf &b, vf &x)
{
    for (int i = 0; i < n; i++)
        x[i] = 0;

    for (int i = 0; i < k; i++)
    {
        GS_itr(n, A, b, x);
    }
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

    GS(n, 10, a, b, x);

    for (int i = 0; i < n; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    return 0;
}
