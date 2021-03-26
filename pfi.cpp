#include <iostream>
#include <vector>
#include "library/debug.h"

using namespace std;
typedef double R;
typedef long long ll;
typedef vector<vector<R>> vvr;
typedef vector<R> vr;

class pfi
{
public:
    ll n;
    vvr I;
    vector<vvr> E;
    vvr B;
    vvr inv_B;
    pfi(ll n)
    {
        this->n = n;
        this->I = vvr(n, vr(n));
        for (ll i = 0; i < n; i++)
        {
            this->I[i][i] = 1;
        }

        this->B = this->I;

        this->E = vector<vvr>(n, this->I);

        this->inv_B = this->I;
    }

    vvr mmult(vvr X, vvr Y)
    {
        vvr XY(X.size(), vr(Y[0].size()));
        for (ll i = 0; i < X.size(); i++)
        {
            for (ll j = 0; j < Y[0].size(); j++)
            {
                for (ll k = 0; k < X[0].size(); k++)
                {
                    XY[i][j] += X[i][k] * Y[k][j];
                }
            }
        }

        return XY;
    }

    void input()
    {
        for (ll i = 0; i < n; i++)
            for (ll j = 0; j < n; j++)
            {
                cin >> B[i][j];
            }
    }

    void print(vvr X)
    {
        for (ll i = 0; i < X.size(); i++)
        {
            for (ll j = 0; j < X[0].size(); j++)
            {
                cout << X[i][j] << " ";
            }
            cout << endl;
        }

        cout << endl;
    }

    void calculate()
    {
        vvr eta_i(n, vr(1));
        for (ll i = 0; i < n; i++)
        {
            debug(i);
            vvr ci(n, vr(1));
            for (ll j = 0; j < n; j++)
            {
                ci[j][0] = B[j][i];
            }
            debug(ci);

            eta_i = this->mmult(inv_B, ci);
            debug(eta_i);

            if (eta_i[i][0] == 0)
            {
                cout << "Division by Zero\n";
                return;
            }
            eta_i[i][0] = (((R)1.0) / eta_i[i][0]);
            for (ll j = 0; j < n; j++)
            {
                if (i != j)
                    eta_i[j][0] =(R)( (-1.0) * (eta_i[j][0] * eta_i[i][0]));

                this->E[i][j][i] = eta_i[j][0];
            }
            debug(E[i]);

            this->inv_B = this->mmult(this->E[i], this->inv_B);
            debug(inv_B);
        }
    }
};

void test()
{
    pfi t(3ll);

    t.input();
    debug(t.B);
    
    t.calculate();

    debug(t.inv_B);
}

int main()
{
    test();
    return 0;
}

/*
2 4 0
0 2 0
0 6 5
*/