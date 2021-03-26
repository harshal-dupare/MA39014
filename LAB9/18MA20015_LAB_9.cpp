/*
18MA20015 | Harshal Dupare | Lab 9 | 17-3-2021 |
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
typedef double R;
typedef long long ll;
typedef vector<vector<R>> vvr;
typedef vector<R> vr;
ll itrcount = 0;
ll wd = 25;
R _EPSILON = 1e-4;
R infinity_R = 1e15;

struct transport_prob
{
    ll n, m;
    vr supply;
    vr demand;

    vr alloc_supply;
    vr alloc_demand;

    vvr c;
    vvr x;
    R net_supply;
    R net_demand;

    ll problem_balanced; // +1 if supply is more than demand -1 if otherwise

    transport_prob()
    {
        this->n = 0;
        this->m = 0;
        this->problem_balanced = 0;
    }

    void input()
    {
        cout << "Number of sources: ";
        cin >> this->n;
        cout << "Number of destinations: ";
        cin >> this->m;
        supply.resize(this->n, (R)0);
        demand.resize(this->m, (R)0);
        alloc_supply.resize(this->n, (R)0);
        alloc_demand.resize(this->m, (R)0);
        c.assign(this->n, vr(this->m));
        x.assign(this->n, vr(this->m));
        cout << "Enter the supply at the i'th source for next " << this->n << " lines:\n ";
        for (ll i = 0; i < this->n; i++)
        {
            cin >> this->supply[i];
            net_supply += this->supply[i];
        }
        cout << "Enter the demand at the i'th destination for next " << this->m << " lines:\n ";
        for (ll i = 0; i < this->m; i++)
        {
            cin >> this->demand[i];
            net_demand += this->demand[i];
        }

        if (net_demand - net_supply > _EPSILON)
        {
            this->problem_balanced = -1ll;
        }
        else if (net_supply - net_demand > _EPSILON)
        {
            this->problem_balanced = 1ll;
        }

        cout << "Enter cost shipping of commodity for i'th source for next " << this->n << " lines in the format:\nc[i][0] c[i][1]... c[i][m-1] \n";
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                cin >> c[i][j];
                x[i][j] = (R)0;
            }
        }
    }

    void northwest_corner()
    {
        ll cell_i = 0, cell_j = 0;
        while ((cell_i + 1ll) <= this->n && (cell_j + 1ll) <= this->m)
        {
            ll ncell_i = cell_i, ncell_j = cell_j;
            R avl_supp = max(this->supply[cell_i] - this->alloc_supply[cell_i], (R)0);
            R req_dem = max(this->demand[cell_j] - this->alloc_demand[cell_j], (R)0);
            R allotment = req_dem;
            if (abs(avl_supp - req_dem) <= _EPSILON)
            {
                if (cell_i + 2ll == this->n)
                {
                    ncell_i++;
                }
                else
                {
                    ncell_j++;
                }
            }
            else if (avl_supp - req_dem > 0)
            {
                ncell_j++;
                allotment = req_dem;
            }
            else
            {
                ncell_i++;
                allotment = avl_supp;
            }

            this->alloc_supply[cell_i] += allotment;
            this->alloc_demand[cell_j] += allotment;
            this->x[cell_i][cell_j] += allotment;

            cell_i = ncell_i;
            cell_j = ncell_j;
        }
    }

    void lowestcost_cell()
    {
        vector<pair<ll, ll>> ps;
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                ps.push_back({i, j});
            }
        }

        sort(ps.begin(), ps.end(), [&](pair<ll, ll> &p, pair<ll, ll> &q) { return this->c[p.first][p.second] < this->c[q.first][q.second]; });
        vector<bool> rdone(this->n, false);
        vector<bool> cdone(this->m, false);
        ll id = 0;
        while (id < ps.size())
        {
            ll cell_i = ps[id].first, cell_j = ps[id].second;
            if (rdone[cell_i] || cdone[cell_j])
            {
                id++;
                continue;
            }
            R avl_supp = max(this->supply[cell_i] - this->alloc_supply[cell_i], (R)0);
            R req_dem = max(this->demand[cell_j] - this->alloc_demand[cell_j], (R)0);
            R allotment;
            if (abs(avl_supp - req_dem) <= _EPSILON)
            {
                allotment = req_dem;
                cdone[cell_j] = true;
            }
            else if (avl_supp - req_dem > 0)
            {
                allotment = req_dem;
                cdone[cell_j] = true;
            }
            else
            {
                allotment = avl_supp;
                rdone[cell_i] = true;
            }

            id++;

            this->alloc_supply[cell_i] += allotment;
            this->alloc_demand[cell_j] += allotment;
            this->x[cell_i][cell_j] += allotment;
        }
    }

    R compute_cost()
    {
        R cost = 0;
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                cost += this->x[i][j] * this->c[i][j];
            }
        }

        return cost;
    }

    void solve_dual(vector<pair<ll, ll>> &vps, vr &u, vr &v, vector<bool> &uf, vector<bool> &vf)
    {
        ll found_cout = 1;
        while (found_cout < u.size() + v.size())
        {
            for (auto p : vps)
            {
                if (uf[p.first] && (!vf[p.second]))
                {
                    vf[p.second] = true;
                    v[p.second] = this->c[p.first][p.second] - u[p.first];
                    found_cout++;
                }
                else if ((!uf[p.first]) && (vf[p.second]))
                {
                    uf[p.first] = true;
                    u[p.first] = this->c[p.first][p.second] - v[p.second];
                    found_cout++;
                }
            }
        }
    }

    // 0 for row 1 for column search
    bool find_path(vector<pair<ll, ll>> &path, ll roc, vector<vector<bool>> &onpath)
    {
        ll i = path[path.size() - 1].first, j = path[path.size() - 1].second;

        if (roc == 0)
        {
            for (ll t = 0; t < this->m; t++)
            {
                if (path.size() > 2 && i == path[0].first && t == path[0].second)
                {
                    return true;
                }

                if (t == j)
                {
                    continue;
                }

                if (this->x[i][t] > (R)0 && (!onpath[i][t]))
                {
                    path.push_back({i, t});
                    onpath[i][t] = true;
                    bool nxt_found = this->find_path(path, 1 - roc, onpath);

                    if (nxt_found)
                    {
                        return true;
                    }
                    else
                    {
                        path.pop_back();
                        onpath[i][t] = false;
                    }
                }
            }
        }
        else
        {

            for (ll t = 0; t < this->n; t++)
            {
                if (path.size() > 2 && t == path[0].first && j == path[0].second)
                {
                    return true;
                }

                if (t == i)
                {
                    continue;
                }

                if (this->x[t][j] > (R)0 && (!onpath[t][j]))
                {
                    path.push_back({t, j});
                    onpath[t][j] = true;
                    bool nxt_found = this->find_path(path, 1 - roc, onpath);

                    if (nxt_found)
                    {
                        return true;
                    }
                    else
                    {
                        path.pop_back();
                        onpath[t][j] = false;
                    }
                }
            }
        }

        return false;
    }

    void modi_method()
    {
        ll num_basic = 0;
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                if (this->x[i][j] > 0.0)
                {
                    num_basic++;
                }
            }
        }
        if (num_basic != (m + n - 1))
        {
            cout << "Method cannot be applied with this initial solution as not enought basic variables\n";
            return;
        }

        vr u(this->n, 0);
        vector<bool> uf(this->n, false);
        uf[0] = true;
        vr v(this->m, 0);
        vector<bool> vf(this->n, false);
        ll itr = 0;
        while (true)
        {
            cout << "||  Iteration number " << itr + 1 << "  ||\nn";
            vector<pair<ll, ll>> vps;
            vr u(this->n, 0);
            vector<bool> uf(this->n, false);
            uf[0] = true;
            vr v(this->m, 0);
            vector<bool> vf(this->n, false);

            for (ll i = 0; i < this->n; i++)
            {
                for (ll j = 0; j < this->m; j++)
                {
                    if (this->x[i][j] > 0.0)
                    {
                        vps.push_back({i, j});
                    }
                }
            }

            this->solve_dual(vps, u, v, uf, vf);
            R maxv_neg = 0;
            ll cell_i = -1, cell_j = -1;

            for (ll i = 0; i < this->n; i++)
            {
                for (ll j = 0; j < this->m; j++)
                {
                    if (this->c[i][j] - u[i] - v[j] < maxv_neg)
                    {
                        cell_i = i;
                        cell_j = j;
                        maxv_neg = this->c[i][j] - u[i] - v[j];
                    }
                }
            }
            cout << "Solution of Dual problem is:\n";
            for (ll kt = 0; kt < this->n; kt++)
            {
                cout << "u_" << kt + 1 << " = " << u[kt] << ", ";
            }
            cout << endl;
            for (ll kt = 0; kt < this->m; kt++)
            {
                cout << "v_" << kt + 1 << " = " << v[kt] << ", ";
            }
            cout << endl;

            if (cell_i == -1)
            {
                cout << "No negative Cij - ui - vj found, hence itteration completed\n";
                break;
            }

            cout << "Cell considered in this itteration is : [" << cell_i << ", " << cell_j << "]\n\n";

            vector<pair<ll, ll>> path;
            path.push_back({cell_i, cell_j});

            vector<vector<bool>> onpath(this->n, vector<bool>(this->m));
            onpath[cell_i][cell_j] = true;

            bool found_path = this->find_path(path, 0, onpath);

            if (!found_path)
            {
                cout << "Path not found\n";
                return;
            }

            R min_shift = infinity_R;
            for (ll i = 1; i < path.size(); i += 2)
            {
                min_shift = min(min_shift, this->x[path[i].first][path[i].second]);
            }

            ll sg = 1;
            for (auto p : path)
            {
                this->x[p.first][p.second] += sg * min_shift;
                sg *= (-1);
            }

            cout << "After " << itr + 1 << " iteration table is:\n\n";
            cout << (*this);
        }

        return;
    }

    friend ostream &operator<<(ostream &os, transport_prob &tp)
    {
        os << endl;
        string s;
        os.width(wd);
        os << "      | ";
        for (ll j = 0; j < tp.m; j++)
        {
            os.width(wd);
            s = to_string(j) + " | ";
            os << s;
        }
        os.width(wd);
        os << "     Supply ";
        os << endl;
        for (ll i = 0; i < tp.n; i++)
        {
            os.width(wd);
            s = to_string(i) + " | ";
            os << s;
            for (ll j = 0; j < tp.m; j++)
            {
                os.width(wd);
                s = to_string(tp.x[i][j]) + " : " + to_string(tp.c[i][j]) + " | ";
                os << s;
            }

            os.width(wd);
            s = to_string(tp.alloc_supply[i]) + "/" + to_string(tp.supply[i]);
            os << s << endl;
        }

        os.width(wd);
        os << "Demand | ";
        for (ll j = 0; j < tp.m; j++)
        {
            os.width(wd);
            os.width(wd);
            s = to_string(tp.alloc_demand[j]) + "/" + to_string(tp.demand[j]) + " | ";
            os << s;
        }
        os << endl;

        return os;
    }
};

void lab_9()
{
    transport_prob tp;

    tp.input();
    cout << tp;

    // tp.northwest_corner();
    tp.lowestcost_cell();
    cout << "After finding Basic Feasible solution we get\n";
    cout << tp;

    tp.modi_method();
    cout << "Final table is\n";
    cout << tp;

    cout << "Final Optimal cost is :" << tp.compute_cost() << endl;
}

int main()
{
    lab_9();
    return 0;
}
