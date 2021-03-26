/*
18MA20015 | Harshal Dupare | Lab 10 | 24-3-2021 |
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

struct assignment_prob
{
    ll n, m;

    vvr c;
    vvr main_c;
    vvr x;

    vector<bool> crossed_row;
    vector<bool> crossed_column;

    vector<bool> applicant_assinged;
    bool printtbl;
    ll problem_balanced; // +1 if jobs are more than applicants -1 if otherwise
    vector<pair<ll, ll>> marked;

    assignment_prob()
    {
        this->n = 0;
        this->m = 0;
        this->problem_balanced = 0;
        this->printtbl = true;
    }

    void input()
    {
        cout << "Number of Applicants: ";
        cin >> this->n;
        cout << "Number of Jobs: ";
        cin >> this->m;
        this->c.assign(this->n, vr(this->m));
        this->x.assign(this->n, vr(this->m));
        this->crossed_row.assign(this->n, false);
        this->crossed_column.assign(this->m, false);
        this->applicant_assinged.assign(this->m, false);

        if (m - n > 0)
        {
            this->problem_balanced = -1ll;
        }
        else if (n - m > 0)
        {
            this->problem_balanced = 1ll;
        }

        cout << "Enter cost of assigning i'th person to j'th job for next " << this->n << " lines in the format:\nc[i][0] c[i][1]... c[i][m-1] \n";
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                cin >> c[i][j];
                x[i][j] = (R)0;
            }
        }

        this->main_c = c;

        this->balance();
    }

    void init()
    {
        for (ll i = 0; i < this->n; i++)
        {
            R minv = infinity_R;
            for (ll j = 0; j < this->m; j++)
            {
                if (this->c[i][j] < minv)
                {
                    minv = this->c[i][j];
                }
            }

            for (ll j = 0; j < this->m; j++)
            {
                this->c[i][j] -= minv;
            }
        }

        for (ll i = 0; i < this->m; i++)
        {
            R minv = infinity_R;
            for (ll j = 0; j < this->n; j++)
            {
                if (this->c[j][i] < minv)
                {
                    minv = this->c[j][i];
                }
            }

            for (ll j = 0; j < this->m; j++)
            {
                this->c[j][i] -= minv;
            }
        }
    }

    string solution()
    {
        string s = "";
        for (ll i = 0; i < this->n; i++)
        {

            for (ll j = 0; j < this->m; j++)
            {

                if (this->x[i][j] == 1)
                    s += "Applicant " + to_string(i) + " Assigned to " + to_string(j) + " Job\n";
            }
        }

        return s;
    }

    R compute_cost()
    {
        R cost = 0;
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                cost += this->x[i][j] * this->main_c[i][j];
            }
        }
        return cost;
    }

    void balance()
    {
        if (this->n == this->m)
            return;

        if (this->n > this->m)
        {
            cout << "Adding Extra zero assignment cost Jobs\n";
            for (ll i = 0; i < this->n; i++)
            {
                this->x[i].resize(this->n, 0);
                this->c[i].resize(this->n, 0);
                this->main_c[i].resize(this->n, 0);
                this->crossed_column.resize(this->n, false);
                this->applicant_assinged.resize(this->n, false);
            }
            this->m = this->n;
        }
        else
        {
            cout << "Adding Extra zero assignment cost Applicant\n";
            this->x.resize(this->m, vr(this->m));
            this->c.resize(this->m, vr(this->m));
            this->main_c.resize(this->m, vr(this->m));
            this->crossed_row.resize(this->m, false);
            this->n = this->m;
        }
    }

    void reset()
    {
        this->x.assign(this->n, vr(this->m, 0));
    }

    bool covers_all()
    {
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                if (((this->crossed_column[j] || this->crossed_row[i]) == false) && (abs(this->c[i][j]) < _EPSILON))
                {
                    return false;
                }
            }
        }
        return true;
    }

    void update_crossed_bruteforce()
    {
        vector<ll> roco(this->n + this->m, 1);
        vector<bool> oco, oro;
        ll mxv = 0;
        bool found = false;

        for (ll max_1 = 1; max_1 < this->n + this->m; max_1++)
        {
            for (ll i = 0; i < max_1; i++)
            {
                roco[i] = 0;
            }
            do
            {
                for (ll i = 0; i < roco.size(); i++)
                {
                    if (i < this->n)
                    {
                        if (roco[i] == 1)
                            this->crossed_row[i] = true;
                        else
                            this->crossed_row[i] = false;
                    }

                    if (i >= this->n)
                    {
                        if (roco[i] == 1)
                            this->crossed_column[i - this->n] = true;
                        else
                            this->crossed_column[i - this->n] = false;
                    }
                }

                ll ctr = count(this->crossed_row.begin(), this->crossed_row.end(), true);
                ll ctc = count(this->crossed_column.begin(), this->crossed_column.end(), true);
                ll ct = ctr + ctc;
                if (ctr == 0 || ctc == 0)
                {
                    continue;
                }
                if (this->covers_all())
                {
                    found = true;
                    ll ctnz = (this->n - ctr) * (this->m - ctc);
                    if (ctnz > mxv)
                    {
                        oco = this->crossed_column;
                        oro = this->crossed_row;
                        mxv = ctnz;
                    }
                }

            } while (next_permutation(roco.begin(), roco.end()));
        }

        if (!found)
        {
            cout << "Coudnt find updated cover\n";
            return;
        }

        this->crossed_column = oco;
        this->crossed_row = oro;

        if (this->printtbl)
        {
            cout << "Crossed Rows:";
            for (ll i = 0; i < this->n; i++)
            {
                if (this->crossed_row[i])
                    cout << i << ", ";
            }
            cout << "\nCrossed Columns:";
            for (ll i = 0; i < this->m; i++)
            {
                if (this->crossed_column[i])
                    cout << i << ", ";
            }
            cout << endl;
            cout << "Table is:\n";
            cout << (*this);
            cout << endl;
        }

        // reducing the non covered cells
        this->reduce();
    }

    void update_crossed_scanning()
    {
        this->crossed_column.assign(this->m, false);
        this->crossed_row.assign(this->n, false);
        this->marked.clear();

        // row scanning
        for (ll i = 0; i < this->n; i++)
        {
            ll zc = 0;
            ll zid = -1;
            for (ll j = 0; j < this->m; j++)
            {
                if (this->crossed_column[j] || this->crossed_row[i])
                {
                    continue;
                }

                if (abs(this->c[i][j]) <= _EPSILON)
                {
                    zc++;
                    zid = j;
                }
            }
            if (zc != 1)
            {
                continue;
            }

            this->crossed_column[zid] = true;
            this->marked.push_back({i, zid});
        }

        // column scanning
        for (ll i = 0; i < this->m; i++)
        {
            ll zc = 0;
            ll zid = -1;
            for (ll j = 0; j < this->n; j++)
            {
                if (this->crossed_column[i] || this->crossed_row[j])
                {
                    continue;
                }

                if (abs(this->c[j][i]) <= _EPSILON)
                {
                    zc++;
                    zid = j;
                }
            }
            if (zc != 1)
            {
                continue;
            }

            this->crossed_row[zid] = true;
            this->marked.push_back({zid, i});
        }

        if (this->printtbl)
        {
            cout << "Selected Zeros are at location (row_id, column_id):\n";
            for (auto mp : marked)
            {
                cout << "(" << mp.first << ", " << mp.second << "), ";
            }
            cout << endl;

            cout << "Crossed Rows:";
            for (ll i = 0; i < this->n; i++)
            {
                if (this->crossed_row[i])
                    cout << i << ", ";
            }
            cout << "\nCrossed Column:";
            for (ll i = 0; i < this->m; i++)
            {
                if (this->crossed_column[i])
                    cout << i << ", ";
            }
            cout << endl;
            cout << "Table is:\n";
            cout << (*this);
            cout << endl;
        }

        if (!this->covers_all())
        {
            cout << "Doesn't cover all Zeros\n";
            return;
        }

        if (marked.size() == this->n)
        {
            cout << "Optimality has reached\n";
            return;
        }
        else
        {
            cout << "Optimality hasn't yet reached\n"
                 << endl;
            this->reduce();
        }
    }

    pair<ll, ll> smallest_uncovered()
    {
        pair<ll, ll> p;
        p.first = -1;
        p.second = -1;
        R minv = infinity_R;
        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                if ((this->crossed_column[j] || this->crossed_row[i]) == false && (this->c[i][j] < minv))
                {
                    minv = this->c[i][j];
                    p.first = i;
                    p.second = j;
                }
            }
        }

        return p;
    }

    bool find_feasible(ll row_i)
    {
        if (row_i >= this->n)
        {
            return true;
        }

        for (ll i = 0; i < this->m; i++)
        {
            if (abs(this->c[row_i][i]) <= _EPSILON && (!this->applicant_assinged[i]))
            {
                this->x[row_i][i] = 1;
                this->applicant_assinged[i] = true;
                bool found = this->find_feasible(row_i + 1);
                if (found)
                {
                    return true;
                }
                this->x[row_i][i] = 0;
                this->applicant_assinged[i] = false;
            }
        }

        return false;
    }

    void reduce()
    {
        auto lowp = this->smallest_uncovered();
        if (lowp.first == -1 || lowp.second == -1)
        {
            cout << "No smallest element found\n";
        }

        R lv = this->c[lowp.first][lowp.second];
        cout << "Smallest uncovered element found at: (" << lowp.first << ", " << lowp.second << ") with value of " << lv << endl;

        for (ll i = 0; i < this->n; i++)
        {
            for (ll j = 0; j < this->m; j++)
            {
                if ((this->crossed_column[j] || this->crossed_row[i]) == false)
                {
                    this->c[i][j] -= lv;
                }

                if ((this->crossed_column[j] && this->crossed_row[i]) == true)
                {
                    this->c[i][j] += lv;
                }
            }
        }
    }

    void hungarian(ll use_brute_force = 0,ll max_itr = 20)
    {
        this->init();
        if (this->printtbl)
        {
            cout << "\n\nAfter reducing the inital rows and columns we get:\n";
            cout << (*this);
        }

        ll i = 1;
        while (max_itr > i)
        {
            this->reset();
            bool found = this->find_feasible(0);

            if (found)
            {
                cout << "Optimal solution found\n";
                break;
            }

            if (this->printtbl)
            {
                cout << "Iteration number " << i << "\n";
            }

            if (use_brute_force)
            {
                this->update_crossed_bruteforce();
            }
            else
            {
                this->update_crossed_scanning();
            }

            i++;
        }

        cout << "Iterations terminated\n";
        this->crossed_column.assign(this->m, false);
        this->crossed_row.assign(this->n, false);
    }

    friend ostream &operator<<(ostream &os, assignment_prob &tp)
    {
        vector<ll> applicants(tp.n, 0);
        vector<ll> jobs(tp.m, 0);
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
        os << "     Applicants ";
        os << endl;
        for (ll i = 0; i < tp.n; i++)
        {
            os.width(wd);
            s = to_string(i) + " | ";
            os << s;
            for (ll j = 0; j < tp.m; j++)
            {
                os.width(wd);
                if (tp.crossed_column[j] || tp.crossed_row[i])
                {
                    s = "* ";
                }
                else
                {
                    s.clear();
                }
                s += to_string((ll)tp.x[i][j]) + " : " + to_string(tp.c[i][j]) + " | ";
                os << s;
                if (tp.x[i][j] == 1)
                {
                    applicants[i] = 1;
                    jobs[j] = 1;
                }
            }

            os.width(wd);
            s = to_string(applicants[i]) + "/" + to_string(1);
            os << s << endl;
        }

        os.width(wd);
        os << "Jobs | ";
        for (ll j = 0; j < tp.m; j++)
        {
            os.width(wd);
            os.width(wd);
            s = to_string((ll)jobs[j]) + "/" + to_string(1) + " | ";
            os << s;
        }
        os << endl;

        return os;
    }
};

void lab_10()
{
    assignment_prob tp;

    tp.input();
    cout << "Initial Table is:\n";
    cout << tp;

    cout << "Input what method to use for calculating min number of crossing line:\n"
         << " \"0\" for Greedy method ( FAST but sometimes doesn't give exact solution)\n \"1\" for Brute Force method ( SLOW but gives exact solution)\n ";
    ll m_type = 0;
    cin >> m_type;
    tp.hungarian(m_type);
    cout << tp;

    cout << "Final Optimal cost is :" << tp.compute_cost() << endl;
    cout << "Solution is:" << endl;
    cout << tp.solution();
}

int main()
{
    lab_10();
    return 0;
}
