#include <iostream>
#include <vector>
#include "linear_eq.h"

using namespace std;


int main()
{
    int ne,nv;
    cin>>ne>>nv;
    linear_eq le(ne,nv);
    le.input();
    cout<<le;
    le.gauss_elemination(true);
    cout<<le;
    return 0;
}