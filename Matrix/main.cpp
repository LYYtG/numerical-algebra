/*
 * @Author: your name
 * @Date: 2021-04-10 13:49:25
 * @LastEditTime: 2021-04-11 22:19:18
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \undefinedd:\numerical-algebra\Matrix\main.cpp
 */
#include <iostream>
#include "Matrix.h"
using namespace std;
#define MAXN 4
int main()
{
    vector<vector<double> >A(MAXN);
    double t;
    int i,j,k;
    /**
    for(i = 0;i<MAXN;i++){
        for(j = 0;j<MAXN;j++){
            cin >> t;
            A[i].push_back(t);
        }
    }
    */
    A = Hilbert(MAXN);
    cout << CondInf(A) << endl;
    return 0;
}