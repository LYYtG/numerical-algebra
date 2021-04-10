/*
 * @Author: your name
 * @Date: 2021-04-10 13:49:25
 * @LastEditTime: 2021-04-10 13:55:55
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \undefinedd:\numerical-algebra\Matrix\main.cpp
 */
#include <iostream>
#include "Matrix.h"
using namespace std;
#define MAXN 3
int main()
{
    vector<vector<double> >A(MAXN),L(MAXN),U(MAXN);
    vector<double> b(MAXN);
    double t;
    int i,j,k;
    for(i = 0;i<MAXN;i++){
        for(j = 0;j<MAXN;j++){
            cin >> t;
            A[i].push_back(t);
        }
    }
    LU(A,L,U);
    Print(L);
    Print(U);
    return 0;
}