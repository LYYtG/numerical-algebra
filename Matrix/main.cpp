/*
 * @Author: your name
 * @Date: 2021-04-10 13:49:25
 * @LastEditTime: 2021-04-13 13:39:02
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
    int n;
    cin >> n;
    vector<vector<double> >A(n),B(n);
    double t;
    int i,j,k;
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            cin >> t;
            A[i].push_back(t);
        }
    }
    cout << CondInf(A) << endl;
    return 0;
}