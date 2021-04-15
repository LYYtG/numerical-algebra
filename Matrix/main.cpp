/**
 * @file main.cpp
 * @author 李杨野 (1300096763@qq.com 3190103519@zju.edu.cn)
 * @brief test file.
 * @version 0.1
 * @date 2021-04-15
 * 
 * @copyright Copyright (c) 2021
 * 
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
    B = Inv(A);
    Print(B);
    return 0;
}