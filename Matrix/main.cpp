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
#include "Eig.h"
using namespace std;
#define MAXN 4
int main()
{
    int n;
    cin >> n;
    vector<vector<double> > A(n, vector<double> (n));
    int i,j;
    for(i = 0;i<n;i++){
        if(i == 0){
            A[i][i] = 2;
            A[i][i+1] = -1;
        }else if(i == n-1){
            A[i][i] = 2;
            A[i][i-1] = -1;
        }else{
            A[i][i+1] = -1;
            A[i][i] = 2;
            A[i][i-1] = -1;
        }
    }
    Print(A);
    double ret = Dich(A,100);
    cout << ret<< endl;
    vector<double> v(n);
    v = Rev_PowerMethod(A,ret);
    Print(v);
    return 0;
}