/**
 * @file Eig.h
 * @author 李杨野 (1300096763@qq.com 3190103519@zju.edu.cn)
 * @brief Functions for Eigenvalue
 * @version 0.1
 * @date 2021-05-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include "Matrix.h"
/**
 * @brief Construct the matrix according to the polynomial.
 * 
 * @param x 
 * @return vector<vector<double > > Matrix
 */
vector<vector<double > > Polynomial(vector<double> x)
{
    int n = x.size();
    vector<vector<double > > A(n,vector<double> (n));
    int i,j,k;
    for(i = 0;i<n-1;i++){
        A[i+1][i] = 1;
        A[i][n-1] = -x[i];
    }
    A[n-1][n-1] = -x[n-1];
    return A;
}
/**
 * @brief Subfunction of power method, to judge if the iteration is about to end.
 * 
 * @param u 
 * @param v 
 * @return true 
 * @return false 
 */
bool uv(vector<double> u,vector<double> v)
{
    int i;
    double t;
    int n = u.size();
    for(i = 0;i<n;i++){
        t = fabs(u[i]/v[i]);
        if(fabs(t-1)>=0.000001){
            return true;
        }
    }
    return false;
}
/**
 * @brief Power Method. Can only apply to specific Matrixs.
 * 
 * @param A 
 * @return vector<double> 
 */
vector<double> PowerMethod(vector<vector<double > > A)
{
    int n = A[0].size();
    vector<double> u(n),v(n),p(n);
    double t;
    int i,j,k = 0;
    for(i = 0;i<n;i++)
    {
        u[i] = 1.0;//Initialize...
    }
    do{
        v = u;//v = u_k-1        
        u = Multiply(A,u);
        t = InfiniteNorm(u);
        for(i = 0;i<n;i++){
            u[i]/=t;
        }
        //Print(u);
        k++;
    }while(uv(u,v)&&k<100);
    if(k>=100)cout <<"Too much iteration!"<<endl;
    cout << t << endl;
    return u;
}
/**
 * @brief Initial QR, with low time efficiency. NOT used in the homework.
 * 
 * @param A 
 */
void QREig(vector<vector<double> > A)
{
    int n = A[0].size();
    vector<vector<double> > Q(n,vector<double> (n));
    int i,j,k;
    for(k = 1;k<50;k++){
        //repeat 50 times… well that's not precise. For some real matrixs that gives a good result, though.
        Q = QR(A);
        A = Multiply(A,Q);
    }
    for(i = 0;i<n;i++)
    {
        cout << A[i][i] << " ";
    }
}