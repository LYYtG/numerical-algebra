/*
 * @Author: 李杨野
 * @Date: 2021-04-07 17:02:05
 * @LastEditTime: 2021-04-13 13:54:05
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \undefinedd:\szds\Hilbert\Matrix.h
 */
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
using namespace std;
vector<double> UpperMatrix(vector<vector<double> > &matrix,vector<double> &b)
{
    int n = matrix[0].size();
    int i,j,k;
    double t;
    vector<double> x(n);
    //回代法
    x[n-1] = b[n-1]/matrix[n-1][n-1];
    for(i = n-2;i>=0;i--){
        t = b[i];
        for(j = i+1;j<n;j++){
            t-=matrix[i][j]*x[j];
        }
        x[i] = t/matrix[i][i];
    }
    return x;
}
vector<double> LowerMatrix(vector<vector<double> > &matrix,vector<double> &b)
{
    int n = matrix[0].size();
    int i,j,k;
    double t;
    vector<double> x(n);
    //前代法
    x[0] = b[0]/matrix[0][0];
    for(i = 1;i<n;i++){
        t = b[i];
        for(j = 0;j<i;j++){
            t-=matrix[i][j]*x[j];
        }
        x[i] = t/matrix[i][i];
    }
    return x;
}
void Transposition(vector<vector<double> > &matrix)
{
    //矩阵原位转置
    int n = matrix[0].size();
    int i,j,k;
    double t;
    for(i = 0;i<n;i++){
        for(j = i;j<n;j++){
            
            t = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = t;
        }
    }
}
void TranspositionOutput(vector<vector<double> > &matrix,vector<vector<double> > &matrixT)
{
    //矩阵转置
    int n = matrix[0].size();
    int i,j,k;
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            matrixT[i].push_back(matrix[j][i]);
        }
    }
}
void ColumnPivot(vector<vector<double> > &matrix,vector<double> &b)
{
    //效果是把矩阵通过Gauss变换转换为上三角阵。
    int n = matrix[0].size();
    int i,j,k;
    int max;
    double t;
    for(k = 0;k<n;k++){
        max = k;
        for(i = k+1;i<n;i++){
            if(fabs(matrix[i][k])>fabs(matrix[max][k]))
            {
                max = i;
            }
        }            
        for(j = 0;j<n;j++)
        {
            t = matrix[max][j];
            matrix[max][j] = matrix[k][j];
            matrix[k][j] = t;
        }           
        t = b[max];
        b[max] = b[k];
        b[k] = t;
        //此时kk个元素应该为k列最大者，不然，本列全部为0，直接进行下一轮循环。
        //输出三角矩阵。
        if(matrix[k][k]!=0){
            for(i = k+1;i<n;i++){
                t = matrix[i][k]/matrix[k][k];
                for(j = k;j<n;j++){
                    matrix[i][j]-=t*matrix[k][j];
                }
                b[i]-=t*b[k];
            }
        }     
    }
}
void Print(vector<vector<double> > &matrix,vector<double> &b)
{
    //同时打印系数矩阵和b。
    int i,j;
    int n = matrix[0].size();
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            cout << matrix[i][j] << " ";
        }
        cout << b[i] << endl;
    }
}
void Print(vector<vector<double> > &matrix)
{
    //同时打印系数矩阵和b。
    int i,j;
    int n = matrix[0].size();
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
void Print(vector<double> &b)
{
    int n = b.size();
    for(int j = 0;j<n;j++){
            cout << b[j] << " ";
        }
        cout << endl;
}
vector<double> GaussColumnPivot(vector<vector<double> > matrix,vector<double> b)
{
    //列主元高斯消元法并输出解x。
    
    int n = matrix[0].size();

    ColumnPivot(matrix,b);
    vector<double> x = UpperMatrix(matrix,b);
    return x;
}
vector<vector<double> > Hilbert(int n)
{
    //生成n阶Hilbert矩阵
    vector<vector<double> > matrix(n);
    double i,j;
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            matrix[i].push_back(1.0/(i+j+1));
        } 
    }
    return matrix;
}
void Reverse(vector<vector<double> > &matrix)
{
    ;//to be done.
}
double OneNorm(vector<double> &b)
{
    //向量的1-范数
    int n = b.size();
    double sum = 0;
    for(int i = 0;i<n;i++){
        sum+=fabs(b[i]);
    }
    return sum;
}
double OneNorm(vector<vector<double> > &matrix)
{
    //矩阵的1-范数
    int n = matrix[0].size();
    int i,j,k;
    double max,sum = 0.0;
    for(i = 0;i<n;i++){
        sum+=fabs(matrix[i][0]);
    }
    max = sum;
    for(j = 1;j<n;j++){
        sum = 0.0;
        for(i = 0;i<n;i++){
            sum+=fabs(matrix[i][j]);
        }
        if(sum > max)max = sum;
    }
    return max;
}
double InfiniteNorm(vector<double> &b)
{
    //向量的∞范数
    int n = b.size();
    double max = fabs(b[0]);
    for(int i = 1;i<n;i++){
        if(fabs(b[i]>max))
            max = fabs(b[i]);
    }
    return max;
}
vector<double> ej(vector<double> &b)
{
    int n = b.size();
    double max = fabs(b[0]);
    int t = 0;
    for(int i = 1;i<n;i++){
        if(fabs(b[i]>max))
        {
            max = fabs(b[i]);
            t = i;
        }
    }
    vector<double> x(n);
    for(int i = 0;i<n;i++){
        if(i == t){
            x[i] = 1;
        }
        else{
            x[i] = 0;
        }
    }
    return x;
}
vector<double> ej(int n,int j)
{
    vector<double> x(n);
    for(int i = 0;i<n;i++){
        if(i!=j){
            x[i] = 0;
        }else{
            x[i] = 1;
        }
    }
    return x;
}
double InfiniteNorm(vector<vector<double> > &matrix)
{
    //矩阵的∞范数
    int n = matrix[0].size();
    int i,j,k;
    double max,sum = 0;
    for(i = 0;i<n;i++){
        sum+=fabs(matrix[0][i]);
    }
    max = sum;
    for(j = 1;j<n;j++){
        sum = 0;
        for(i = 0;i<n;i++){
            sum+=fabs(matrix[j][i]);
        }
        if(sum > max)max = sum;
    }
    return max;
}
int Sign(double x)
{
    //sgn(x)
    int ret;
    if(x>0){
        ret = 1;
    }else if(x == 0){
        ret = 0;
    }else if(x < 0){
        ret = -1;
    }
    return ret;
}
vector<double> Multiply(vector<vector<double> > &matrix,vector<double> &b)
{
    //矩阵和向量乘法
    int n = matrix[0].size();
    if(b.size()!=n){
        std::cerr << "Dimension no Match!";
    }
    vector<double > x(n);
    int i,j;
    for(i = 0;i<n;i++){
        x[i] = 0;
        for(j = 0;j<n;j++){
            x[i]+=matrix[i][j]*b[j];
        }
    }
    return x;
}
double Multiply(vector<double> &a,vector<double> &b)
{
    double sum = 0;
    int n = a.size();
    for(int i = 0;i<n;i++){
        sum+=a[i]*b[i];
    }
    return sum;
}
vector<vector<double> > CertainMatrix(int n)
{
    //生成n阶题目所示矩阵
    vector<vector<double> > matrix(n);
    double i,j;
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            if(i == j||j == n-1){
                matrix[i][j] = 1;
            }
            else if(i>j){
                matrix[i][j] =  0;
            }else{
                matrix[i][j] = -1;
            }
        } 
    }
    return matrix;
}
void LU(vector<vector<double> >&matrix,vector<vector<double> >&L,vector<vector<double> >&U)
{
    int n = matrix[0].size();
    int i,j,k;
    double t;
    L = matrix;
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            U[i].push_back(0);
        }
    }
    for(k = 0;k<n;k++){
        if(L[k][k] == 0){
            for(i = 0;i<n;i++){
                if(L[i][k]!=0)break;
            }
            if(i == n){
                //标志这一列全部是0，直接进行下一列的gauss变换。 
                continue;
            }else{
                //对非0行进行交换。
                for(j = 0;j<n;j++){
                    t = L[i][j];
                    L[i][j] = L[k][j];
                    L[k][j] = t;
                }
            }
        }
        //此时k列所有元素应该不为0。
        //用k+1:n区域存储矩阵L。
        for(i = k+1;i<n;i++){
            L[i][k] /= L[k][k];
        }
        //用k+1:n,k+1:n区域存储矩阵U。
        for(i = k+1;i<n;i++){
            for(j = k+1;j<n;j++){
                L[i][j]-=L[i][k]*L[k][j];
            }
        }
    }
    for(i = 0;i<n;i++){
        for(j = i;j<n;j++){
            U[i][j] = L[i][j];
            if(j == i){
                L[i][j] = 1;
            }else{
                L[i][j] = 0;
            }
        }
    }
}
double CondInf(vector<vector<double> > &A)
{
    //计算某些病态矩阵（Hilbert）不太精确。其他的尚可。
    int n = A[0].size();
    vector<vector<double> > AT(n);
    TranspositionOutput(A,AT);
    vector<double> x(n),w(n),v(n),z(n);
    double t;
    int i,j,k;
        for(j = 0;j<n;j++){
            x[j] = 1.0/n; 
        }
        k = 1;
        while(k == 1){
            w = GaussColumnPivot(AT,x);
            for(j = 0;j<n;j++){
                v[j] = Sign(w[j]);
            }
            z = GaussColumnPivot(A,v);
            if(InfiniteNorm(z)<=Multiply(z,x)){
                t = OneNorm(w);//A-1的无穷范数估计
                k = 0;
            }else{
                x = ej(z);
                k = 1;
            }
        }
        return t*InfiniteNorm(A);
}
void Inv(vector<vector<double> > &matrix,vector<vector<double> > &invmatrix)
{
    //对于病态矩阵求无穷范数条件数相比CondInf较为精确，仅作测试使用！计算量较大
    int n = matrix[0].size();
    int i,j,k;
    for(i = 0;i<n;i++){
        invmatrix[i] = GaussColumnPivot(matrix,ej(n,i));
    }
    return Transposition(invmatrix);
}