/**
 * @file Matrix.h
 * @author 李杨野 (1300096763@qq.com 3190103519@zju.edu.cn)
 * @brief Header file for matrix, for test only. Completely written by myself.
 * @version 0.1
 * @date 2021-04-15
 * 
 * @copyright Copyright (c) 2021 
 * 
 */
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
using namespace std;
/**
Class Matrix(){
    private:
        vector<vector<double> > matrix;
        int n;
        int m;
    public:
        Matrix();
        ~Matrix();
};
*/

/**I was intended to construct a class to simplify my code before. but it would take too much time to rewrite it.
 * so the class was abandoned.
 */

/**
 * @brief Solution of Upper Matrix equations.
 * 
 * @param matrix 
 * @param b 
 * @return vector<double> the solution.
 */
vector<double> UpperMatrix(vector<vector<double> > &matrix,vector<double> &b)
{
    int n = matrix[0].size();
    int i,j,k;
    double t;
    vector<double> x(n);
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
/**
 * @brief Solution of Lower Matrix equations
 * 
 * @param matrix 
 * @param b 
 * @return vector<double> the solution.
 */
vector<double> LowerMatrix(vector<vector<double> > &matrix,vector<double> &b)
{
    int n = matrix[0].size();
    int i,j,k;
    double t;
    vector<double> x(n);
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
/**
 * @brief in-place matrix Transposition.
 * 
 * @param matrix 
 */
void Transposition_In_Place(vector<vector<double> > &matrix)
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
/**
 * @brief not in-place matrix Transpostion
 * 
 * @param matrix 
 * @param matrixT where the output lies.
 */
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
/**
 * @brief another Transposition.
 * 
 * @param matrix 
 * @return vector<vector<double> > 
 */
vector<vector<double> > Transposition(vector<vector<double> > &matrix)
{
    int m = matrix.size();
    int n = matrix[0].size();
    vector<vector<double> > matrixT(n,vector<double> (m));
    int i,j;
    for(i = 0;i<n;i++){
        for(j = 0;j<m;j++){
            matrixT[i][j] = matrix[j][i];
        }
    }
    return matrixT;
}
/**
 * @brief Common Gauss column pivot decomposition. Subfunction of the GaussColumnPivot function below.
 * 
 * @param matrix 
 * @param b 
 */
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
/**
 * @brief Print the Matrix and vector. Useless function, I will delete this soon.
 * 
 * @param matrix 
 * @param b 
 */
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
/**
 * @brief Print the matrix.
 * 
 * @param matrix 
 */
void Print(vector<vector<double> > &matrix)
{
    int i,j;
    int m = matrix.size();
    int n = matrix[0].size();
    for(i = 0;i<m;i++){
        for(j = 0;j<n;j++){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
/**
 * @brief Print the vector.
 * 
 * @param b 
 */
void Print(vector<double> &b)
{
    int n = b.size();
    for(int j = 0;j<n;j++){
            cout << b[j] << " ";
        }
        cout << endl;
}
/**
 * @brief Gauss Column Pivot way to solve equations. It's worth mention that the matrix must be non-singular.
 * 
 * @param matrix 
 * @param b 
 * @return vector<double> the soluton vector.
 */
vector<double> GaussColumnPivot(vector<vector<double> > matrix,vector<double> b)
{
    //列主元高斯消元法并输出解x。
    
    int n = matrix[0].size();

    ColumnPivot(matrix,b);
    vector<double> x = UpperMatrix(matrix,b);
    return x;
}
/**
 * @brief create a Hilbert Matrix.
 * 
 * @param n order of matrix.
 * @return vector<vector<double> > n-ordered Hilbert Matrix.
 */
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
/**
 * @brief 1-norm of a vector.
 * 
 * @param b 
 * @return double 1-norm.
 */
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
/**
 * @brief 1-norm of a matrix, defined by 1-norm of vector.
 * 
 * @param matrix 
 * @return double 1-norm
 */
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
/**
 * @brief inf-norm of a matrix.
 * 
 * @param b 
 * @return double int-norm.
 */
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
/**
 * @brief according to the input vector generate a Ej. Subfunction of CondInf.
 * 
 * @param b 
 * @return vector<double> Ej
 */
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
/**
 * @brief common definition of ej, aka index vector.
 * 
 * @param n the order of vector.
 * @param j 
 * @return vector<double> Ej
 */
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
/**
 * @brief the inf-norm of a matrix.
 * 
 * @param matrix 
 * @return double inf norm.
 */
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
/**
 * @brief sgn(x) in math.
 * 
 * @param x 
 * @return int sgn(x)
 */
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
/**
 * @brief A*b, A belongs to Rn×n, b belongs to Rn. 
 * 
 * @param matrix 
 * @param b 
 * @return vector<double> A*b
 */
vector<double> Multiply(vector<vector<double> > &matrix,vector<double> &b)
{
    //矩阵和向量乘法
    int m = matrix.size();
    int n = matrix[0].size();
    if(b.size()!=n){
        std::cerr << "Dimension no Match!";
    }
    vector<double > x(n);
    int i,j;
    for(i = 0;i<m;i++){
        x[i] = 0;
        for(j = 0;j<n;j++){
            x[i]+=matrix[i][j]*b[j];
        }
    }
    return x;
}
/**
 * @brief matrix multipulication.
 * 
 * @param A 
 * @param B 
 * @return vector<vector<double> > matrix
 */
vector<vector<double> > Multiply(vector<vector<double> > A,vector<vector<double> > B)
{
    int m,n,i,j,k,n1;
    m = A.size();
    n = A[0].size();
    n1 = B[0].size();
    if(n!=B.size()){
        cerr << "Dimensions no match!";
    }
    vector<vector<double> > AB(m,vector<double> (n1));
    for(i = 0;i<m;i++){
        for(j = 0;j<n1;j++){
            AB[i][j] = 0.0;
            for(k = 0;k<n;k++){
                AB[i][j]+=A[i][k]*B[k][j];
            }
        }
    }
    return AB;
}
/**
 * @brief vector multiplies vector
 * 
 * @param a 
 * @param b 
 * @return double a*b
 */
double Multiply(vector<double> &a,vector<double> &b)
{
    double sum = 0;
    int n = a.size();
    for(int i = 0;i<n;i++){
        sum+=a[i]*b[i];
    }
    return sum;
}
/**
 * @brief Generate a certain matrix in 2.17.2
 * 
 * @param n 
 * @return vector<vector<double> > 
 */
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
/**
 * @brief LU decomposition of non-singular matrix.
 * 
 * @param matrix 
 * @param L due to 2 outputs, directly use &L as input. 
 * @param U due to 2 outputs, directly use &U as input. 
 */
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
/**
 * @brief equals to cond(A,inf) in MATLAB, but not very precise.
 * 
 * @param A 
 * @return double cond(A,inf)
 */
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
/**
 * @brief use Gauss Column Pivot to generate the inverse matrix. Equals to inv(A) in MATLAB
 * 
 * @param matrix 
 * @param invmatrix 
 */
vector<vector<double> > Inv(vector<vector<double> > &matrix)
{
    //对于病态矩阵求无穷范数条件数相比CondInf较为精确(也可能是我上面Cond写太烂了)，仅作测试使用！计算量较大，实用价值不高。
    int n = matrix[0].size();
    vector<vector<double> > invmatrix(n,vector<double> (n));
    int i,j,k;
    for(i = 0;i<n;i++){
        invmatrix[i] = GaussColumnPivot(matrix,ej(n,i));
    }
    Transposition_In_Place(invmatrix);
    return invmatrix;
}
/**
 * @brief Cholesky decomposition. Applies to symmetrical matrix only.
 * 
 * @param matrix 
 * @return vector<vector<double> > 
 */
vector<vector<double> > Cholesky(vector<vector<double> > matrix)
{
    //使用平方根法的矩阵分解，将原矩阵转化为L。
    int n = matrix[0].size();
    int i,j,k;
    for(k = 0;k<n;k++){
        matrix[k][k] = sqrt(matrix[k][k]);
        for(i = k+1;i<n;i++){
            matrix[i][k]/=matrix[k][k];
        }
        for(j = k+1;j<n;j++){
            for(i = j;i<n;i++){
                matrix[i][j]-=matrix[i][k]*matrix[j][k];
            }
        }
    }
    for(i = 1;i<n;i++){
        for(j = 0;j<i;j++){
            matrix[j][i] = 0.0;
        }
    }
    return matrix;  
}
/**
 * @brief Cholesky alternative for LDL decomposition. The diag(matrix) is D, and the remaining is L.
 * 
 * @param matrix 
 * @return vector<vector<double> > 
 */
vector<vector<double> > Cholesky_No_Sqrt(vector<vector<double> > matrix)
{
    //不作开方运算，采用平方根法的矩阵分解，将矩阵转化为LDL。
    int n = matrix[0].size();
    int i,j,k;
    double v[n];
    for(j = 0;j<n;j++){
        for(i = 0;i<j;i++){
            v[i] = matrix[j][i]*matrix[i][i];
        }
        for(i = 0;i<j;i++){
            matrix[j][j]-=matrix[j][i]*v[i];
        }
        for(k = j+1;k<n;k++){
            for(i = 0;i<j;i++){
                matrix[k][j]-=matrix[k][i]*v[i];
            }
            matrix[k][j]/=matrix[j][j];
        }
    }
    for(i = 1;i<n;i++){
        for(j = 0;j<i;j++){
            matrix[j][i] = 0.0;
        }
    }   
    return matrix;
}
/**
 * @brief Least Squares(in 2-norm)
 * 
 * @param matrix 
 * @param b 
 * @return vector <double> LS x
 */
vector <double> LS(vector<vector<double> > matrix,vector<double> b)
{
    //最小二乘解 估计
    int m = matrix.size();
    int n = matrix[0].size();
    int i,j,k;
    double sum = 0;
    vector<vector<double> > matrixT(n,vector<double >(m));
    matrixT = Transposition(matrix);
    //进行矩阵乘法，时间复杂度过高 matrixTmatrix 即书本上C
    vector<vector<double> > matrixTmatrix(n,vector<double> (n)),L(n,vector<double> (n)),LT(n,vector<double> (n));
    matrixTmatrix = Multiply(matrixT,matrix);
    vector<double> x(n),y(n),d(n);
    d = Multiply(matrixT,b);
    L = Cholesky(matrixTmatrix);
    LT = Transposition(L);
    y = LowerMatrix(L,d);
    x = UpperMatrix(LT,y);
    return x;
}
//还有一些有待改进的问题：1.对奇异矩阵没有预警机制，导致NaN出现 2.可以用符号重载简化操作