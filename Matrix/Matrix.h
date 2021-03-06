/**
 * @file Matrix.h
 * @author 李杨野 (1300096763@qq.com 3190103519@zju.edu.cn)
 * @brief Header file for matrix, for test only. Completely written by myself.
 * @version 4.2
 * @date 2021-05-16
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
        if(fabs(b[i])>max)
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
 * @brief 
 * 
 * @param A 
 * @return double 
 */
double TwoNorm(vector<double> A)
{
    int n = A.size();
    int i;
    double sum = 0;
    for(i = 0;i<n;i++){
        sum+=A[i]*A[i];
    }
    return sqrt(sum);
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
/**
 * @brief Householder transformation
 * 
 * @param x 
 * @param beta 
 * @return vector<double> 
 */
vector<double> Householder(vector<double> x, double &beta)
{
    int n = x.size();
    vector<double> v(n);
    int i,j,k;
    double sum = 0;
    //double beta;
    double xinf = InfiniteNorm(x);
    for(i = 0;i<n;i++){
        x[i]/=xinf;//除以x的无穷范数，防止溢出
    }
    for(i = 1;i<n;i++){
        sum+=x[i]*x[i];
        v[i] = x[i];
    }
    if(sum == 0){
        beta = 0;
    }else{
        double alpha = sqrt(x[0]*x[0]+sum);
        if(x[0]<=0){
            v[0] = x[0] - alpha;
        }else{
            v[0] = -sum/(x[0]+alpha);
        }
        beta = 2*v[0]*v[0]/(sum + v[0]*v[0]);
        double t = v[0];
        for(j = 0;j<n;j++){
            v[j]/=t;
        }
    }
    return v;
}
/**
 * @brief Generates Householder Matrixs. Subfunction of lsqr.
 * 
 * @param v 
 * @param beta 
 * @param M 
 * @return vector<vector<double> > Householder Matrix.
 */
vector<vector<double> > HouseholderMatrix(vector<double> v,double beta,int M)
{
    int n = v.size();//v的长度，和矩阵高度不等，需要用I补齐
    int i,j,k;
    vector<vector<double> >H(M,vector<double > (M));
    for(i = 0;i<M-n;i++){
        H[i][i] = 1;
    }
    for(i = M-n;i<M;i++){
        for(j = M-n;j<M;j++){
            if(i == j){
                H[i][j] = 1-beta*v[i-M+n]*v[j-M+n];
            }else{
                H[i][j] = -beta*v[i-M+n]*v[j-M+n];
            }
        }
    }
    return H;
}
/**
 * @brief QR composition. Subfunction of lsqr.
 * 
 * @param matrix Used to store R in the later procession.
 * @return vector<vector<double> > Q
 */
vector<vector<double> > QR(vector<vector<double> > &matrix)
{
    int i,j,k;
    int m = matrix.size();
    int n = matrix[0].size();
    vector<double> d(n);
    double beta;
    for(j = 0;j<n;j++){
        if(j<m){
            vector<double> v(m-j),x(m-j);
            for(i = 0;i<m-j;i++){
                x[i] = matrix[j+i][j];
            }
            v = Householder(x,beta);
            //delete [] x;
            vector<vector<double> >temp(m-j,vector<double > (n-j));
            for(i = j;i<m;i++){
                for(k = j;k<n;k++){
                    temp[i-j][k-j] = matrix[i][k];
                }
            }
            vector<vector<double> >H(m-j,vector<double > (m-j));
            for(i = 0;i<m-j;i++){
                for(k = 0;k<m-j;k++){
                    if(i == k){
                        H[i][k] = 1-beta*v[i]*v[k];
                    }else{
                        H[i][k] = -beta*v[i]*v[k];
                    }
                }
            }
            temp = Multiply(H,temp);
            for(i = j;i<m;i++){
                for(k = j;k<n;k++){
                    matrix[i][k] = temp[i-j][k-j];
                }
            }
            //delete [] H;
            //delete [] temp;
            d[j] = beta;
            for(k = j+1;k<m;k++){
                matrix[k][j] = v[k-j];
            }
        }
    }
    vector<vector<double> >Q(m,vector<double > (m));
    for(j = 0;j<n;j++){
        vector<double> r(m-j);
        r[0] = 1;
        for(i = j+1;i<m;i++){
            r[i-j] = matrix[i][j];
            matrix[i][j] = 0;
        }
        vector<vector<double> >H(m,vector<double > (m));
        H = HouseholderMatrix(r,d[j],m);
        if(j == 0){
            Q = H;
        }else{
            Q = Multiply(Q,H);
        }
    }
    return Q;
}
/**
 * @brief Equals to lsqr(A,b) in MATLAB.
 * 
 * @param matrix 
 * @param b 
 * @return vector <double> 
 */
vector <double> QRLS(vector<vector<double> > matrix,vector<double> b)
{
    //最小二乘解 估计
    int m = matrix.size();
    int n = matrix[0].size();
    vector<vector<double> >Q(m,vector<double > (m));
    Q = QR(matrix);
    vector<double> x(n),c(n);
    int i,j,k;
    Transposition_In_Place(Q);
    b = Multiply(Q,b);
    for(i = 0;i<n;i++){
        c[i] = b[i];
    }
    x = UpperMatrix(matrix,c);
    return x;
}
/**
 * @brief Jacobi Iteration
 * 
 * @param matrix 
 * @param b 
 * @return vector <double> 
 */
vector <double> Jacobi(vector<vector<double> > matrix,vector<double> b)
{
    //非奇异，暂定m=n
    int n = matrix[0].size();
    vector<double> x(n),x1(n),t(n);
    int i,j,k = 0;
    do{
        x = x1;//记录x_{k-1}的值，便于误差估计
        for(i = 0;i<n;i++){
            double sum = 0;
            for(int j=0;j<n;j++){
                if(i!=j){
                    sum += matrix[i][j]*x[j];
                }
            }
            x1[i]=(b[i]-sum)/matrix[i][i];
        }
        for(i = 0;i<n;i++){
            t[i] = x[i]-x1[i];
        }
        k++;
    }while(k<100&&InfiniteNorm(t)>=0.00001);
    if(k>100){
        cout << "迭代次数过多，请确认矩阵是否收敛！"<< endl;
    }
    return x;
}
/**
 * @brief Gauss-Seidel Iteration
 * 
 * @param matrix 
 * @param b 
 * @return vector <double> 
 */
vector <double> GaussSeidel(vector<vector<double> > matrix,vector<double> b)
{
    //非奇异，暂定m=n
    int n = matrix[0].size();
    vector<double> x(n),x1(n),t(n);
    int i,j,k = 0,l;
    do{
        x = x1;//记录x_{k-1}的值，便于误差估计
        for(i = 0;i<n;i++){
            double sum = 0;
            for(int j=0;j<n;j++){
                if(i!=j){
                    sum += matrix[i][j]*x1[j];//这里用x1而不是x，和Jacobi有所不同
                }
            }
            x1[i]=(b[i]-sum)/matrix[i][i];
        }
        for(i = 0;i<n;i++){
            t[i] = x[i]-x1[i];
        }
        k++;
    }while(k<100&&InfiniteNorm(t)>=0.00001);
    if(k>100){
        cout << "迭代次数过多，请确认矩阵是否收敛！"<< endl;
    }
    return x1;
}
/**
 * @brief SOR Iteration
 * 
 * @param matrix 
 * @param b 
 * @param omega 
 * @return vector <double> 
 */
vector <double> SOR(vector<vector<double> > matrix,vector<double> b,double omega)
{
    //omega由系数矩阵决定，取最佳松弛因子
    int n = matrix[0].size();
    vector<double> x(n),x1(n),t(n);
    int i,j,k = 0,l;
    do{
        x = x1;//记录x_{k-1}的值，便于误差估计
        for(i = 0;i<n;i++){
            double sum = 0;
            for(int j=0;j<n;j++){
                if(i!=j){
                    sum += matrix[i][j]*x1[j];
                }
            }
            x1[i]=omega*(b[i]-sum)/matrix[i][i]+(1-omega)*x[i];//应该是这样吧…？
        }
        for(i = 0;i<n;i++){
            t[i] = x[i]-x1[i];
        }
        k++;
    }while(k<100&&InfiniteNorm(t)>=0.00001);
    if(k>100){
        cout << "迭代次数过多，请确认矩阵是否收敛！"<< endl;
    }
    return x1;
}
/**
 * @brief double the sum of b[i]^2
 * 
 * @param b 
 * @return 
 */
double Square(vector<double> b)
{
    double sum = 0;
    int n = b.size();
    for(int i = 0;i<n;i++){
        sum+=b[i]*b[i];
    }
    return sum;
}
/**
 * @brief Conjugate Gradient method.
 * 
 * @param A 
 * @param b 
 * @return vector <double> 
 */
vector <double> CG(vector<vector<double> > A,vector<double> b)
{
    int n = A[0].size();
    int i,j;
    double alpha,beta;
    vector <double> x(n),r(n),t(n),p(n),temp(n);//初值为0
    r = b;//r_0 = b-Ax_0 = b.
    t = b;
    int k = 0;
    while(OneNorm(t)>0.0001 && k<1000){
        k++;
        if(k == 1){
            p = t;
        }else{
            beta = Square(t)/Square(r);
            for(i = 0;i<n;i++){
                p[i] = t[i]+beta*p[i];//p_k-1 -> p_k-2
            }
        }
        r = t;
        temp = Multiply(A,p);
        alpha = Square(t)/Multiply(p,temp);
        for(i = 0;i<n;i++){
            x[i] += alpha*p[i];
            for(j = 0;j<n;j++){
                t[i] -= alpha*A[i][j]*p[j];//r_k-1 -> r_k
            }
        }//t = rk
    }
    if(k>=1000)cout << "too much iteration!";
    return x;
}