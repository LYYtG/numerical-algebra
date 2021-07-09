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
#include <complex>
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

vector<double> Rev_PowerMethod(vector<vector<double > > A,double lambda)
{
    int n = A[0].size();
    vector<double> u(n),v(n);
    double t;
    int i,j,k = 0;
    for(i = 0;i<n;i++)
    {
        u[i] = 1.0;//Initialize...
        A[i][i]-=lambda;
    }
    do{
        v = u;//v = u_k-1        
        u = GaussColumnPivot(A,u);
        t = TwoNorm(u);
        for(i = 0;i<n;i++){
            u[i]/=t;
        }
        //Print(u);
        k++;
    }while(uv(u,v)&&k<100);
    if(k=100){
        //cout << "???";
    }
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
/**
 * @brief Convert A to an upper Hessenberg Matrix
 * 
 * @param A 
 */
void Hessenberg(vector<vector<double> > &A)
{
    int n = A[0].size();
    double u = 1e-15;
    vector<vector<double> > H(n,vector<double> (n));
    vector<vector<double> > Q(n,vector<double> (n));
    int i,j,k;
    double beta;
    k = 0;
    for(k = 0;k<n-2;k++)
    {
        vector<double> h(n-k-1),v(n-k-1);
        for(i = 0;i<n-k-1;i++){
            h[i] = A[k+i+1][k];
        }
        v = Householder(h,beta);
        H = HouseholderMatrix(v,beta,n);
        A = Multiply(H,A);
        A = Multiply(A,H);
    }
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++)
        {
            if(fabs(A[i][j])<u){
                A[i][j] = 0;
            }
        }
    }
    //Print(A);
}
vector<vector<double> > Francis(vector<vector<double> > H)//not &H!
{
    //H is an upper Hessenberg matrix.
    int n = H[0].size();
    if(n==1){
        vector<vector<double> > I(1,vector<double> (1));
        I[0][0] = 1;
        return I;
    }
    int m = n-1;
    int i,j,k,q,r;
    double u = 1e-15;
    double s = H[m-1][m-1] + H[n-1][n-1];
    double t = H[m-1][m-1]*H[n-1][n-1]-H[m-1][n-1]*H[n-1][m-1];
    double beta;
    vector<double> x(3);
    vector<double> v(3);
    vector<vector<double> > Q(n,vector<double> (n));
    x[0] = H[0][0]*H[0][0]+H[0][1]*H[1][0]-s*H[0][0]+t;
    x[1] = H[1][0]*(H[0][0]+H[1][1]-s);
    if(n>=3){
        x[2] = H[1][0]*H[2][1];
    }
    for(k = -1;k<n-3;k++)
    {
        v = Householder(x,beta);//Householder vector.
        if(k>0){
            q = k;//first row.
        }else{
            q = 0;
        }
        vector<vector<double> > H1(3,vector<double> (3));
        H1 = HouseholderMatrix(v,beta,3);//Householder Matrix.
        vector<vector<double> > P(n,vector<double> (n));
        for(i = k+1;i<=k+3;i++){
            for(j = k+1;j<=k+3;j++){
                P[i][j] = H1[i-k-1][j-k-1];
            }
        }
        for(i = 0;i<k+1;i++){
            P[i][i] = 1;
        }
        for(i = k+4;i<n;i++){
            P[i][i] = 1;
        }
        //Print(P);
        if(k == -1){
            Q = P;
        }else{
            Q = Multiply(Q,P);
        }
        H = Multiply(P,H);
        H = Multiply(H,P);
        x[0] = H[k+2][k+1];
        x[1] = H[k+3][k+1];
        if(k<n-4){
            x[2] = H[k+4][k+1];
        }
        /*
        vector<vector<double> > T(3,vector<double > (n-q));//temp Matrix.
        for(i = k+1;i<k+4;i++){
            for(j = q;j<n;j++){
                T[i-k-1][j-q] = H[i][j];
            }
        }
        T = Multiply(H1,T);
        for(i = k+1;i<k+4;i++){
            for(j = q;j<n;j++){
                H[i][j] = T[i-k-1][j-q];
            }
        }
        if(k+4<n-1){
            r = k+4;//last column.
        }else{
            r = n-1;
        }
        vector<vector<double> > T2(r+1,vector<double > (3));
        for(i = 0;i<=r;i++){
            for(j = k+1;j<k+4;j++){
                T2[i][j-k-1] = H[i][j];
            }
        }
        T2 = Multiply(T2,H1);
        for(i = 0;i<=r;i++){
            for(j = k+1;j<k+4;j++){
                 H[i][j] = T2[i][j-k-1];
            }
        }
        x[0] = H[k+2][k+1];
        x[1] = H[k+3][k+1];
        if(k<n-4){
            x[2] = H[k+4][k+1];
        }
        */
    }
    
    vector<double> x1(2);
    x1[0] = x[0];
    x1[1] = x[1];
    vector<double> v1(2);
    v1 = Householder(x1,beta);
    vector<vector<double> > H2(2,vector<double> (2));
    H2 = HouseholderMatrix(v1,beta,2);
    vector<vector<double> > P1(n,vector<double> (n));
    for(i = n-2;i<n;i++){
        for(j = n-2;j<n;j++){
            P1[i][j] = H2[i-n+2][j-n+2];
        }
    }
    for(i = 0;i<n-2;i++){
        P1[i][i] = 1;
    }
        H = Multiply(P1,H);
        H = Multiply(H,P1);
        if(n == 2){
            Q = P1;
        }else{
            Q = Multiply(Q,P1);
        }
        
    /*
    vector<vector<double> > Tr(2,vector<double > (3));
    vector<vector<double> > Tc(n,vector<double > (2));

    for(j = 0;j<2;j++){
        for(i = 0;i<n;i++){        
            Tc[i][j] = H[i][n-2+j];
        }
        for(i = 0;i<3;i++){
            Tr[j][i] = H[n-2+j][n-3+i];
        }
    }   
    Tr = Multiply(H2,Tr);
    Tc = Multiply(Tc,H2);
    for(j = 0;j<2;j++){
        for(i = 0;i<n;i++){        
            H[i][n-2+j] = Tc[i][j];
        }
        for(i = 0;i<3;i++){
            H[n-2+j][n-3+i] = Tr[j][i];
        }
    }
    */
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++)
        {
            if(fabs(H[i][j])<u){
                H[i][j] = 0;
            }
        }
    }
    //Print(H);
    return Q;
}
/**
 * @brief subfunction of implicit QR
 * 
 * @param A 
 * @param m 
 * @return true 
 * @return false 
 */
bool IsUpper(vector<vector<double> > &A,int m)
{
    int n = A[0].size();
    if(A[n-m][n-m-1]!=0)return false;
    if(m<=2)return true;

    int i,j,k;
    for(i = n-m+1;i<n;i++){
        if(A[i][i-1]!=0){
            if(i == n-m+1&&A[i+1][i]!=0){
                return false;
            }else if(i == n-1&&A[i-2][i-1]!=0){
                return false;
            }else if(A[i+1][i]!=0||A[i-2][i-1]!=0){
                return false;
            }
        }
    }
    return true;
}
vector<double> EigVector(vector<vector<double> > A,double lambda)
{
    int i,j,k = 0;
    int n = A[0].size();
    vector<double> b(n),v(n);
    for(i = 0;i<n;i++){
        A[i][i]-=lambda;
    }
    v = GaussColumnPivot(A,b);
    return v;
}//I can't figure out how to use C++ to manipulate complex...if I use class "complex", I'm forced to rewrite whole Matrix.h\
     to get the complex eigenvector...so I simply gave up.  
/**
 * @brief Implicit QR
 * 
 * @param A 
 */
vector<complex<double>> ImplicitQR(vector<vector<double> > &A)
{
    int i,j,k = 0;
    int m = 0,l = 0;
    int n = A[0].size();
    vector<complex<double>> eig(n);
    complex<double> temp;
    Hessenberg(A);
    double u = 1e-22;
    while((m!=n||l!=n)&&k<100){
        for(i = 1;i<n;i++){
            if(fabs(A[i][i-1])<u*(fabs(A[i][i])+fabs(A[i-1][i-1]))){
                A[i][i-1] = 0;
            }
        }
        m = 1;
        while(IsUpper(A,m)){
            m++;
        }
        m--;
        for(i = n-m-1;i>=1;i--){
            if(A[i][i-1] == 0){
                break;
            }
        }

        l = i;
        //cout << "l:"<< l<< ",m:"<< m<<endl;
        vector<vector<double> > H22(n-m-l,vector<double> (n-m-l));
        vector<vector<double> > P1(n-m-l,vector<double> (n-m-l));
        vector<vector<double> > P(n,vector<double> (n));
        vector<vector<double> > PT(n,vector<double> (n));
        //vector<vector<double> > H12(l,vector<double> (n-m-l));
        //vector<vector<double> > H23(n-m-l,vector<double> (l));
        for(i = 0;i<n-m-l;i++){
            for(j = 0;j<n-m-l;j++){
                H22[i][j] = A[i+l][j+l];
            }
        }
        //Print(H22);
        P1 = Francis(H22);
        for(i = l;i<n-m;i++){
            for(j = l;j<n-m;j++){
                P[i][j] = P1[i-l][j-l];
            }
        }
        for(i = 0;i<l;i++){
            P[i][i] = 1;
        }
        for(i = n-m;i<n;i++){
            P[i][i] = 1;
        }
        //Print(P);
        PT = Transposition(P);
        A = Multiply(PT,A);
        A = Multiply(A,P);
        for(i = 0;i<n;i++){
        for(j = 0;j<n;j++)
        {
            if(fabs(A[i][j])<u){
                A[i][j] = 0;
            }
        }
    }
        //Print(A);
        k++;
    }  
    
   //Print(A);
   double a = 1,b,c,R,C;
   double g;
   for(i = 0;i<n-1;i++){
       if(A[i+1][i] == 0){
           R = A[i][i];
           eig[i] = complex<double>(R,0);
       }else{
           b = -A[i][i]-A[i+1][i+1];
           c = A[i][i]*A[i+1][i+1]-A[i+1][i]*A[i][i+1];
           g = b*b-4*a*c;
           if(g>=0){
               R = 0.5*(-b+sqrt(g));
               eig[i] = complex<double>(R,0);
               R = 0.5*(-b-sqrt(g));
               eig[++i] = complex<double>(R,0);
           }else{
               g = -g;
               R = -0.5*b;
               C = 0.5*sqrt(g);
               eig[i] = complex<double>(R,C);
               eig[++i] = complex<double>(R,-C);
           }
       }
   }
   if(i == n-1){
       eig[i] = complex<double>(A[i][i],0);
   }
   return eig;
}
vector<vector<double> > Givens(double a,double b ,int j,int k,int n)
{
    double c,s,t;
    if(b == 0){
        c = 1;
        s = 0;
    }else{
        if(fabs(b)>fabs(a)){
            t = a/b;s = 1/(sqrt(1+t*t));c = s*t;
        }else{
            t = b/a;s = 1/sqrt(1+t*t);s = c*t;
        }
    }
    vector<vector<double> > A(n,vector<double> (n));
    int i;
    for(i = 0;i<n;i++){
        A[i][i] = 1;
    }
    A[j][j] = c;
    A[j][k] = s;
    A[k][j] = -s;
    A[k][k] = c;
    return A;
}
vector<vector<double> > Wilkinson(vector<vector<double> > T)
{
    int n = T[0].size();
    if(n==1){
        vector<vector<double> > I(1,vector<double> (1));
        I[0][0] = 1;
        return I;
    }
    vector<vector<double> > Q(n,vector<double> (n));
    int i,j,k;
    double d = (T[n-2][n-2]-T[n-1][n-1])/2;
    double u = T[n-1][n-1]-T[n-1][n-2]*T[n-1][n-2]/(d+Sign(d)*sqrt(d*d+T[n-1][n-2]*T[n-1][n-2]));
    double x = T[0][0]-u;
    double z = T[1][0];
    for(k = 0;k<n-1;k++){
        vector<vector<double> > G(n,vector<double> (n));
        vector<vector<double> > GT(n,vector<double> (n));
        G = Givens(x,z,k,k+1,n);
        GT = Transposition(G);
        T = Multiply(G,T);
        T = Multiply(T,GT);
        if(k == 0){
            Q = G;
        }else{
            Q = Multiply(G,Q);
        }
        if(k<n-2){
            x = T[k+1][k];
            z = T[k+2][k];
        }
    }
    return Q;
}
bool IsUpper2(vector<vector<double> > &A,int m)
{
    int n = A[0].size();
    int i,j,k;
    for(i = n-m+1;i<n;i++){
        if(A[i][i-1]!=0){
            return false;
        }
    }
    return true;
}
vector<complex<double>> ImplicitQR2(vector<vector<double> > &A)
{
    int i,j,k = 0;
    int m = 0,l = 0;
    int n = A[0].size();
    vector<complex<double>> eig(n);
    complex<double> temp;
    double u = 1e-14;
    while((m!=n||l!=n)&&k<100){
        for(i = 1;i<n;i++){
            if(fabs(A[i][i-1])<u*(fabs(A[i][i])+fabs(A[i-1][i-1]))){
                A[i][i-1] = 0;
                A[i-1][i] = 0;
            }
        }
        m = 1;
        while(IsUpper2(A,m)){
            m++;
        }
        m--;
        for(i = n-m-1;i>=1;i--){
            if(A[i][i-1] == 0){
                break;
            }
        }

        l = i;
        //cout << "l:"<< l<< ",m:"<< m<<endl;
        vector<vector<double> > H22(n-m-l,vector<double> (n-m-l));
        vector<vector<double> > P1(n-m-l,vector<double> (n-m-l));
        vector<vector<double> > P(n,vector<double> (n));
        vector<vector<double> > PT(n,vector<double> (n));
        //vector<vector<double> > H12(l,vector<double> (n-m-l));
        //vector<vector<double> > H23(n-m-l,vector<double> (l));
        for(i = 0;i<n-m-l;i++){
            for(j = 0;j<n-m-l;j++){
                H22[i][j] = A[i+l][j+l];
            }
        }
        //Print(H22);
        P1 = Wilkinson(H22);
        for(i = l;i<n-m;i++){
            for(j = l;j<n-m;j++){
                P[i][j] = P1[i-l][j-l];
            }
        }
        for(i = 0;i<l;i++){
            P[i][i] = 1;
        }
        for(i = n-m;i<n;i++){
            P[i][i] = 1;
        }
        //Print(P);
        PT = Transposition(P);
        A = Multiply(P,A);
        A = Multiply(A,PT);
        
        for(i = 0;i<n;i++){
        for(j = 0;j<n;j++)
        {
            if(fabs(A[i][j])<u){
                A[i][j] = 0;
            }
        }
    }
        //Print(A);
        k++;
    }  
    
   //Print(A);
   double a = 1,b,c,R,C;
   double g;
   for(i = 0;i<n-1;i++){
       if(A[i+1][i] == 0){
           R = A[i][i];
           eig[i] = complex<double>(R,0);
       }else{
           b = -A[i][i]-A[i+1][i+1];
           c = A[i][i]*A[i+1][i+1]-A[i+1][i]*A[i][i+1];
           g = b*b-4*a*c;
           if(g>=0){
               R = 0.5*(-b+sqrt(g));
               eig[i] = complex<double>(R,0);
               R = 0.5*(-b-sqrt(g));
               eig[++i] = complex<double>(R,0);
           }else{
               g = -g;
               R = -0.5*b;
               C = 0.5*sqrt(g);
               eig[i] = complex<double>(R,C);
               eig[++i] = complex<double>(R,-C);
           }
       }
   }
   if(i == n-1){
       eig[i] = complex<double>(A[i][i],0);
   }
   return eig;
}

vector<vector<double> > JacobiEig(vector<vector<double> > A)
{
    int i,j,k,p,q;
    int n = A[0].size();
    vector<double> temp1(n),temp2(n);
    vector<vector<double> > Q(n,vector<double> (n));
    double c,s,delta = 0,t1,t2;
    for(i = 0;i<n;i++){
        for(j = 0;j<n;j++){
            if(i!=j){
                delta+=A[i][j]*A[i][j];
            }else{
                Q[i][i] = 1;
            }           
        }
    }
    delta = sqrt(delta);
    while(delta>=1e-13){
        int flag = 0;
        for(p = 0;p<=n-1;p++){
            for(q = 0;q<=n-1;q++){
                if(fabs(A[p][q])>=delta&&p!=q){
                    double tao = (A[q][q]-A[p][p])/(2*A[p][q]);
                    double t = Sign(tao)/(fabs(tao)+sqrt(1+tao*tao));
                    c = 1.0/sqrt(1+t*t);
                    s = c*t;
                    //cout << "c:"<< c<<" s:"<<s<< endl;
                    for(i = 0;i<n;i++){
                        if(i == p){
                            temp1[i] = c*c*A[p][p]-2*s*c*A[p][q]+s*s*A[q][q];
                            temp2[i] = 0;
                        }else if(i == q){
                            temp1[i] = 0;
                            temp2[i] = s*s*A[p][p]+2*s*c*A[p][q]+c*c*A[q][q];
                        }else{
                            temp1[i] = c*A[i][p]-s*A[i][q];
                            temp2[i] = s*A[i][p]+c*A[i][q];
                        }
                    }
                    for(i = 0;i<n;i++){
                        A[i][p] = temp1[i];
                        A[p][i] = temp1[i];
                        A[i][q] = temp2[i];
                        A[q][i] = temp2[i];
                        t1 = Q[i][p];
                        t2 = Q[i][q];
                        Q[i][p] = c*t1-s*t2;
                        Q[i][q] = s*t1+c*t2;//construct Q.
                    }
                    flag = 1;
                }
            }
        }
        if(!flag){
            delta/=n;
        }
    }
    for(i = 0;i<n;i++){
        cout << A[i][i]<< " ";
    }
    cout << endl;
    return Q;
}
int AltCount(double u,vector<double> x,vector<double> y)
{
    int n = x.size();
    double q = x[0]-u;
    int s = 0;
    int i,j,k;
    for(k = 0;k<n;k++){
        if(q<0){
            s++;
        }
        if(k<n-1){
            if(q == 0){
                q = fabs(y[k+1])+1e-14;
            }
            q = x[k+1]-u-y[k+1]*y[k+1]/q;
        }
    }
    return s;
}
double Dich(vector<vector<double > > A,int m)
{
    int n = A[0].size();
    int i,j,k;
    vector<double> x(n);
    vector<double> y(n);
    double u = 1e-14;
    double upper = InfiniteNorm(A)+u;
    double lower = -upper;//Eigenvalue belongs to [l,u].
    for(i = 0;i<n;i++){
        x[i] = A[i][i];
        if(i == 0){
            y[i] = 0;
        }else{
            y[i] = A[i-1][i];
        }
    }//x = {a_1,a_2,...,a_n} y = {0,b_1,...,b_n}
    while(fabs(upper-lower)>=u){
        //cout << upper << " "<<lower<<" " <<AltCount(upper,x,y)<<endl;
        double mid = (upper+lower)/2;
        if(AltCount(mid,x,y)>=m){
            upper = mid;
        }else{
            lower = mid;
        }
        
    }
    return (upper+lower)/2;
}