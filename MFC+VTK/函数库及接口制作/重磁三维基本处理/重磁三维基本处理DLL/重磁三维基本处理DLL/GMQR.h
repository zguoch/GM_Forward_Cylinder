//3GMQR.CPP
  //用Householder变换求解线性最小二乘问题
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;

#ifndef GMQR
#define GMQR

  class  gmqr
  {
  private: 
           int m, n;
		   double  **a, **q, *b;
  public:
	       gmqr (int mm, int nn)
		   {
			   int i;
			   m = mm; n = nn;
	           a = new double*[m];   //动态分配内存空间
	           for (i=0; i<m; i++) a[i] = new double[n];
	           q = new double*[m];
	           for (i=0; i<m; i++) q[i] = new double[m];
	           b = new double[m];
		   }
	       void input ();    //由文件读入系数矩阵A与常数向量B
		   void input(double** A,double* B);
		   void qr ();       //QR分解
		   void a_gmqr ();   //计算最小二乘解
           void output ();   //输出Q和R矩阵以及最小二乘解到文件并显示
		   void output(double* X);
		   ~gmqr ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<m; i++) { delete [] q[i]; }
			   delete [] q;
			   delete [] b;
		   }
  };

#endif