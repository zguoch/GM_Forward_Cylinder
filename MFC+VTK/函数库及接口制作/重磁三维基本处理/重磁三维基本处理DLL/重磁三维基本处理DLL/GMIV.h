//3GMIV.CPP
  //用广义逆法求解最小二乘问题
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;

#ifndef GMIV
#define GMIV
  class  gmiv
  {
  private: 
           int m, n, ka;
		   double  **a, **u, **v, **aa, *b, *x, eps;
		   double  *s, *e, *w, fg[2], cs[2];
  public:
	       gmiv (int mm, int nn)
		   {
			   int i;
			   m = mm;  n = nn;
	           a = new double*[m];   //动态分配内存空间
	           for (i=0; i<m; i++) a[i] = new double[n];
	           u = new double*[m];
	           for (i=0; i<m; i++) u[i] = new double[m];
	           v = new double*[n];
	           for (i=0; i<n; i++) v[i] = new double[n];
	           aa = new double*[n];
	           for (i=0; i<n; i++) aa[i] = new double[m];
	           b = new double[m];
	           x = new double[n];
               ka = m + 1;
	           if (m < n)  ka = n + 1;
               s = new double[ka];
               e = new double[ka];
               w = new double[ka];
		   }
	       void input ();   //由文件读入矩阵A与常数向量B以及eps
		   void input(double** A,double* B,double error);
		   void uav ();          //奇异值分解
           void ppp();            //奇异值分解的子程序
           void sss();            //奇异值分解的子程序
		   void ginv ();         //求A的广义逆
		   void a_gmiv ();       //计算最小二乘解
		   void output ();  //输出广义逆以及最小二乘解到文件并显示
		   void output(double* X);
		   ~gmiv ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<m; i++) { delete [] u[i]; }
			   delete [] u;
			   for (i=0; i<n; i++) { delete [] v[i]; }
			   delete [] v;
			   for (i=0; i<n; i++) { delete [] aa[i]; }
			   delete [] aa;
			   delete [] b, x;
			   delete [] s, e, w;
		   }
  };

#endif