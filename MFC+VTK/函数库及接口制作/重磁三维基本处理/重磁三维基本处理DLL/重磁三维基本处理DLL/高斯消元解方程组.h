  //3GJDN.CPP
  //全选主元Gauss-Jordan消去法求解实系数方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;

#ifndef GSXY
#define GSXY

  class  gjdn
  {
  private: 
           int n, m;
		   double  **a, **b;
  public:
	       gjdn (int nn, int mm)
		   {
			   int i;
			   n = nn; m = mm;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           b = new double*[n];
	           for (i=0; i<n; i++) b[i] = new double[m];
		   }
	       void input ();        //从文件读入系数矩阵A以及m组常数向量B
		   void gauss_jordan ();      //执行Gauss-Jordan消去法
           void output ();       //输出结果到文件并显示
		   ~gjdn ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   delete [] b;
		   }

		   void InPut(double** A, double* B)
		   {
			   for (int i = 0; i < n; i++)
			   {
				   for (int j = 0; j < n; j++)
				   {
					   a[i][j]=A[i][j];
				   }
				   b[i][0]=B[i];
			   }
		   }

		   void OutPut(double* X)
		   {
			   for(int i=0;i<n;i++)
			   {
				   X[i]=b[i][0];
			   }
		   }
  };

#endif