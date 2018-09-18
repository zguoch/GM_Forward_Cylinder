//Ext.h
/*               Extending by folding boundary data               */

/*               a...first address of  data array                */
/*               m/n...number of line/point;                      */
/*               mx/nx...extended number of lines/points;         */
#ifndef EXT
#define EXT

void DTfunc(double **a,int m,int n,double A,double B,int sw);	
void DTari(double **a,int m,int n,double A,int sw);	
void DTfilter(double **a,int m,int n,double A,double B,double Q,int sw);	
void ext(double **a,int m,int n,int mx,int nx);
void MaxMin(double **a,int m,int n,double *MxMn);
int VectorMax(double* data,int datanum);
void intp(double **a,int m,int n,int mIns,int nIns);
double int4p(double x,double d,
			double y1,double y2,double y3,double y4);
double prep(double **zr,double **zi,int m,int n,int m1,int n1,int np,int nq);
void fft1(double *ar,double *ai,int m,int n,int key);
void dct1(double *indata,double* outdata,int n,bool key=true);
void fft(double **zr,double **zi,int m1,int n1,int me,int ne,int key);
void fact(double **s,double *x,double *y,double *g,double *b,int mt,int m,int n);
void fact1(double **s,double *x,double *y,double *g,double *b,int mt,int n,int xLt,int yLt);
void tsv(double *s,double x,double y,int m,int n);
void tsv1(double *s,double x,double y,int n,int xLt,int yLt);
void tran(double **a,int n);
void metrimult(double **a,double **b,int m,int n,int l,double **c);
void svd(double **a,int m,int n,double **u,double **v,double *q,int Index);
void Extension_Polynomial_1D(double* indatax,double* indatay,double* outdatax,double* outdatay,const int indatanum,const int outdatanum,const int N);
void Extension_Polynomial_1D(double* indatay,double* outdatay,const double dx,const int indatanum,const int outdatanum,const int N);
void quickSort(double* arr,int startPos, int endPos);
void Qsort(double a[],int low,int high);

#endif

