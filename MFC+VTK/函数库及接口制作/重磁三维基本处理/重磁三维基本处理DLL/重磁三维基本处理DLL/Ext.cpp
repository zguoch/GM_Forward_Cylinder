//Ext.cpp
/*               Extending by folding boundary data               */

/*               a...first address of  data array                */
/*               m/n...number of line/point;                      */
/*               mx/nx...extended number of lines/points;         */

//#include "stdafx.h"
#include "proj.h"

void DTfilter(double **a,int m,int n,
			  double A,double B,double Q,int sw)	
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			switch(sw)
			{
			case 0:
				{
					if(a[i][j]>=A)
						a[i][j]=Q;
					break;
				}
			case 1:
				{
					if(a[i][j]<=A)
						a[i][j]=Q;
					break;
				}
			case 2:
				{
					if((a[i][j]<=A)||(a[i][j]>=B))
						a[i][j]=Q;
					break;
				}
			case 3:
				{
					if((a[i][j]>=A)&&(a[i][j]<=B))
						a[i][j]=Q;
					break;
				}
			}
		}
	}
}

void DTari(double **a,int m,int n,double A,int sw)	
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			switch(sw)
			{
			case 0:	{   a[i][j]+=A;	break;	}
			case 1:	{   a[i][j]*=A;	break;	}
			case 2:	{   a[i][j]=A/(double)(a[i][j]+1e-40);	break;	}
			case 3:	{   a[i][j]*=a[i][j];   break;	}
			}
		}
	}
}

void DTfunc(double **a,int m,int n,double A,double B,int sw)	
{
	int i,j;
	double min,p;
	if(sw==2)
	{
		min=(double)fabs(a[0][0]);
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				if(fabs(a[i][j])<min) min=(double)fabs(a[i][j]);
			}
		}
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			switch(sw)
			{
			case 0:	{a[i][j]=A+B*a[i][j];break;}
			case 1:	{a[i][j]=A*(double)exp((double)(B*a[i][j]));break;}
			case 2:	
				{	
					p=(double)fabs(a[i][j])/(double)(a[i][j]+1e-80);
					a[i][j]=A*p*(double)log((double)(B*(fabs(a[i][j])+1e-80)/(min+1e-80)));
					break;
				}
			case 3:	{a[i][j]=A*(double)pow((double)(a[i][j]+1e-80),(double)B);break;}
			case 4:	{a[i][j]=A*(double)atan((double)(B*a[i][j]));break;}
			}
		}
	}
}


///////////////////////////////////////////////////////
//                Preparation & Extending Data                    
//
//    zr/zi...first address of real/image data array              
//   m/n...number of line/point; m1/n1...extended m/n             
//   me/ne...{m1/n1=pow(2,me/ne)}; np/nq...{m1/n1=2*(np/nq)+m/n}  
double prep(double **zr,double **zi,int m,int n,int m1,int n1,int np,int nq)
{  
	int i,j; double pjg;
	pjg=(double)0;
	if(m==1)				//计算数据边界均值
	{	pjg=(zr[0][0]+zr[0][n-1])/(double)2.;
	}
	else
	{	for(j=0;j<n;j++)
		{   pjg+=(zr[0][j]+zr[m-1][j]);
		}
		for(i=1;i<(m-1);i++)
		{   pjg+=(zr[i][0]+zr[i][n-1]);
		}
		pjg=pjg/(double)(2*m+2*n-4);			//四个角点重复计算了一次
	}
	for(i=(m-1);i>=0;i--)
	{   for(j=(n-1);j>=0;j--)
		{   zr[i+np][j+nq]=zr[i][j]-pjg;
	}   }
//---------------------------------------------
//            extending data                   
//---------------------------------------------
	
	if(m>1)
	{   for(i=np;i<(m+np);i++)
		{   zr[i][n1-1]=(double)0;
			for(j=0;j<nq;j++)
			{	zr[i][j]=zr[i][nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
				zr[i][j+n+nq]=zr[i][n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/2;
		}   }
		for(j=0;j<n1;j++)
		{ 	zr[m1-1][j]=(double)0;
			for(i=0;i<np;i++)
			{   zr[i][j]=zr[np][j]*(double)(1-cos(((double)PI*(double)i/(double)np)))/(double)2;
				zr[i+m+np][j]=zr[m+np-1][j]*(double)(1+cos(((double)PI*(double)(i+1)/(double)(m1-m-np))))/(double)2;
		}   }
		for(i=0;i<m1;i++)
		{   for(j=0;j<n1;j++)
			{   zi[i][j]=(double)0;
		}   }
	}
	else
	{   for(j=0;j<nq;j++)
		{   zr[0][j]=zr[0][nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
			zr[0][j+n+nq]=zr[0][n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
		}
		for(j=0;j<n1;j++)
		{  zi[0][j]=(double)0; }
	}

	////4. 把二维数组转换为一位数组
	//int datanum=m1*n1;
	//double *mm=new double[datanum];
	//for (int i = 0; i < m1; i++)
	//{
	//	for (int j = 0; j < n1; j++)
	//	{
	//		mm[j+i*n1]=zr[i][j];
	//	}
	//}
	//WriteGrd("扩边后.grd",mm,datanum,m1,n1,-2000,2000,-2000,2000);

	return pjg;
}
///////////////////////////////////////////
////To extend the data in array a[][] from
//// a[m][n] to a[m+2*mx][n+2*nx]
void ext(double **a,int m,int n,int mx,int nx)
{   int i,j,i1,j1;

	for(i=(m+mx-1);i>=mx;i--)
		for(j=(n+nx-1);j>=nx;j--)
			a[i][j]=a[i-mx][j-nx];
	for(i=mx;i<(m+mx);i++)
	{   for(j=0;j<nx;j++)
		{   j1=n+nx+j;
			a[i][j]=a[i][2*nx-j]; a[i][j1]=a[i][n+nx-2-j];
	}   }
	for(j=0;j<n+2*nx;j++)
	{   for(i=0;i<mx;i++)
		{   i1=m+mx+i;
			a[i][j]=a[2*mx-i][j]; a[i1][j]=a[m+mx-2-i][j];
	}   }
 	return ;
}

////////////////////////////////////////////
//To find the Maximun and Minimun of present array data 
void MaxMin(double **a,int m,int n,double *MxMn)
{
	double p;
 	MxMn[0]=MxMn[1]=a[0][0]; p=(double)0;
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			p+=a[i][j];
			if(MxMn[1]<a[i][j])
				MxMn[1]=a[i][j];
			else if(MxMn[0]>a[i][j])
				MxMn[0]=a[i][j];
		}
	}
	MxMn[2]=p/(double)(m*n+1e-40);
}
int VectorMax(double* data,int datanum)
{
	int index=0;
	double max=data[index];
	for (int i = 0; i < datanum; i++)
	{
		if (max<data[i])
		{
			max=data[i];
			index=i;
		}
	}
	return index;
}
void intp(double **a,int m,int n,int mi,int ni)
 {	int i,j,k; double d,x;
//	mi=(int)(dx0/dx+.0001)
// ---number of inserted rows + 1 
//	ni=(int)(dy0/dy+.0001)
// ---number of inserted columns + 1 
 	ext(a,m,n,1,1);
	for(i=m+1;i>=0;i--)
		for(j=n+1;j>=0;j--)
			a[i*mi][j*ni]=a[i][j];
 	d=(double)ni; //(0, d, 2d, 3d,) are y_values knawn
	for(i=0;i<=(m+1)*mi;i+=mi)
	{   for(j=0;j<=(n-2)*ni;j+=ni)
		{	for(k=1;k<ni;k++)
			{	x=(double)k;	// y=f(x) is the value required
		a[i][j+ni+k]=int4p(x,d,a[i][j],a[i][j+ni],a[i][j+2*ni],a[i][j+3*ni]);
	}
	}	}
	d=(double)mi; //(0, d, 2d, 3d,) are y_values knawn
	for(j=ni;j<=n*ni;j++)
	{   for(i=0;i<=(m-2)*mi;i+=mi)
		{	for(k=1;k<mi;k++)
			{	x=(double)k;	// y=f(x) is the value required
		a[i+mi+k][j]=int4p(x,d,a[i][j],a[i+mi][j],a[i+2*mi][j],a[i+3*mi][j]);
	}	}	}
	for(i=0;i<(m-1)*mi+1;i++)
		for(j=0;j<(n-1)*ni+1;j++)
			a[i][j]=a[i+mi][j+ni];
	return;
}


double int4p(double x,double d,double y1,double y2,double y3,double y4)
 {	double z,t;
	t=y1*x*(x-d)/(2*d*d)-y2*(x+d)*(x-d)/(d*d)+y3*(x+d)*x/(2*d*d);
	z=y2*(x-d)*(x-2*d)/(2*d*d)-y3*x*(x-2*d)/(d*d)+y4*x*(x-d)/(2*d*d);
	z=(z+t)/2;
	return(z);
}
///////////////////////////////////////////////////////
/*                   Fast-Fourier_Transform                        */

/*    zr/zi...first address of real/image data array               */
/*    m/n...number of line/point; m1/n1...extended m/n             */
/*    me/ne...{m1/n1=pow(2,me/ne)};                                */
/*    key... 0 for +FFT, 1 for -FFT                                */


void fft(double **zr,double **zi,int m1,int n1,int me,int ne,int key)
{   double *cr=new double [m1]; 
	double *ci=new double [m1];
	int i,j;
	double *ir0,*ii0,*ir1,*ii1;

	if(m1==1)
	{	ir0=*(zr+0)+0; ii0=*(zi+0)+0;
		fft1(ir0,ii0,ne,n1,key);
		return ;
	}
	for(i=0;i<m1;i++)
		{	ir0=*(zr+i)+0; ii0=*(zi+i)+0;  ir1=*(zr+0)+i+8;
			fft1(ir0,ii0,ne,n1,key);
		}
	for(j=0;j<n1;j++)
	{   for(i=0;i<m1;i++)
		{	cr[i]=*(*(zr+i)+j); ci[i]=*(*(zi+i)+j);
		}
		fft1(cr,ci,me,m1,key);
		for(i=0;i<m1;i++)
		{	ir1=*(zr+i)+j; ii1=*(zi+i)+j;
			*ir1=cr[i]; *ii1=ci[i];
	}	}
	delete cr;
	delete ci;
	return ;
}
/*  1-Dimension Fast Fourier Transform (FFT)  */

void fft1(double *ar,double *ai,int m,int n,int key)
{   int i,j,k,k1,l1,l2,ip;
	double *ir,*ii,*i1,*i2;
	double f1,t1,t2,v,u,w1,w2,d,vu,pi1;
	f1=(double)-1;
	if(key==0) f1=(double)1;
	j=0;
	for(i=0,ir=ar,ii=ai;i<n;i++,ir++,ii++)
		{	if(i<j)
			{	i1=ar+j; i2=ai+j;
				t1=*i1; t2=*i2;	*i1=*ir; *i2=*ii; *ir=t1; *ii=t2;
			}
			k1=n/2;
			while((j>=k1)&&(k1!=0))
			{	j=j-k1; k1=k1/2;
			}
			j=j+k1;
		}
	for(j=1;j<m+1;j++)
	{	l1=2;
		for (i=1;i<j;i++) l1=2*l1;
		l2=l1/2; pi1=(double)PI/(double)l2; v=(double)1; u=(double)0;
		w1=(double)cos(pi1); w2=(double)sin(pi1);
		for(k=1;k<l2+1;k++)
		{	d=f1*u;
			for(i=k-1;i<n-l1+k;i+=l1)
			{	ip=i+l2; i1=ar+ip; i2=ai+ip; ir=ar+i; ii=ai+i;
				t1=*i1*v+*i2*d; t2=*i2*v-*i1*d;
				*i1=*ir-t1; *i2=*ii-t2; *ir=*ir+t1; *ii=*ii+t2;
			}
			vu=v*w1-u*w2; u=v*w2+u*w1; v=vu;
	}	}
	if(key==0)
	for(i=0,ir=ar,ii=ai;i<n;i++,ir++,ii++)
		{	*ir/=n; *ii/=n;
		}
	return ;
}

////离散余弦变换DCT////////////////////////////////////////////////////////
/*=========================================================================
indata/outdata:输入数组；输出数组
n:一维数组元素个数
key:true表示正变换；false表示你变换
==========================================================================*/
void dct1(double *indata,double* outdata,int n,bool key)
{
	if (key==true)
	{
		double ck=1.0/sqrt(2)*sqrt(2.0/n);
		double m=PI/(2*n);
		//k=0
		outdata[0]=0;
		for (int i = 0; i < n; i++)
		{
			outdata[0]+=indata[i];
		}
		outdata[0]=outdata[0]*ck;
		//k≠0
		for (int j = 1; j < n; j++)
		{
			outdata[j]=0;
			for (int i = 0; i < n; i++)
			{
				outdata[j]+=indata[i]*cos((2*i+1)*j*m);
			}
			outdata[j]=outdata[j]*sqrt(2.0/n);
		}
	}else
	{
		double m=PI/(2*n);
		for (int j = 0; j < n; j++)
		{
			outdata[j]=0;
			for (int i = 0; i < n; i++)
			{
				if (i==0)
				{ 
					outdata[j]+=indata[i]/sqrt(2.0)*cos((2*j+1)*i*m);
				}else
				{
					outdata[j]+=indata[i]*cos((2*j+1)*i*m);
				}
			}
			outdata[j]=outdata[j]*sqrt(2.0/n);
		}
	}
	
}

//////////////////////////////////////////////////////////////
// To form the metrix of factors /////////////////////////////
// x,y sould are positive
void fact(double **s,double *x,double *y,double *g,double *b,
		  int mt,int n)
		  //////////// n is order of f(x,y) 
{	int i,j,l,l1,l2,lx,ly,k,k1,k2,kx,ky,kk;
	double p,f1,f2,f3,f4;
	for(i=0;i<=n;i++)
	{	l1=(i+1)*i/2+1; l2=l1+i; lx=i; ly=0; 
		for(l=l1;l<=l2;l++)
		{	for(j=0;j<=n;j++)
			{	k1=(j+1)*j/2+1;	k2=k1+j; kx=j; ky=0;
				for(k=k1;k<=k2;k++)
				{	p=0.;
					for(kk=0;kk<mt;kk++)
					{	f1=pow(x[kk],lx);
						f2=pow(y[kk],ly);
                        f3=pow(x[kk],kx);
                        f4=pow(y[kk],ky);
						p+=f1*f2*f3*f4;
					}
					s[l-1][k-1]=p;										
					kx=kx-1; ky=ky+1;
				}
			}
			lx=lx-1; ly=ly+1;
		}
	}
	for(j=0;j<=n;j++)
	{	k1=(j+1)*j/2+1; k2=k1+j; kx=j; ky=0;
		for(k=k1;k<=k2;k++)
		{	p=0.;
			for(kk=0;kk<mt;kk++)
			{	
				f1=pow(x[kk],kx);
				f2=pow(y[kk],ky);
				p+=f1*f2*g[kk];
			}
			b[k-1]=p;
			kx=kx-1; ky=ky+1;
		}
	}
	return;
}

//////////////////////////////////////////////////////////////
// To form the metrix of factors /////////////////////////////
// x,y sould are positive
void fact1(double **s,double *x,double *y,double *g,double *b,
		  int mt,int n,int xLt,int yLt)
		  //////////// m is heigher order of f(x,y) along between x and y, 
{	int i,j,l,l1,l2,lx,ly,k,k1,k2,kx,ky,kk,ks,ls;
	double p,f1,f2,f3,f4;
	ls=0;
	for(i=0;i<=n;i++)
	{	
		l1=(i+1)*i/2+1; l2=l1+i; lx=i; ly=0; 
		for(l=l1;l<=l2;l++)
		{	
			if((lx>xLt)||(ly>yLt))
			{	
				ls+=1;	lx-=1;	ly+=1;	continue; 
			}
			for(j=0;j<=n;j++)
			{	
				k1=(j+1)*j/2+1;	k2=k1+j; kx=j; ky=0;
				ks=0;
				for(k=k1;k<=k2;k++)
				{	
					if((kx>xLt)||(ky>yLt))
					{	
						ks+=1;	kx-=1;	ky+=1;	continue; 
					}
					p=0.;
					for(kk=0;kk<mt;kk++)
					{
						f1=pow(x[kk],lx);
						f2=pow(y[kk],ly);
						f3=pow(x[kk],kx);
						f4=pow(y[kk],ky);
						p+=f1*f2*f3*f4;
					}
					s[l-1-ls][k-1-ks]=p;										
					kx=kx-1; ky=ky+1;
				}
			}
			lx=lx-1; ly=ly+1;
		}
	}
	ks=0;
	for(j=0;j<=n;j++)
	{	
		k1=(j+1)*j/2+1; k2=k1+j; kx=j; ky=0;
		for(k=k1;k<=k2;k++)
		{	
			if((kx>xLt)||(ky>yLt))
			{	
				ks+=1;	kx-=1;	ky+=1;	continue; 
			}
			p=0.;
			for(kk=0;kk<mt;kk++)
			{	
				f1=pow(x[kk],kx);
				f2=pow(y[kk],ky);
				p+=f1*f2*g[kk];
			}
			b[k-1-ks]=p;
			kx=kx-1; ky=ky+1;
		}
	}
	return;
}

//////////////////////////////////////////////////////
//To calculate the multiply values of x and y/////////////
// x,y sould are positive
void tsv(double *s,double x,double y,int n)
		  //////////// n is order of f(x,y) in y,
{	int j,k,k1,k2,kx,ky;
	s[0]=(double)1;
	for(j=1;j<=n;j++)
	{	k1=(j+1)*j/2+1; k2=k1+j; kx=j; ky=0;
		for(k=k1;k<=k2;k++)
		{	s[k-1]=pow(x,kx)*pow(y,ky);
			kx=kx-1; ky=ky+1;
		}
	}
	return;
}
//////////////////////////////////////////////////////
//To calculate the multiply values of x and y/////////////
// x,y sould are positive
void tsv1(double *s,double x,double y,int n,int xLt,int yLt)
//////////// n is higher order of f(x,y) along that between x and y
{	int j,k,k1,k2,kx,ky,ks;
	ks=0;
	s[0]=(double)1;
	for(j=1;j<=n;j++)
	{	k1=(j+1)*j/2+1; k2=k1+j; kx=j; ky=0;
		for(k=k1;k<=k2;k++)
		{	if((kx>xLt)||(ky>yLt))
			{	ks+=1;	kx-=1;	ky+=1;	continue; }
			s[k-1-ks]=pow(x,kx)*pow(y,ky);
			kx=kx-1; ky=ky+1;
		}
	}
	return;
}
/////////////////////////////////
void tran(double **a,int n)
{	int i,j; double t;
	for(i=0;i<n-1;i++)
	{	for(j=i+1;j<n;j++)
		{	t=a[i][j]; a[i][j]=a[j][i]; a[j][i]=t; }
	}
	return;
}


void metrimult(double **a,double **b,int m,int n,int l,double **c)
{	
	int i,j,k; double p;
	for(i=0;i<m;i++)
	{	for(j=0;j<l;j++)
		{	p=0.;
			for(k=0;k<n;k++)
				p+=a[i][k]*b[k][j];
		   c[i][j]=p;
		}
	}
	return;
}

  ///////////////////////////

//C*********************************************************************C
//C								      C
//C	        Calculation for Singular Values and Vecters	      C
//C								      C
//C*********************************************************************C
//C [ SVD.V2 ]
//C
//C******************************************************************C
//C
//        SUBROUTINE SVD(MD,ND,M,N,A,U,V,Q,INDEX)
//C
//C    CALLS NO OTHER ROUTINES
//C    (  SINGULAR VALUE DECOMPOSITION  )
//C  FOR ALGO PROGRAM SEE WILKINSON+REINSC HANDBOOK
//C  FOR AUTOMATIC COMPUTATION VOL 2 - LINEAR ALGEBRA, PG140-144,
//C  TRANSLATED FROM ALGOL BY R.L.PARKER.
//C    THE MATRIX A(M,N) IS DECOMPOSED. SINGULAR VALUES IN Q(N),
//C  PRE-MATRIX IN U(M,N), POST-MATRIX IN V(N,N).
//C    INDEX MAY BE 1,2,3 OR 4.
//C       IF 1, FIND U,V;
//C       IF 2, FIND ONLY U;
//C       IF 3, FIND ONLY V;
//C       IF 4, FIND NEITHER.
//C	 Translated from FORTRAN to C by Chen Chao in 1997
void svd(double **a,int m,int n,double **u,double **v,double *q,int Index)
{	
	int i,j,k,l;
	double c,f,g,h,s,x,y,z;
	double *e=new double [n];
	double eps=(double)(1e-40);
	double tol=(double)(1e-60);
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			u[i][j]=a[i][j];

//C  HOUSEHOLDER REDUCTION TO BI-DIAGONAL FORM
	g=(double)0; x=(double)0;
	for(i=0;i<n;i++)
	{	
		l=i+1;	e[i]=g;  g=(double)0;
		if(i<m)
		{	
			s=(double)0;
			for(j=i;j<m;j++)
				s+=u[j][i]*u[j][i];
			if(s>=tol)
			{	
				f=u[i][i]; 
				g=-sqrt(s)*f/(fabs(f)+(double)1e-40);
				h=f*g-s;
				u[i][i]=f-g;
				if(l<n)
				{	
					for(j=l;j<n;j++)
					{
						s=(double)0;
						for(k=i;k<m;k++)
							s+=u[k][i]*u[k][j];
						f=s/h;
						for(k=i;k<m;k++)
							u[k][j]+=f*u[k][i];
					}
				}
			}
		}
//C/////////////////////////
		q[i]=g;  g=(double)0;
		if(l<n)
		{
			s=(double)0;
			for(j=l;j<n;j++)
				s+=u[i][j]*u[i][j];
			if(s>=tol)
			{
				f=u[i][i+1];
				g=-sqrt(s)*f/(fabs(f)+(double)1e-40);
				h=f*g-s;
				u[i][i+1]=f-g;
				for(j=l;j<n;j++)
					e[j]=u[i][j]/h;
				if(l<m)
				{
					for(j=l;j<m;j++)
					{
						s=(double)0;
						for(k=l;k<n;k++)
							s+=u[j][k]*u[i][k];
						for(k=l;k<n;k++)
							u[j][k]+=s*e[k];
					}
				}
			}
		}
		y=fabs(q[i])+fabs(e[i]);
		if(y>x)
			x=y;
	}
//C
//C  ACCUMULATION OF RIGHT-HAND TRANSFORMS (V)
//C
	if((Index==1)||(Index==3))
	{
		v[n-1][n-1]=(double)1;
		g=e[n-1];
		for(i=n-2;i>=0;i--)
		{
			l=i+1;
			if(g!=(double)0)
			{
				h=u[i][i+1]*g;
				for(j=l;j<n;j++)
					v[j][i]=u[i][j]/h;
				for(j=l;j<n;j++)
				{
					s=(double)0;
					for(k=l;k<n;k++)
						s+=u[i][k]*v[k][j];
					for(k=l;k<n;k++)
						v[k][j]+=s*v[k][i];
				}
			}
			for(j=l;j<n;j++)
			{
				v[j][i]=(double)0;
				v[i][j]=(double)0;
			}
			v[i][i]=(double)1;  g=e[i];
		}
	}
//C
//C  ACCUMULATION OF LEFT-HAND TRANSFORMS ( U )
//C
	if((Index==1)||(Index==2))
	{
		for(i=(n<m)?(n-1):(m-1);i>=0;i--)
		{
			l=i+1;  g=q[i];
			if(l<n)
				for(j=l;j<n;j++)
					u[i][j]=(double)0;
			if(g!=(double)0)
			{
				h=u[i][i]*g;
				if(l<n)
				{
					for(j=l;j<n;j++)
					{
						s=(double)0;
						for(k=l;k<m;k++)
							s+=u[k][i]*u[k][j];
						f=s/h;
						for(k=i;k<m;k++)
							u[k][j]+=f*u[k][i];
					}
				}
				for(j=i;j<m;j++)
					u[j][i]=u[j][i]/g;
			}
			else
			{
				for(j=i;j<m;j++)
					u[j][i]=(double)0;
			}
			u[i][i]+=(double)1;
		}
	}
//C
//C  DIAGONALIZATION OF BI-DIAGONAL FORM
//C
	eps*=x;
	for(k=n-1;k>=0;k--)
	{
label_380: /////////////C  TEST F-SPLITTING
		for(l=k;l>=0;l--)	   
		{
			if(fabs(e[l])<=eps)
				goto label_440;
			if(l==0)
				continue;
			if(fabs(q[l-1])<=eps)
				goto label_410;
		}

label_410:  /////  CANCELLATION OF E(L), IF L .GT. 1
		c=(double)0; 
		s=(double)1;
		for(i=l;i<=k;i++)
		{
			f=s*e[i];  e[i]*=c;
			if(fabs(f)<=eps)
				goto label_440;
			g=q[i];  q[i]=sqrt(f*f+g*g);  
			h=q[i];  c=g/h;  s=-f/h;  
			if((Index==1)||(Index==2))
			{
				for(j=0;j<m;j++)
				{
					y=u[j][l-1];  z=u[j][i];
					u[j][l-1]=y*c+z*s;
					u[j][i]=-y*s+z*c;
				}
			}
		}

label_440: ////  TEST F-CONVERGENCE
		z=q[k];
		if(l==k)
			goto label_480;
///////////  SHIFT FROM BOTTOM 2 X 2 MINOR			
		x=q[l];    
		y=q[k-1];
		g=e[k-1];
		h=e[k];
		f=(double)0.5*((y-z)*(y+z)+(g-h)*(g+h))/h/y;
		g=sqrt(f*f+(double)1);
		f=((x-z)*(x+z)+h*(y/(f+g*f/(fabs(f)+(double)1e-40))-h))/x;
/////////////  NEXT Q-R TRANSFORMATION
        c=(double)1;   
		s=(double)1;
		for(i=l+1;i<=k;i++)
		{
			g=e[i];  y=q[i];  h=s*g;  g=c*g;  
			z=sqrt(f*f+h*h);
			e[i-1]=z;  c=f/z;  s=h/z;  
			f=x*c+g*s;  g=-x*s+g*c;
			h=y*s;  y=y*c;  
			if((Index==1)||(Index==3))
			{
				for(j=0;j<n;j++)
				{
					x=v[j][i-1];  z=v[j][i];
					v[j][i-1]=x*c+z*s;
					v[j][i]=-x*s+z*c;
				}
			}
			z=sqrt(f*f+h*h);
			q[i-1]=z;  c=f/z;  s=h/z;  
			f=c*g+s*y;  x=-s*g+c*y;
			if((Index==1)||(Index==2))
			{
				for(j=0;j<m;j++)
				{
					y=u[j][i-1];  z=u[j][i];  
					u[j][i-1]=y*c+z*s;
					u[j][i]=-y*s+z*c;
				}
			}
		}
		e[l]=(double)0;
		e[k]=f;
		q[k]=x;
		goto label_380;

label_480: 	///////////C  CONVERGENCE///////////////
		if(z<(double)0)
		{
			q[k]=-z;
			if((Index==1)||(Index==3))
			{
				for(j=0;j<n;j++)
					v[j][k]=-v[j][k];
			}
		}
	}
	delete e;
	return ;
}

/*=======================局部多项式扩边==========================
输入参数:
		indatax/indatay:输入待扩边的x坐标y坐标
		outdatax/outdatay:扩边以后的x坐标y坐标
		indatanum/outdatanum:待扩边的数据个数；扩边后的数据个数；
		N:多项式阶数
==============================================================*/
void Extension_Polynomial_1D(double* indatax,double* indatay,double* outdatax,double* outdatay,const int indatanum,const int outdatanum,const int N)
{
	//const int N=2;
	int coefficientNum=N+1;				//N阶多项式扩边
	int extentNum=outdatanum-indatanum;
	int leftnum=extentNum/2;
	int rightnum=extentNum-leftnum;
	double* Left=new double[leftnum],*right=new double[rightnum];
	double dx=fabs(indatax[1]-indatax[0]);

	//1. 求扩边后的横坐标
	for (int i = 0; i < leftnum; i++)
	{
		outdatax[i]=indatax[0]-(leftnum-i)*dx;
	}
	for (int i = 0; i < indatanum; i++)
	{
		outdatax[leftnum+i]=indatax[i];
	}
	for (int i = 0; i < rightnum; i++)
	{
		outdatax[leftnum+indatanum+i]=indatax[indatanum-1]+(i+1)*dx;
	//扩边端点值
	Left[0]=(indatay[0]+indatay[indatanum-1])/2.0;
	right[rightnum-1]=Left[0];
	}
//2. 先求左边
	
	//求解多项式系数:AX=b,求X
	double* b=new double[coefficientNum],*X=new double[coefficientNum];
	double** A=new double* [coefficientNum];
	for (int i = 0; i < coefficientNum; i++)
	{
		A[i]=new double[coefficientNum];
	}
	for (int i = 1; i < coefficientNum; i++)
	{
		for (int j = 0; j < coefficientNum; j++)
		{
			A[i][j]=pow(outdatax[leftnum+i-1],j);
		}
		b[i]=indatay[i-1];
	}b[0]=Left[0];
	for (int i = 0; i < coefficientNum; i++)
	{
		A[0][i]=pow(outdatax[0],i);
	}
	//利用高斯消元法求X
	gjdn  c(coefficientNum, 1);
	c.InPut(A,b);
	c.gauss_jordan ();  //执行Gauss-Jordan消去法
	c.OutPut(X);
	//利用上面求得的X计算扩边值
	for (int i = 1; i < leftnum; i++)
	{
		Left[i]=0;
		for (int j = 0; j < coefficientNum; j++)
		{
			Left[i]+=X[j]*pow(outdatax[i],j);
		}
	}

//3. 求右边扩边
	for (int i = 0; i < coefficientNum-1; i++)
	{
		for (int j = 0; j < coefficientNum; j++)
		{
			A[i][j]=pow(outdatax[indatanum+leftnum-coefficientNum+i],j);
		}
		b[i]=indatay[indatanum-coefficientNum+i];
	}b[coefficientNum-1]=right[rightnum-1];
	for (int i = 0; i < coefficientNum; i++)
	{
		A[coefficientNum-1][i]=pow(outdatax[outdatanum-1],i);
	}
	c.InPut(A,b);
	c.gauss_jordan ();  //执行Gauss-Jordan消去法
	c.OutPut(X);
	//利用上面求得的X计算扩边值
	for (int i = 0; i < rightnum-1; i++)
	{
		right[i]=0;
		for (int j = 0; j < coefficientNum; j++)
		{
			right[i]+=X[j]*pow(outdatax[i+leftnum+indatanum],j);
		}
	}
	
	//整合扩边值
	for (int i = 0; i < leftnum; i++)
	{
		outdatay[i]=Left[i];
	}
	for (int i = 0; i < indatanum; i++)
	{
		outdatay[leftnum+i]=indatay[i];
	}
	for (int i = 0; i < rightnum; i++)
	{
		outdatay[leftnum+indatanum+i]=right[i];
	}

	//销毁指针
	delete[] right,Left,A,X,b;
}
/*=======================局部多项式扩边==========================
输入参数:
		outdatax/outdatay:扩边以后的x坐标y坐标
		dx:x间距
		indatanum/outdatanum:待扩边的数据个数；扩边后的数据个数；
		N:多项式阶数
==============================================================*/
void Extension_Polynomial_1D(double* indatay,double* outdatay,const double dx,const int indatanum,const int outdatanum,const int N)
{
	//const int N=2;
	int coefficientNum=N+1;				//N阶多项式扩边
	int extentNum=outdatanum-indatanum;
	int leftnum=extentNum/2;
	int rightnum=extentNum-leftnum;
	double* Left=new double[leftnum],*right=new double[rightnum];
	//double dx=fabs(indatax[1]-indatax[0]);
	double* indatax=new double[indatanum],*outdatax=new double[outdatanum];
	for (int i = 0; i < indatanum; i++)
	{
		indatax[i]=i*dx;
	}
	//1. 求扩边后的横坐标
	for (int i = 0; i < leftnum; i++)
	{
		outdatax[i]=indatax[0]-(leftnum-i)*dx;
	}
	for (int i = 0; i < indatanum; i++)
	{
		outdatax[leftnum+i]=indatax[i];
	}
	for (int i = 0; i < rightnum; i++)
	{
		outdatax[leftnum+indatanum+i]=indatax[indatanum-1]+(i+1)*dx;
	//扩边端点值
	Left[0]=(indatay[0]+indatay[indatanum-1])/2.0;
	//Left[0]=0;
	right[rightnum-1]=Left[0];
	}
//2. 先求左边
	
	//求解多项式系数:AX=b,求X
	double* b=new double[coefficientNum],*X=new double[coefficientNum];
	double** A=new double* [coefficientNum];
	for (int i = 0; i < coefficientNum; i++)
	{
		A[i]=new double[coefficientNum];
	}
	for (int i = 1; i < coefficientNum; i++)
	{
		for (int j = 0; j < coefficientNum; j++)
		{
			A[i][j]=pow(outdatax[leftnum+i-1],j);
		}
		b[i]=indatay[i-1];
	}b[0]=Left[0];
	for (int i = 0; i < coefficientNum; i++)
	{
		A[0][i]=pow(outdatax[0],i);
	}
	//利用高斯消元法求X
	gjdn  c(coefficientNum, 1);
	c.InPut(A,b);
	c.gauss_jordan ();  //执行Gauss-Jordan消去法
	c.OutPut(X);
	//利用上面求得的X计算扩边值
	for (int i = 1; i < leftnum; i++)
	{
		Left[i]=0;
		for (int j = 0; j < coefficientNum; j++)
		{
			Left[i]+=X[j]*pow(outdatax[i],j);
		}
	}

//3. 求右边扩边
	for (int i = 0; i < coefficientNum-1; i++)
	{
		for (int j = 0; j < coefficientNum; j++)
		{
			A[i][j]=pow(outdatax[indatanum+leftnum-coefficientNum+i],j);
		}
		b[i]=indatay[indatanum-coefficientNum+i];
	}b[coefficientNum-1]=right[rightnum-1];
	for (int i = 0; i < coefficientNum; i++)
	{
		A[coefficientNum-1][i]=pow(outdatax[outdatanum-1],i);
	}
	c.InPut(A,b);
	c.gauss_jordan ();  //执行Gauss-Jordan消去法
	c.OutPut(X);
	//利用上面求得的X计算扩边值
	for (int i = 0; i < rightnum-1; i++)
	{
		right[i]=0;
		for (int j = 0; j < coefficientNum; j++)
		{
			right[i]+=X[j]*pow(outdatax[i+leftnum+indatanum],j);
		}
	}
	
	//整合扩边值
	for (int i = 0; i < leftnum; i++)
	{
		outdatay[i]=Left[i];
	}
	for (int i = 0; i < indatanum; i++)
	{
		outdatay[leftnum+i]=indatay[i];
	}
	for (int i = 0; i < rightnum; i++)
	{
		outdatay[leftnum+indatanum+i]=right[i];
	}

	//销毁指针
	delete[] right,Left,A,X,b,indatax,outdatax;
}


//-------------快速排序算法(来自百度百科)--------------------------
//输入数组arr
//开始序号low
//结束序号high
//----------------------------------------------------------------
void Qsort(double a[],int low,int high)
{
	if(low>=high)
	{
	return;
	}
	int first=low;
	int last=high;
	double key=a[first];//用字表的第一个记录作为枢轴
	while(first<last)
	{
	while(first<last&&a[last]>=key)--last;
	a[first]=a[last];//将比第一个小的移到低端
	while(first<last&&a[first]<=key)++first;
	a[last]=a[first];//将比第一个大的移到高端
	}
	a[first]=key;//枢轴记录到位
	Qsort(a,low,first-1);
	Qsort(a,last+1,high);
}