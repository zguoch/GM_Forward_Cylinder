

//#include "stdafx.h"
#include "proj.h"
#include "Filter.h"

bool DoFilter(double **raw,double **a,int mRow,int mClm,
					  double dx,double dy) 
{
	/*c=raw;*/ 
	int m=mRow; int n=mClm;
	int me,ne,m1,n1,np,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		if(m==1)
		{   me=0; m1=1; np=0;	}
		else
		{   me=4; m1=(int)(pow(2,me)+.001);
			while(m>=(m1-12))
			{   me++; m1=(int)(pow(2,me)+.001);	}
		}
		ne=4; n1=(int)(pow(2,ne)+.001);
		while(n>=(n1-12))
		{   ne++; n1=(int)(pow(2,ne)+.001);
		}
		np=(int)((m1-m)/2+.001); nq=(int)((n1-n)/2+.001);

		int mj=m1; int nj=n1;
		double **br,**bi;
		br=new double* [mj];	bi=new double* [mj];	
		br[0]=new double [mj*nj]; bi[0]=new double [mj*nj]; 
		for(int i=1;i<mj;i++)
		{	br[i]=br[i-1]+nj; bi[i]=bi[i-1]+nj; }
//////////////////////////////////////////////////////////////
		for(int i=0;i<m;i++)
			for(int j=0;j<n;j++)
				br[i][j]=raw[i][j];
 		filter(br,bi,m,n,m1,n1,me,ne,np,nq,dx,dy,
			0/*mContUpDown 向上延拓：-1；下延：1*/,
			0/*m_DerivAngleNum 求导方向角*/,
			0/*mFiltLPCtrl 低通滤波*/,
			0/*mFiltHPCtrl 高通滤波*/,
			0/*mFiltBPCtrl 带通滤波*/,
			0/*mContCtrl 延拓开关*/,
			0/*mDerivZCtrl 垂向一次导数*/,
			0/*mDerivZZCtrl 垂向二次导数*/,
			1/*mDerivHCtrl 水平一次导数*/,
			0,//分量换算开关
			3,//分量换算类型：1,2,3分别表示由Ta换算x，y，z分量
			3000/*m_ContNum 延拓高度*/,
			0/*m_FiltStart 开始波长*/,
			0/*m_FiltStop 结束波长*/);	
  		for(int i=0;i<m;i++)
			for(int j=0;j<n;j++)
				a[i][j]=br[i][j];
		delete br;
		delete bi;
		return TRUE;
}
void filter(double **zr,double **zi,int m,int n,
				 int m1,int n1,int me,int ne,
				 int np,int nq,double dx,double dy,
				 int iP1,int iP2,int iP3,int iP4,int iP5,int iP6,
				 int iP7,int iP8,int iP9,int iP10,int iP11,double fP1,double fP2,double fP3)	
{	
	int i,j,i1,j1; 
	double u1,u2,v1,v2,x11,x12,x21,x22,y11,y12,y21,y22,w11,w12,w21,w22;
	double co,si,em,en,emn,h,df,fb1,fb2,fe1,fe2,pjg,tp,temp;
	double AC11,AC12,AC21,AC22;
	double p11,p12,p21,p22,q11,q12,q21,q22;		//分量换算时使用
	//===============
	int mContUpDown=iP1;
	int m_DerivAngleNum=iP2;
	int mFiltLPCtrl=iP3;
	int mFiltHPCtrl=iP4;
	int mFiltBPCtrl=iP5;
	int mContCtrl=iP6;
	int mDerivZCtrl=iP7;
	int mDerivZZCtrl=iP8;
	int mDerivHCtrl=iP9;
	int mComponent=iP10;
	int mComponentStyle=iP11;
	double m_ContNum=fP1;
	double m_FiltStart=fP2;
	double m_FiltStop=fP3;	
	//===============
	h=(double)(mContUpDown*m_ContNum);
	co=(double)cos(PI*m_DerivAngleNum/180.); 
	si=(double)sin(PI*m_DerivAngleNum/180.);
	em=dx*(double)(m1-1); en=dy*(double)(n1-1); tp=(double)(2*PI);
	emn=(double)(0.5*sqrt(m1*m1/(em*em)+n1*n1/(en*en)));
	df=(double)(5.*sqrt(1./(em*em)+1./(en*en)));
	if(m_FiltStart>m_FiltStop)
	{	
		fb1=1/(double)(m_FiltStart+1e-6)-df; 
		fb2=1/(double)(m_FiltStart+1e-6); 
		fe1=1/(double)(m_FiltStop+1e-6); 
		fe2=1/(double)(m_FiltStop+1e-6)+df; 
	}
	else
	{	
		fb1=1/(double)(m_FiltStop+1e-6)-df; 
		fb2=1/(double)(m_FiltStop+1e-6); 
		fe1=1/(double)(m_FiltStart+1e-6); 
		fe2=1/(double)(m_FiltStart+1e-6)+df; 
	}
	pjg=prep(zr,zi,m,n,m1,n1,np,nq);
	fft(zr,zi,m1,n1,me,ne,0);
	if(m>1)
	{
		for(i=0;i<(m1/2);i++)
		{	
			u1=(double)i/em+1E-10; u2=(double)(i-m1/2)/em+1E-10; i1=m1/2+i;
			for(j=0;j<(n1/2);j++)
			{	
				v1=(double)j/en+1E-10; v2=(double)(j-n1/2)/en+1E-10; j1=n1/2+j;
				x11=tp*(u1*co+v1*si); x12=tp*(u1*co+v2*si);
				x21=tp*(u2*co+v1*si); x22=tp*(u2*co+v2*si);
				y11=(double)sqrt(u1*u1+v1*v1); y12=(double)sqrt(u1*u1+v2*v2);
				y21=(double)sqrt(u2*u2+v1*v1); y22=(double)sqrt(u2*u2+v2*v2);
				w11=tp*y11; w12=tp*y12;
				w21=tp*y21; w22=tp*y22;
				AC11=AC12=AC21=AC22=(double)1;	
				if(mFiltLPCtrl==1)
				{
					AC11*=(double)(0.5*(1+cos(PI*y11/emn))); AC12*=(double)(0.5*(1+cos(PI*y12/emn)));
					AC21*=(double)(0.5*(1+cos(PI*y21/emn))); AC22*=(double)(0.5*(1+cos(PI*y22/emn)));
				}
				if(mFiltHPCtrl==1)
				{
					AC11*=(double)(0.5*(1-cos(PI*y11/emn))); AC12*=(double)(0.5*(1-cos(PI*y12/emn)));
					AC21*=(double)(0.5*(1-cos(PI*y21/emn))); AC22*=(double)(0.5*(1-cos(PI*y22/emn)));
				}
				if(mFiltBPCtrl==1)
				{
					if(y11<=fb1)
						AC11*=(double)0;
					else if((y11>fb1)&&(y11<=fb2))  
						AC11*=(double)(0.5*(1-cos(PI*(y11-fb2)/df)));
					else if((y11>fb2)&&(y11<=fe1))  
						AC11*=(double)1;
					else if((y11>fe1)&&(y11<=fe2))  
						AC11*=(double)(0.5*(1+cos(PI*(y11-fe2)/df)));
					else                            
						AC11*=(double)0;				
					if(y12<=fb1)                    
						AC12*=(double)0;
					else if((y12>fb1)&&(y12<=fb2))  
						AC12*=(double)(0.5*(1-cos(PI*(y12-fb2)/df)));
					else if((y12>fb2)&&(y12<=fe1))  
						AC12*=(double)1;
					else if((y12>fe1)&&(y12<=fe2))  
						AC12*=(double)(0.5*(1+cos(PI*(y12-fe2)/df)));
					else                            
						AC12*=(double)0;				
					if(y21<=fb1)                    
						AC21*=(double)0;
					else if((y21>fb1)&&(y21<=fb2))  
						AC21*=(double)(0.5*(1-cos(PI*(y21-fb2)/df)));
					else if((y21>fb2)&&(y21<=fe1))  
						AC21*=(double)1;
					else if((y21>fe1)&&(y21<=fe2))  
						AC21*=(double)(0.5*(1+cos(PI*(y21-fe2)/df)));
					else                            
						AC21*=(double)0;				
					if(y22<=fb1)                    
						AC22*=(double)0;
					else if((y22>fb1)&&(y22<=fb2))  
						AC22*=(double)(0.5*(1-cos(PI*(y22-fb2)/df)));
					else if((y22>fb2)&&(y22<=fe1))  
						AC22*=(double)1;
					else if((y22>fe1)&&(y22<=fe2))  
						AC22*=(double)(0.5*(1+cos(PI*(y22-fe2)/df)));
					else                            
						AC22*=(double)0;				
				}
				if(mContCtrl==1)
				{ 
					AC11*=(double)exp(w11*h); AC12*=(double)exp(w12*h);
					AC21*=(double)exp(w21*h); AC22*=(double)exp(w22*h);
				}
				if(mDerivZCtrl==1)
				{	
					AC11*=w11; AC12*=w12; AC21*=w21; AC22*=w22;
				}
				if(mDerivZZCtrl==1)
				{	
					AC11*=w11*w11; AC12*=w12*w12; AC21*=w21*w21; AC22*=w22*w22;
				}
				zr[i][j]=zr[i][j]*AC11; zi[i][j]=zi[i][j]*AC11;
				zr[i][j1]=zr[i][j1]*AC12; zi[i][j1]=zi[i][j1]*AC12;
				zr[i1][j]=zr[i1][j]*AC21; zi[i1][j]=zi[i1][j]*AC21;
				zr[i1][j1]=zr[i1][j1]*AC22; zi[i1][j1]=zi[i1][j1]*AC22;
				if(mDerivHCtrl==1)
				{	
					//========================水平一次导数============================
					temp=zr[i][j]; zr[i][j]=-zi[i][j]*x11; zi[i][j]=temp*x11;
					temp=zr[i][j1]; zr[i][j1]=-zi[i][j1]*x12; zi[i][j1]=temp*x12;
					temp=zr[i1][j]; zr[i1][j]=-zi[i1][j]*x21; zi[i1][j]=temp*x21;
					temp=zr[i1][j1]; zr[i1][j1]=-zi[i1][j1]*x22; zi[i1][j1]=temp*x22;
					//====================================================================
					
				}
				if (mComponent==1)
				{
					double angleI=-90.0/180.0*PI,angleD=10.0/180*PI;
					double aerfa0=cos(angleI)*cos(angleD),beta0=cos(angleI)*sin(angleD),gama0=sin(angleI);
					double kx1=aerfa0*u1,kx2=aerfa0*u2;
					double ky1=beta0*v1,ky2=beta0*v2;
					double kz11=gama0*sqrt(u1*u1+v1*v1),kz22=gama0*sqrt(u2*u2+v2*v2),kz12=gama0*sqrt(u1*u1+v2*v2),kz21=gama0*sqrt(u2*u2+v1*v1);
					switch (mComponent)
					{
					case 1:			//X分量
						q11=kx1*kz11/aerfa0/(pow(kx1+ky1,2.0)+pow(kz11,2.0));
						q12=kx1*kz12/aerfa0/(pow(kx1+ky2,2.0)+pow(kz12,2.0));
						q21=kx2*kz21/aerfa0/(pow(kx2+ky1,2.0)+pow(kz21,2.0));
						q22=kx2*kz22/aerfa0/(pow(kx2+ky2,2.0)+pow(kz22,2.0));
						p11=kx1*(kx1+ky1)/aerfa0/(pow(kx1+ky1,2.0)+pow(kz11,2.0));
						p12=kx1*(kx1*ky2)/aerfa0/(pow(kx1+ky2,2.0)+pow(kz12,2.0));
						p21=kx2*(kx2*ky1)/aerfa0/(pow(kx2+ky1,2.0)+pow(kz21,2.0));
						p22=kx2*(kx2*ky2)/aerfa0/(pow(kx2+ky2,2.0)+pow(kz22,2.0));
						break;
					case 2:			//Y分量
							break;
					case 3:				//Z分量
						p11=kz11*kz11/gama0/(pow(kx1+ky1,2.0)+pow(kz11,2.0));
						p12=kz12*kz12/gama0/(pow(kx1+ky2,2.0)+pow(kz12,2.0));
						p21=kz21*kz21/gama0/(pow(kx2+ky1,2.0)+pow(kz21,2.0));
						p22=kz22*kz22/gama0/(pow(kx2+ky2,2.0)+pow(kz22,2.0));
						q11=-kz11*(kx1+ky1)/gama0/(pow(kx1+ky1,2.0)+pow(kz11,2.0));
						q12=-kz12*(kx1*ky2)/gama0/(pow(kx1+ky2,2.0)+pow(kz12,2.0));
						q21=-kz21*(kx2*ky1)/gama0/(pow(kx2+ky1,2.0)+pow(kz21,2.0));
						q22=-kz22*(kx2*ky2)/gama0/(pow(kx2+ky2,2.0)+pow(kz22,2.0));
						break;
					default:
						break;
					}
					temp=zr[i][j]; zr[i][j]=(temp*p11-zi[i][j]*q11); zi[i][j]=temp*q11+zi[i][j]*p11;
					temp=zr[i][j1]; zr[i][j1]=(temp*p12-zi[i][j1]*q12); zi[i][j1]=temp*q12+zi[i][j1]*p12;
					temp=zr[i1][j]; zr[i1][j]=(temp*p21-zi[i1][j]*q21); zi[i1][j]=temp*q21+zi[i1][j]*p21;
					temp=zr[i1][j1]; zr[i1][j1]=(temp*p22-zi[i1][j1]*q22); zi[i1][j1]=temp*q22+zi[i1][j1]*p22;
					////====================================================================
				}
			}
		}	
	}
	else
	{	
		emn=n1/(double)(2*en);	df=(double)5./en;
		for(j=0;j<(n1/2);j++)
		{	
			v1=j/en; v2=(j-n1/2)/en; j1=n1/2+j;
			y11=tp*v1; y22=tp*v2;
			w11=tp*(double)fabs(v1); w22=tp*(double)(v2);
			AC11=AC22=(double)1;
			if(mFiltLPCtrl==1)
			{
				AC11*=(double)(0.5*(1+cos(PI*y11/emn)));  AC22*=(double)(0.5*(1+cos(PI*y22/emn)));
			}
			if(mFiltHPCtrl==1)
			{
				AC11*=(double)(0.5*(1-cos(PI*y11/emn)));  AC22*=(double)(0.5*(1-cos(PI*y22/emn)));
			}
			if(mFiltBPCtrl==1)
			{
				if(y11<=fb1)                    AC11*=(double)0;
				else if((y11>fb1)&&(y11<=fb2))  AC11*=(double)(0.5*(1-cos(PI*(y11-fb2)/df)));
				else if((y11>fb2)&&(y11<=fe1))  AC11*=(double)1;
				else if((y11>fe1)&&(y11<=fe2))  AC11*=(double)(0.5*(1+cos(PI*(y11-fe2)/df)));
				else                            AC11*=(double)0;				
				if(y22<=fb1)                    AC22*=(double)0;
				else if((y22>fb1)&&(y22<=fb2))  AC22*=(double)(0.5*(1-cos(PI*(y22-fb2)/df)));
				else if((y22>fb2)&&(y22<=fe1))  AC22*=(double)1;
				else if((y22>fe1)&&(y22<=fe2))  AC22*=(double)(0.5*(1+cos(PI*(y22-fe2)/df)));
				else                            AC22*=(double)0;				
			}
			if(mContCtrl==1)
			{ 	
				AC11*=(double)exp(w11*h); AC22*=(double)exp(w22*h);	
			}
			if(mDerivZCtrl==1)
			{	
				AC11*=w11; AC22*=w22;
			}
			if(mDerivZZCtrl==1)
			{	
				AC11*=w11*w11; AC22*=w22*w22;
			}
			zr[0][j]=zr[0][j]*AC11; zi[0][j]=zi[0][j]*AC11;
			zr[0][j1]=zr[0][j1]*AC22; zi[0][j1]=zi[0][j1]*AC22;
			if(mDerivHCtrl==1)
			{	
				temp=zr[0][j]; zr[0][j]=-zi[0][j]*y11; zi[0][j]=temp*y11;
				temp=zr[0][j1]; zr[0][j1]=-zi[0][j1]*y22; zi[0][j1]=temp*y22;
			}
		}
	}
	fft(zr,zi,m1,n1,me,ne,1);
	int FAC;
	FAC=(1-mDerivZCtrl)*(1-mDerivZZCtrl)*(1-mDerivHCtrl);
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			zr[i][j]=zr[i+np][j+nq]+FAC*pjg;
	return ;
}