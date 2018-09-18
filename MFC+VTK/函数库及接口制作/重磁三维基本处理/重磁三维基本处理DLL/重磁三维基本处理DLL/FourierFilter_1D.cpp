
//#include "stdafx.h"
#include "proj.h"
#include "FourierFilter_1D.h"


/*=============================================
函数功能：频率域求一位位场数据向上延拓
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距
	rph：延拓高度
=============================================*/
bool UpWardContinuation(double* OriginalData,double* TransData,const int datanum,const double dx,const double rph)
{

	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=-fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//====================================
		for(int j=0;j<(n1/2);j++)
		{	
			v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
			w11=tp*(double)fabs(v1); w22=tp*(double)fabs((v2));
			AC11=AC22=(double)1;
			AC11*=(double)exp(w11*h); AC22*=(double)exp(w22*h);	
			zr[j]=zr[j]*AC11; zi[j]=zi[j]*AC11;
			zr[j1]=zr[j1]*AC22; zi[j1]=zi[j1]*AC22;
		}
	//傅里叶逆变换
	fft1(zr,zi,ne,n1,1);
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=zr[j+nq];

	//销毁指针
	delete zr,zi;
	return true;
}


/*=============================================
函数功能：频率域求一位位场数据向上延拓(利用离散余弦变换）
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距
	rph：延拓高度
=============================================*/
bool UpWardContinuation_Dct(double* OriginalData,double* TransData,const int datanum,const double dx,const double rph)
{

	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=-fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//离散余弦变换
	dct1(zr,zi,n1);
	//====================================
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		w11=tp*(double)fabs(v1); w22=tp*(double)fabs((v2));
		AC11=AC22=(double)1;
		AC11*=(double)exp(w11*h); AC22*=(double)exp(w22*h);	
		zi[j]=zi[j]*AC11/v1;
		zi[j1]=zi[j1]*AC22/v2;
	}
	//傅里叶逆变换
	dct1(zi,zr,n1,false);
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=zr[j+nq];

	//销毁指针
	delete zr,zi;
	return true;
}

/*===============================================================================
函数功能：频率域求一位位场数据归一化总梯度（主要是重力，磁场数据有待研究）
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum:谐波数
	rph：延拓高度
===============================================================================*/
bool NormalFullGradient(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	double* Vxzr=new double[n1],*Vxzi=new double[n1],*Vzzr=new double[n1],*Vzzi=new double[n1];
	double* Vzr=new double[n1],*Vzi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//====================================
	double temp;
	double a=1;
	double x11,x22;
	double P=5;
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=tp*(double)v1; x22=tp*(double)v2;
		w11=fabs(x11); w22=fabs(x22);
		AC11=AC22=(double)1;
		//AC11*=(double)pow(exp(-w11*h),P)*pow((4-exp(-w11*h)-2*cos(w11*h)),P); AC22*=(double)pow(exp(-w22*h),P)*pow((4-exp(-w22*h)-2*cos(w22*h)),P);//差分法向下延拓回返
		AC11*=(double)exp(-w11*h)/(exp(-2*w11*h)+TRP); AC22*=(double)exp(-w22*h)/(exp(-2*w22*h)+TRP);	
		////计算x方向导数加延拓
		Vxzr[j]=-zi[j]*AC11*x11;Vxzi[j]=zr[j]*AC11*x11;
		Vxzr[j1]=-zi[j1]*AC22*x22;Vxzi[j1]=zr[j1]*AC22*x22;
		////计算z方向导数加延拓
		Vzzr[j]=zr[j]*AC11*w11;Vzzi[j]=zi[j]*AC11*w11;
		Vzzr[j1]=zr[j1]*AC22*w22;Vzzi[j1]=zi[j1]*AC22*w22;
		//只计算延拓
		Vzr[j]=-zi[j]*AC11;Vxzi[j]=zr[j]*AC11;
		Vzi[j1]=-zi[j1]*AC22;Vxzi[j1]=zr[j1]*AC22;
	}
	//傅里叶逆变换
	fft1(Vxzr,Vxzi,ne,n1,1);
	fft1(Vzzr,Vzzi,ne,n1,1);
	fft1(Vzr,Vzi,ne,n1,1);
	//..缩边
	for(int j=0;j<n;j++)
	{
		Vxzr[j]=/*Vzzr[j+nq]*1E4*/Vxzr[j+nq];
		Vzzr[j]=/*Vzzr[j+nq]*1E4*/Vzzr[j+nq];
	}
	//计算归一化总梯度
	double Gxz=0;
	double *GH=new double[n1],*GH0=new double[n1];
	for (int i = 0; i < n1; i++)
	{
		GH0[i]=sqrt(Vxzr[i]*Vxzr[i]+Vzzr[i]*Vzzr[i]);
		//GH0[i]=Vzr[i];
	}
	//===========窗口归一化===============================
	/*int WindowSize=datanum/2;
	int HalfWindowSize=WindowSize/2;
	double GHmax;
	Gxz=0;
	GHmax=GH0[0];
	for (int j = 0; j < WindowSize; j++)
	{
		if (GHmax<GH0[j])
		{
			GHmax=GH0[j];
		}	
	}
	for (int i = 0; i < HalfWindowSize; i++)
	{
		GH[i]=GH0[i]/GHmax;
	}

	for (int i = 0; i < datanum-WindowSize; i++)
	{
		GHmax=GH0[i];
		for (int j = 0; j < WindowSize; j++)
		{
			if (GHmax<GH0[i+j])
			{
				GHmax=GH0[i+j];
			}	
		}
		GH[i+HalfWindowSize]=GH0[i+HalfWindowSize]/GHmax;
	}

	for (int j = datanum-HalfWindowSize; j < datanum; j++)
	{
		GH[j]=GH0[j]/GHmax;
	}*/
	//===================================================
	//=========算数平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz+=GH0[i];
	//}Gxz=Gxz/(n);//均值
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================

	//=========调和平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz+=1.0/GH0[i];
	//}Gxz=n/Gxz;//
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================

	//=========几何平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz*=GH0[i];
	//}Gxz=pow(Gxz,1.0/n);//
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================
	switch (NormalizeType)
	{
	case 0:			//算术平均归一
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=GH0[i];
			}Gxz=Gxz/(n);//均值
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}

		}
		break;
	case 1:		//中值归一
		{
			double* GH0sort=new double[datanum];
			for (int k = 0; k < datanum; k++)
			{
				GH0sort[k]=GH0[k];
			}
			Qsort(GH0sort,0,datanum-1);//排序
			if (datanum%2==0)			//取中值
			{
				Gxz=(GH0sort[datanum/2]+GH0sort[datanum/2-1])/2.0;
			}else
			{
				Gxz=GH0sort[datanum/2];
			}
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
			delete GH0sort;
		}
		break;
	case 2:
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=1.0/GH0[i];
			}Gxz=n/Gxz;//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
	case 3:
		{
			for (int i = 0; i < n; i++)
			{
				Gxz*=GH0[i];
			}Gxz=pow(Gxz,1.0/n);//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
	default:
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=1.0/GH0[i];
			}Gxz=n/Gxz;//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
		break;
	}
	//..变换后数据赋值
	for(int j=0;j<n;j++)
	{
		TransData[j]=/*Vxzr[j]*1E4*/GH[j];
	}
	//销毁指针
	delete zr,zi,Vzzi,Vzzr,Vxzi,Vxzr,GH,GH0,Vzr,Vzi;
	return true;
}
/*=============================================
函数功能：频率域求一位位场数据归一化总梯度（主要是重力，磁场数据有待研究）
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组(列数与originaldata相同，行数为断面层数）
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum_max:谐波数最大值（最小值默认为1）
	h2/dh/h1：断面下层z坐标；相邻两层间距；顶层z坐标
=============================================*/
bool NormalFullGradient(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const int hNum,const double dh,const double h1,int NormalizeType)
{
	double *GHdata=new double[datanum];
	double* TransData2=new double[datanum];
	double rph=0;
	double trptemp;
	for (int j = 0; j < hNum; j++)
	{
		trptemp=TRP/**(j+1)*/;
		rph=h1+dh*j;
		NormalFullGradient(OriginalData,GHdata,datanum,dx,trptemp,rph,NormalizeType);
		//================滑动窗口归一化==============================
		int WinSize=20;
		
	//Separation_EntropyFilter(GHdata,TransData2,datanum,dx,WinSize);
	//=====================================================================
		for (int i = 0; i < datanum; i++)
		{
			TransData[j*datanum+i]=GHdata[i];
		}
	}
	delete GHdata,TransData2;
	return true;
}

/*=============================================
函数功能：泰勒级数迭代法求取归一化总梯度
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组(列数与originaldata相同，行数为断面层数）
	datanum：数据个数
	dx：点距;
	N/L/K：谐波数；展开项数；迭代次数
	h2/dh/h1：断面下层z坐标；相邻两层间距；顶层z坐标
=============================================*/
bool NormalFullGradient_TSI(double* OriginalData,double* TransData,const int datanum,const double dx,const int N,const int L,const int K,const double rph)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	double* Vxzr=new double[n1],*Vxzi=new double[n1],*Vzzr=new double[n1],*Vzzi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//====================================
	double pn1,pn2;
	double temp;
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		pn1=pow(sin(v1*en*PI/N)/(v1*en*PI/N),2);pn2=pow(sin(v2*en*PI/N)/(v2*en*PI/N),2);
		w11=tp*(double)fabs(v1); w22=tp*(double)fabs((v2));
		AC11=AC22=(double)1;
		AC11*=(double)exp(w11*pn1*h)*pow(1-exp(-w11*pn1*h),K+1); AC22*=(double)exp(w22*pn2*h)*pow(1-exp(-w22*pn2*h),K+1);	
		////计算x方向导数加延拓
		temp=zr[j];Vxzr[j]=-zi[j]*AC11*w11;Vxzi[j]=temp*AC11*w22;
		temp=zr[j1];Vxzr[j1]=-zi[j1]*AC22*w11;Vxzi[j1]=temp*AC22*w22;
		////计算z方向导数加延拓
		Vzzr[j]=zr[j]*AC11*w11;Vzzi[j]=zi[j]*AC11*w22;
		Vzzr[j1]=zr[j1]*AC22*w11;Vzzi[j1]=zi[j1]*AC22*w22;
	}
	//傅里叶逆变换
	fft1(Vxzr,Vxzi,ne,n1,1);
	fft1(Vzzr,Vzzi,ne,n1,1);

	//计算归一化总梯度
	double Gxz=0;
	double *GH=new double[n1];
	for (int i = 0; i < n1; i++)
	{
		Gxz+=sqrt(Vxzr[i]*Vxzr[i]+Vzzr[i]*Vzzr[i]);
	}Gxz=Gxz/(n1);//均值
	for (int i = 0; i < n1; i++)
	{
		GH[i]=sqrt(Vxzr[i]*Vxzr[i]+Vzzr[i]*Vzzr[i])/Gxz;
	}
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=GH[j+nq];

	//销毁指针
	delete zr,zi,Vzzi,Vzzr,Vxzi,Vxzr,GH;
	return true;
}
bool NormalFullGradient_TSI(double* OriginalData,double* TransData,const int datanum,const double dx,const int N,const int L,const int K,const int hNum,const double dh,const double h1)
{
	return true;
}


/*=============================================
函数功能：频率域求一位位场数据向下延拓-吉洪诺夫正则化
参数：
	输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距
	rph：延拓高度
	TRP：正则化参数
=============================================*/
bool DownWardContinuation_Tik(double* OriginalData,double* TransData,const int datanum,const double dx,const double rph,double TRP)
{

	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	int nn=100;
	//====================================
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		w11=tp*(double)fabs(v1); w22=tp*(double)fabs((v2));
		AC11=AC22=(double)1;
		AC11*=(double)exp(-w11*h)/(exp(-2*w11*h)+TRP); AC22*=(double)exp(-w22*h)/(exp(-2*w22*h)+TRP);	
		//AC11*=(double)exp(w11*h)*(1-pow(TRP/exp(-2*w11*h)+TRP,nn)); AC22*=(double)exp(w22*h)*(1-pow(TRP/exp(-2*w22*h)+TRP,nn));	
		zr[j]=zr[j]*AC11; zi[j]=zi[j]*AC11;
		zr[j1]=zr[j1]*AC22; zi[j1]=zi[j1]*AC22;
	}
	//傅里叶逆变换
	fft1(zr,zi,ne,n1,1);
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=zr[j+nq];

	//销毁指针
	delete zr,zi;
	return true;
}

/*=============================================
函数功能：频率域求一位位场数据z方向一阶导数
参数：
	输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距
=============================================*/
bool Derivative_Z(double* OriginalData,double* TransData,const int datanum,const double dx)
{

	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,pjg,tp;
	double v1,v2;
	double x11,x22;
	int j1;

	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	for(int j=0;j<n1;j++)
	{ 
		zi[j]=(double)0; 
	}
	//=====扩边=================================================================================================

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	double temp;
	int order=1;
	//====================================
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=tp*fabs(v1);x22=tp*fabs(v2);
		zr[j]=zr[j]*pow(x11,order); zi[j]=zi[j]*pow(x11,order);
		zr[j1]=zr[j1]*pow(x22,order); zi[j1]=zi[j1]*pow(x22,order);
	}
	//傅里叶逆变换
	fft1(zr,zi,ne,n1,1);
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=zr[j+nq]*1E4;

	//销毁指针
	delete zr,zi;
	return true;
}

/*=============================================
函数功能：频率域求一位位场数据x方向一阶导数
参数：
	输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距
=============================================*/
bool Derivative_X(double* OriginalData,double* TransData,const int datanum,const double dx,int Order)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,pjg,tp;
	double v1,v2;
	double x11,x22;
	int j1;

	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	for(int j=0;j<n1;j++)
	{ 
		zi[j]=(double)0; 
	}
	//=====扩边=================================================================================================

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	double temp;
	//====================================
	/*switch (Order%4)
	{
	case 1:
		{
			for(int j=0;j<(n1/2);j++)
			{	
				v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
				x11=pow(tp*v1,Order);x22=pow(tp*v2,Order);
				temp=zr[j];zr[j]=zi[j]*x11; zi[j]=temp*x11;
				temp=zr[j1];zr[j1]=zi[j1]*x22; zi[j1]=temp*x22;
			}
		}
		break;
	case 2:
		{
			for(int j=0;j<(n1/2);j++)
			{	
				v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
				x11=pow(tp*v1,Order);x22=pow(tp*v2,Order);
				temp=zr[j];zr[j]=zi[j]*x11; zi[j]=temp*x11;
				temp=zr[j1];zr[j1]=zi[j1]*x22; zi[j1]=temp*x22;
			}
		}
		break;
	case 0:
		{
			for(int j=0;j<(n1/2);j++)
			{	
				v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
				x11=pow(tp*v1,Order);x22=pow(tp*v2,Order);
				zr[j]=zi[j]*x11; zi[j]*=x11;
				zr[j1]=zi[j1]*x22; zi[j1]=zi[j1]*x22;
			}
		}
		break;
	default:
		break;
	}*/
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=pow(tp*v1,Order);x22=pow(tp*v2,Order);
		temp=zr[j];zr[j]=-zi[j]*x11; zi[j]=temp*x11;
		temp=zr[j1];zr[j1]=-zi[j1]*x22; zi[j1]=temp*x22;
	}
	//傅里叶逆变换
	fft1(zr,zi,ne,n1,1);
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=zr[j+nq]*1E4;

	//销毁指针
	delete zr,zi;
	return true;
}

/*=============================================
函数功能：差分法向下延拓回返迭代反演场源深度
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距;
	P:回返次数
	rph：延拓高度
	参考文献：杨辉等（1999），重力异常视深度滤波及应用
=============================================*/
bool SourceDepthInversion(double* OriginalData,double* TransData,const int datanum,const double dx,const double P,const double rph)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//====================================
	double temp;
	double x11,x22;
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=tp*(double)v1; x22=tp*(double)v2;
		w11=fabs(x11); w22=fabs(x22);
		AC11=AC22=(double)1;
		AC11*=(double)pow(exp(-w11*h),P)*pow((4-exp(-w11*h)-2*cos(w11*h)),P+1); AC22*=(double)pow(exp(-w22*h),P)*pow((4-exp(-w22*h)-2*cos(w22*h)),P+1);

		zr[j]=zr[j]*AC11;zi[j]=zi[j]*AC11;
		zr[j1]=zr[j1]*AC22;zi[j1]=zi[j1]*AC22;
	}
	//傅里叶逆变换
	fft1(zr,zi,ne,n1,1);


	//..变换后数据赋值
	for(int j=0;j<n;j++)
	{
		TransData[j]=/*Vxzr[j]*1E4*/zr[j+nq];
	}
	//销毁指针
	delete zr,zi;
	return true;
}

/*=============================================
函数功能：差分法向下延拓
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组(列数与originaldata相同，行数为断面层数）
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum_max:谐波数最大值（最小值默认为1）
	h2/dh/h1：断面下层z坐标；相邻两层间距；顶层z坐标
=============================================*/
bool SourceDepthInversion(double* OriginalData,double* TransData,const int datanum,const double dx,const double P,const int hNum,const double dh,const double h1)
{
	double *GHdata=new double[datanum];
	double rph=0;
	for (int j = 0; j < hNum; j++)
	{
		rph=h1+j*dh;
		SourceDepthInversion(OriginalData,GHdata,datanum,dx,P,rph);
		for (int i = 0; i < datanum; i++)
		{
			TransData[j*datanum+i]=GHdata[i];
		}
	}
	delete GHdata;
	return true;
}


/*=============================================
函数功能：输出功率谱
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组(列数与originaldata相同，行数为断面层数）
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum_max:谐波数最大值（最小值默认为1）
	h2/dh/h1：断面下层z坐标；相邻两层间距；顶层z坐标
=============================================*/
bool AmplitudeSpectrum(double* OriginalData,const int datanum,const double dx,char* outfilename)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	//Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);
	//=====扩边=================================================================================================
	pjg=(double)0;
	pjg=(zr[0]+zr[n-1])/(double)2.;
	
	for(int j=(n-1);j>=0;j--)
	{ 
		zr[j+nq]=zr[j]-pjg;
	}

	for(int j=0;j<nq;j++)
		{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
			zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
		}zr[n1-1]=zr[0];
	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//fft1(zr,zi,ne,n1,1);
	//====================================
	double f=0,Amplitude;
	FILE* fp;
	if ((fp=fopen(outfilename,"w"))==NULL)
	{
		printf("打开%s失败\n");
		return false;
	}
	for(int j=0;j<n1/2;j++)
	{
		f=j/en+1E-10; 
		Amplitude=sqrt(zr[j]*zr[j]+zi[j]*zi[j]);
		fprintf(fp,"%lf\t%lf\n",f,log(Amplitude));
	}fclose(fp);

	//销毁指针
	delete zr,zi;
	return true;
}

bool Derivative_X_space(double* OriginalData,double* TransData,const int datanum,const double dx,int Order)
{
	double *temptransdata=new double[datanum],*ExtensionData=new double[datanum];
	int extennum=datanum+20;
	int halfextennum=(extennum-datanum)/2;
	Extension_Polynomial_1D(OriginalData,ExtensionData,dx,datanum,extennum,2);
	for (int k = 0; k < Order; k++)
	{
		for (int i = 0; i < extennum-1; i++)
		{
			temptransdata[i]=(ExtensionData[i+1]-ExtensionData[i])/dx;
		}
		for (int i = 0; i < extennum-1; i++)
		{
			ExtensionData[i]=temptransdata[i];
		}
	}
	
	for (int i = 0; i < datanum; i++)
	{
		TransData[i]=temptransdata[i];
	}
	return true;
}

bool Derivative_Z(double* OriginalData,double* TransData,const int datanum,const double dx,const int hNum,const double dh,const double h1)
{
	double *GHdata=new double[datanum],*contidata=new double[datanum];
	double rph=0;
	for (int j = 0; j < hNum; j++)
	{
		rph=h1+j*dh;
		UpWardContinuation(OriginalData,contidata,datanum,dx,rph);
		//DownWardContinuation_Tik(OriginalData,contidata,datanum,dx,rph,0.001);
		//Derivative_Z(contidata,GHdata,datanum,dx);
		Derivative_X(contidata,GHdata,datanum,dx,1);
		for (int i = 0; i < datanum; i++)
		{
			TransData[j*datanum+i]=GHdata[i];
		}
	}
	delete GHdata;
	return true;
}

/*===============================================================================
函数功能：频率域求一位位场数据归一化增强的解析信号振幅（主要是重力，磁场数据有待研究）
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum:谐波数
	rph：延拓高度
===============================================================================*/
bool NormalEAS(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	double* Vxzr=new double[n1],*Vxzi=new double[n1],*Vzzr=new double[n1],*Vzzi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//====================================
	double temp;
	double x11,x22;
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=tp*(double)v1; x22=tp*(double)v2;
		w11=fabs(x11); w22=fabs(x22);
		AC11=AC22=(double)1;
		AC11*=(double)exp(-w11*h)/(exp(-2*w11*h)+TRP); AC22*=(double)exp(-w22*h)/(exp(-2*w22*h)+TRP);	
		////计算x方向导数加延拓
		temp=zr[j];Vxzr[j]=-zi[j]*AC11*x11*pow(w11,2);Vxzi[j]=temp*AC11*x11*pow(w11,2);
		temp=zr[j1];Vxzr[j1]=-zi[j1]*AC22*x22*pow(w22,2);Vxzi[j1]=temp*AC22*x22*pow(w22,2);
		////计算z方向导数加延拓
		Vzzr[j]=zr[j]*AC11*w11;Vzzi[j]=zi[j]*AC11*w11;
		Vzzr[j1]=zr[j1]*AC22*w22;Vzzi[j1]=zi[j1]*AC22*w22;
	}
	//傅里叶逆变换
	fft1(Vxzr,Vxzi,ne,n1,1);
	fft1(Vzzr,Vzzi,ne,n1,1);

	//..缩边
	for(int j=0;j<n;j++)
	{
		Vxzr[j]=/*Vzzr[j+nq]*1E4*/Vxzr[j+nq];
		Vzzr[j]=/*Vzzr[j+nq]*1E4*/Vzzr[j+nq];
	}
	//计算归一化总梯度
	double Gxz=0;
	double *GH=new double[n1],*GH0=new double[n1];
	for (int i = 0; i < n1; i++)
	{
		GH0[i]=sqrt(Vxzr[i]*Vxzr[i]+Vzzr[i]*Vzzr[i]);
	}
	//===========窗口归一化===============================
	/*int WindowSize=datanum/2;
	int HalfWindowSize=WindowSize/2;
	double GHmax;
	Gxz=0;
	GHmax=GH0[0];
	for (int j = 0; j < WindowSize; j++)
	{
		if (GHmax<GH0[j])
		{
			GHmax=GH0[j];
		}	
	}
	for (int i = 0; i < HalfWindowSize; i++)
	{
		GH[i]=GH0[i]/GHmax;
	}

	for (int i = 0; i < datanum-WindowSize; i++)
	{
		GHmax=GH0[i];
		for (int j = 0; j < WindowSize; j++)
		{
			if (GHmax<GH0[i+j])
			{
				GHmax=GH0[i+j];
			}	
		}
		GH[i+HalfWindowSize]=GH0[i+HalfWindowSize]/GHmax;
	}

	for (int j = datanum-HalfWindowSize; j < datanum; j++)
	{
		GH[j]=GH0[j]/GHmax;
	}*/
	//===================================================
	//=========算数平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz+=GH0[i];
	//}Gxz=Gxz/(n);//均值
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================

	//=========调和平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz+=1.0/GH0[i];
	//}Gxz=n/Gxz;//
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================

	//=========几何平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz*=GH0[i];
	//}Gxz=pow(Gxz,1.0/n);//
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================
	switch (NormalizeType)
	{
	case 0:			//算术平均归一
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=GH0[i];
			}Gxz=Gxz/(n);//均值
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
		break;
	case 1:		//调和平均归一
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=1.0/GH0[i];
			}Gxz=n/Gxz;//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
		break;
	default:
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=1.0/GH0[i];
			}Gxz=n/Gxz;//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
		break;
	}
	//..变换后数据赋值
	for(int j=0;j<n;j++)
	{
		TransData[j]=/*Vxzr[j]*1E4*/GH0[j];
	}
	//销毁指针
	delete zr,zi,Vzzi,Vzzr,Vxzi,Vxzr,GH,GH0;
	return true;
}

/*=============================================
函数功能：频率域求一位位场数据归一EAS（主要是重力，磁场数据有待研究）
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组(列数与originaldata相同，行数为断面层数）
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum_max:谐波数最大值（最小值默认为1）
	h2/dh/h1：断面下层z坐标；相邻两层间距；顶层z坐标
=============================================*/
bool NormalEAS(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const int hNum,const double dh,const double h1,int NormalizeType)
{
	double *GHdata=new double[datanum];
	double rph=0;
	for (int j = 0; j < hNum; j++)
	{
		rph=h1+j*dh;
		NormalEAS(OriginalData,GHdata,datanum,dx,TRP,rph,NormalizeType);
		for (int i = 0; i < datanum; i++)
		{
			TransData[j*datanum+i]=GHdata[i];
		}
	}
	delete GHdata;
	return true;
}


/*=============================================
函数功能：利用Hilbert变换求一维位场数据z方向一阶导数
参数：
	输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距
=============================================*/
bool Derivative_Z_DHT(double* OriginalData,double* TransData,const int datanum,const double dx)
{

	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,pjg,tp;
	double v1,v2;
	double x11,x22;
	int j1;

	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	for(int j=0;j<n1;j++)
	{ 
		zi[j]=(double)0; 
	}
	//=====扩边=================================================================================================

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	double temp;
	int order=1;
	//====================================
	//求z方向导数的频谱
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=tp*fabs(v1);x22=tp*fabs(v2);
		zr[j]=zr[j]*pow(x11,order); zi[j]=zi[j]*pow(x11,order);
		zr[j1]=zr[j1]*pow(x22,order); zi[j1]=zi[j1]*pow(x22,order);
	}
	//利用频谱和DHT计算X方向导数
	double* HilbertTemp=new double[n1];
	for(int j=0;j<(n1);j++)
	{	
		HilbertTemp[j]=0;
		for (int i = 0; i < (n1/2); i++)
		{
			v1=i/en+1E-10; 
			x11=tp*fabs(v1);
			HilbertTemp[j]+=zi[i]*cos(x11*j*dx)+zr[i]*sin(x11*j*dx);
		}
	}
	//..变换后数据赋值
	for(int j=0;j<n;j++)
		TransData[j]=HilbertTemp[j+nq]*1E4;

	//销毁指针
	delete zr,zi,HilbertTemp;
	return true;

}


/*===============================================================================
函数功能：频率域求一位位场数据归一化局部波数
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距;
	HarmonicWaveNum:谐波数
	rph：延拓高度
===============================================================================*/
bool NormalLocalWaveNumber(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType)
{
	int n=datanum;
	int ne,n1,nq;
//To determine the size of arrays for FFT ///////////////////////////////////////
		
	ne=4; n1=(int)(pow(2,ne)+.001);
	while(n>=(n1-12))
	{   ne++; n1=(int)(pow(2,ne)+.001);
	}
	nq=(int)((n1-n)/2+.001);//半边扩变数
		
	int nj=n1;		//扩边后的总数
	double en,h,pjg,tp;
	double AC11,AC22;
	double v1,v2;
	double w11,w22;
	int j1;

	h=fabs(rph);
	en=dx*(double)(n1-1); 
	tp=(double)(2*PI);

	//---
	double* zr=new double[n1],*zi=new double[n1];
	//double* Vxzr=new double[n1],*Vxzi=new double[n1],*Vzzr=new double[n1],*Vzzi=new double[n1];
	double *Txzr=new double[n1],*Txzi=new double[n1],*Txr=new double[n1],*Txi=new double[n1],*Txxr=new double[n1],
		*Txxi=new double[n1],*Tzr=new double[n1],*Tzi=new double[n1],*Tzzr=new double[n1],*Tzzi=new double[n1];
	double* Kx=new double[n1],*Kz=new double[n1];
	for (int i = 0; i < datanum; i++)
	{
		zr[i]=OriginalData[i];
	}
	//=====扩边=================================================================================================
	//pjg=(double)0;
	//pjg=(zr[0]+zr[n-1])/(double)2.;
	//
	//for(int j=(n-1);j>=0;j--)
	//{ 
	//	zr[j+nq]=zr[j]-pjg;
	//}

	//for(int j=0;j<nq;j++)
	//	{   zr[j]=zr[nq]*(double)(1-cos(((double)PI*(double)j/(double)nq)))/(double)2;
	//		zr[j+n+nq]=zr[n+nq-1]*(double)(1+cos(((double)PI*(double)(j+1)/(double)(n1-n-nq))))/(double)2;
	//	}zr[n1-1]=zr[0];
	for(int j=0;j<n1;j++)
		{ 
			zi[j]=(double)0; 
		}

	Extension_Polynomial_1D(OriginalData,zr,dx,datanum,nj,2);

	//==========================================================================================================
	//傅里叶变换
	fft1(zr,zi,ne,n1,0);
	//====================================
	double temp;
	double a=1;int IterationNum=100;
	double x11,x22;
	for(int j=0;j<(n1/2);j++)
	{	
		v1=j/en+1E-10; v2=(j-n1/2)/en+1E-10; j1=n1/2+j;
		x11=tp*(double)v1; x22=tp*(double)v2;
		w11=fabs(x11); w22=fabs(x22);
		AC11=AC22=(double)1;
		AC11*=(double)exp(-w11*h)/(exp(-2*w11*h)+TRP); AC22*=(double)exp(-w22*h)/(exp(-2*w22*h)+TRP);	
		//AC11*=(double)exp(w11*h)*(1-pow(1-a*exp(-2*h*w11),IterationNum)); AC22*=(double)exp(w22*h)*(1-pow(1-a*exp(-2*h*w22),IterationNum));
		////计算x方向导数加延拓
		temp=zr[j];Txzr[j]=-zi[j]*AC11*x11*w11;Txzi[j]=temp*AC11*x11*w11;
		temp=zr[j1];Txzr[j1]=-zi[j1]*AC22*x22*w22;Txzi[j1]=temp*AC22*x22*w22;

		temp=zr[j];Txr[j]=-zi[j]*AC11*x11;Txi[j]=temp*AC11*x11;
		temp=zr[j1];Txr[j1]=-zi[j1]*AC22*x22;Txi[j1]=temp*AC22*x22;

		temp=zr[j];Txxr[j]=-zi[j]*AC11*x11*x11;Txxi[j]=temp*AC11*x11*x11;
		temp=zr[j1];Txxr[j1]=-zi[j1]*AC22*x22*x22;Txxi[j1]=temp*AC22*x22*x22;
		////计算z方向导数加延拓
		Tzr[j]=zr[j]*AC11*w11;Tzi[j]=zi[j]*AC11*w11;
		Tzr[j1]=zr[j1]*AC22*w22;Tzi[j1]=zi[j1]*AC22*w22;

		Tzzr[j]=zr[j]*AC11*w11*w11;Tzzi[j]=zi[j]*AC11*w11*w11;
		Tzzr[j1]=zr[j1]*AC22*w22*w22;Tzzi[j1]=zi[j1]*AC22*w22*w22;
	}
	//傅里叶逆变换
	fft1(Txzr,Txzi,ne,n1,1);
	fft1(Txr,Txi,ne,n1,1);
	fft1(Txxr,Txxi,ne,n1,1);
	fft1(Tzr,Tzi,ne,n1,1);
	fft1(Tzzr,Tzzi,ne,n1,1);

	//..缩边
	for(int j=0;j<n;j++)
	{
		Txzr[j]=/*Vzzr[j+nq]*1E4*/Txzr[j+nq];
		Txr[j]=/*Vzzr[j+nq]*1E4*/Txr[j+nq];
		Txxr[j]=/*Vzzr[j+nq]*1E4*/Txxr[j+nq];
		Tzr[j]=/*Vzzr[j+nq]*1E4*/Tzr[j+nq];
		Tzzr[j]=/*Vzzr[j+nq]*1E4*/Tzzr[j+nq];
	}
	//计算局部波数
	double Gxz=0;
	double *GH=new double[n1],*GH0=new double[n1];
	for (int i = 0; i < n1; i++)
	{
		Kx[i]=(Txzr[i]*Txr[i]-Txxr[i]*Tzr[i])/sqrt(Txr[i]*Txr[i]+Tzr[i]*Tzr[i]);
		Kz[i]=-(Txzr[i]*Tzr[i]-Tzzr[i]*Txr[i])/sqrt(Txr[i]*Txr[i]+Tzr[i]*Tzr[i]);
		GH0[i]=sqrt(Kx[i]*Kx[i]+Kz[i]*Kz[i]);
	}
	//===========窗口归一化===============================
	/*int WindowSize=datanum/2;
	int HalfWindowSize=WindowSize/2;
	double GHmax;
	Gxz=0;
	GHmax=GH0[0];
	for (int j = 0; j < WindowSize; j++)
	{
		if (GHmax<GH0[j])
		{
			GHmax=GH0[j];
		}	
	}
	for (int i = 0; i < HalfWindowSize; i++)
	{
		GH[i]=GH0[i]/GHmax;
	}

	for (int i = 0; i < datanum-WindowSize; i++)
	{
		GHmax=GH0[i];
		for (int j = 0; j < WindowSize; j++)
		{
			if (GHmax<GH0[i+j])
			{
				GHmax=GH0[i+j];
			}	
		}
		GH[i+HalfWindowSize]=GH0[i+HalfWindowSize]/GHmax;
	}

	for (int j = datanum-HalfWindowSize; j < datanum; j++)
	{
		GH[j]=GH0[j]/GHmax;
	}*/
	//===================================================
	//=========算数平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz+=GH0[i];
	//}Gxz=Gxz/(n);//均值
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================

	//=========调和平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz+=1.0/GH0[i];
	//}Gxz=n/Gxz;//
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================

	//=========几何平均归一化==================================
	//for (int i = 0; i < n; i++)
	//{
	//	Gxz*=GH0[i];
	//}Gxz=pow(Gxz,1.0/n);//
	//for (int i = 0; i < n; i++)
	//{
	//	GH[i]=GH0[i]/Gxz;
	//	//GH[i]=atan2(Vxzr[i],Vzzr[i]);
	//}
	//====================================================
	switch (NormalizeType)
	{
	case 0:			//算术平均归一
		{
			//int WindowSize=n/2;//总长度的十分之一
			//int HalfWindowSize=WindowSize/2;
			//WindowSize=HalfWindowSize*2;
			//for (int i = 0; i < n; i++)
			//{
			//	Gxz+=GH0[i];
			//}Gxz=Gxz/(n);//均值
			////左边归一
			//for (int i = 0; i < HalfWindowSize; i++)
			//{
			//	GH[i]=0;
			//}
			////右边归一
			//for (int i = n-HalfWindowSize; i < n; i++)
			//{
			//	GH[i]=0;
			//}
			////中间滑动窗口归一
			//for (int i = HalfWindowSize; i < n-HalfWindowSize; i++)
			//{
			//	Gxz=0;
			//	for (int k = i-HalfWindowSize; k < i+HalfWindowSize; k++)
			//	{
			//		Gxz+=GH0[k];
			//	}Gxz=Gxz/(WindowSize);
			//	GH[i]=GH0[i]/Gxz;
			//}

			for (int i = 0; i < n; i++)
			{
				Gxz+=GH0[i];
			}Gxz=Gxz/(n);//均值
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}

		}
		break;
	case 1:		//调和平均归一
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=1.0/GH0[i];
			}Gxz=n/Gxz;//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
		break;
	default:
		{
			for (int i = 0; i < n; i++)
			{
				Gxz+=1.0/GH0[i];
			}Gxz=n/Gxz;//
			for (int i = 0; i < n; i++)
			{
				GH[i]=GH0[i]/Gxz;
				//GH[i]=atan2(Vxzr[i],Vzzr[i]);
			}
		}
		break;
	}
	//..变换后数据赋值
	for(int j=0;j<n;j++)
	{
		TransData[j]=/*Vxzr[j]*1E4*/Kz[j];
	}
	//销毁指针
	delete zr,zi,Txr,Txi,Txxr,Txxi,Txzr,Txzi,Tzr,Tzi,Tzzr,Tzzi,GH,GH0;
	return true;
}



/*===============================================================================
函数功能：频率域求一位位场数据熵滤波分场
参数：
	OriginalData/TransData：输入位场数据数组；转换后的数组
	datanum：数据个数
	dx：点距;
	Winsize：窗口大小
===============================================================================*/
bool Separation_EntropyFilter(double* OriginalData,double* TransData,const int datanum,const double dx,int WinSize)
{
	int np=datanum+WinSize;
	double* ExtensionData=new double[datanum+np];
	Extension_Polynomial_1D(OriginalData,ExtensionData,dx,datanum,np,1);		//扩边
	
	double mindata,maxdata,meandata,rmax;
	double *Entropy=new double[WinSize],*EntropyCoe=new double[WinSize];
	double SumEntropy=0,RegData=0;
	for (int i = 0; i < datanum; i++)
	{
		double SumTemp=0;
		mindata=ExtensionData[i],maxdata=ExtensionData[i];
		for (int j = i; j < i+WinSize; j++)		//计算均值，最大值，最小值
		{
			SumTemp+=ExtensionData[j];
			if (ExtensionData[j]>maxdata)
			{
				maxdata=ExtensionData[j];
			}
			if (ExtensionData[j]<mindata)
			{
				mindata=ExtensionData[j];
			}
		}meandata=SumTemp/WinSize;
		rmax=fabs(ExtensionData[i]-meandata);
		for (int j = i; j < i+WinSize; j++)		//计算rmax
		{
			if (fabs(ExtensionData[j]-meandata)>maxdata)
			{
				rmax=fabs(ExtensionData[j]-meandata);
			}
		}
		SumEntropy=0;
		for (int j = i; j < i+WinSize; j++)		//计算熵值
		{
			Entropy[j-i]=fabs(ExtensionData[j]-meandata)/rmax;
			printf("%lf\n",Entropy[j-i]);
			SumEntropy+=1-Entropy[j-i];
		}
		RegData=0;
		for (int j = i; j < i+WinSize; j++)		//计算熵系数
		{
			EntropyCoe[j-i]=(1-Entropy[j-i])/SumEntropy;
			RegData+=ExtensionData[j]*EntropyCoe[j-i];
		}
		
		TransData[i]=OriginalData[i]/meandata;
		printf("原始值:%lf\n平均值:%lf\n熵滤波:%lf\n熵系数:%lf\n",OriginalData[i],meandata,RegData,SumEntropy);
		return true;
	}

	delete ExtensionData,Entropy;
	return true;
}
