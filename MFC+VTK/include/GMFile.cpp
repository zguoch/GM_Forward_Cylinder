#include "GM3DForward.h"
#include "GMFile.h"
#include "comutil.h "
CGMFile::CGMFile(void)
	: m_FileName(_T(""))
{
}


CGMFile::~CGMFile(void)
{
}


int CGMFile::SetFileName(CString filename)
{
	m_FileName=filename;
	return 0;
}


int CGMFile::ReadGrd(CGM_GRDDATAINFO& datainfo,vector<double>& GMData)
{
	FILE* FP;
	if ((FP = fopen((_bstr_t)m_FileName, "r")) == NULL)
	{
		//AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
		return 0;
	}
	//1.1读取文件头信息
	char dsaa[6];
	fscanf(FP,"%s",dsaa);
	CString dsaaStr=(PCTSTR)dsaa;
	if (dsaaStr.Left(4)!=_T("DSAA"))
	{
		//AfxMessageBox(_T("无法读取\n请检查Grd文件头，此文件有可能是二进制文件"));
		fclose(FP);
		return 0;
	}else
	{
		double tempdata;
		int number_x,number_y;
		fscanf(FP,"%d%d%",&number_x,&number_y);
		datainfo.SetNumber(number_x,number_y);
		for (int i = 0; i < 4; i++)
		{
			fscanf(FP,"%lf",&datainfo.m_AxisBounds[i]);
		}
		fscanf(FP,"%lf%lf",&datainfo.m_AxisBounds[4],&datainfo.m_AxisBounds[5]);
		//1.2读取重力数据
		GMData.clear();
		int datanum=datainfo.GetNumber_x()*datainfo.GetNumber_y();
		for (int i=0;i<datanum;i++)
		{
			fscanf(FP,"%lf",&tempdata);
			GMData.push_back(tempdata);
		}fclose(FP);//关闭文件
	}
	return 1;
}
int CGMFile::ReadGrd(CGM_GRDDATA& GMData)
{
	/*首先判断GMData数组是否为空*/

	FILE* FP;
	if ((FP = fopen((_bstr_t)m_FileName, "r")) == NULL)
	{
		MessageBox(NULL,_T("打开")+m_FileName+_T("失败!"),_T("打开文件提示"),MB_OK);
		return 0;
	}
	//1.1读取文件头信息
	char dsaa[6];
	double AxisBounds[6];
	fscanf(FP, "%s", dsaa);
	CString dsaaStr = (PCTSTR)(_bstr_t)dsaa;
	if (dsaaStr.Left(4) != _T("DSAA"))
	{
		MessageBox(NULL, _T("无法读取\n请检查Grd文件头，此文件有可能是二进制文件"), _T("打开文件提示"), MB_OK);
		fclose(FP);
		return 0;
	}
	else
	{
		double tempdata;
		int number_x, number_y;
		fscanf(FP, "%d%d%", &number_x, &number_y);
		for (int i = 0; i < 6; i++)
		{
			fscanf(FP, "%lf", &AxisBounds[i]);
		}
		//1.2读取重力数据
		GMData.Initialize(number_x, number_y, AxisBounds);
		for (int i = 0; i<number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				fscanf(FP, "%lf", &tempdata);
				GMData.SetData(i,j,tempdata);
			}
		}fclose(FP);//关闭文件
	}
	return 1;
}

//int CGMFile::ReadGrd(double* GMData, CGM_GRDDATAINFO& datainfo)
//{
//	FILE* FP;
//	if ((FP=fopen(m_FileName,"r"))==NULL)
//	{
//		//AfxMessageBox(("打开")+m_FileName+("失败!"));
//		return 0;
//	}
//	//1.1读取文件头信息
//	char dsaa[6];
//	fscanf(FP,"%s",dsaa);
//	CString dsaaStr=(PCTSTR)dsaa;
//	if (dsaaStr.Left(4)!=_T("DSAA"))
//	{
//		//AfxMessageBox(_T("无法读取\n请检查Grd文件头，此文件有可能是二进制文件"));
//		fclose(FP);
//		return 0;
//	}else
//	{
//		double tempdata;
//		fscanf(FP,"%d%d%lf%lf%lf%lf%lf%lf",&datainfo.columns,&datainfo.rows,&datainfo.gridXmin,&datainfo.gridXmax,&datainfo.gridYmin,&datainfo.gridYmax,
//		&datainfo.gridZmin,&datainfo.gridZmax);
//		//1.2读取重力数据
//		int datanum=datainfo.rows*datainfo.columns;
//		if(GMData)delete GMData;
//		GMData=new double[datanum];
//		for (int i=0;i<datanum;i++)
//		{
//			fscanf(FP,"%lf",&tempdata);
//			GMData[i]=tempdata;
//		}fclose(FP);//关闭文件
//	}
//	datainfo.deltX=(datainfo.gridXmax-datainfo.gridXmin)/(datainfo.columns-1);
//	datainfo.deltY=(datainfo.gridYmax-datainfo.gridYmin)/(datainfo.rows-1);
//	return 1;
//}
//
//int CGMFile::ReadDat(vector<double>& dataX, vector<double>& dataY)
//{
//	FILE* FP;
//	double tempx,tempy;
//	if ((FP=fopen((_bstr_t)m_FileName,"r"))==NULL)
//	{
//		AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
//		return 0;
//	}	
//	dataX.clear();dataY.clear();
//	while (!feof(FP))
//	{
//		fscanf(FP,"%lf%lf",&tempx,&tempy);
//		dataX.push_back(tempx);
//		dataY.push_back(tempy);
//	}
//	if (dataX[dataX.size()-1]==dataX[dataX.size()-2] && dataY[dataY.size()-1]==dataY[dataY.size()-2])
//	{
//		dataX.pop_back();
//		dataY.pop_back();
//	}
//	fclose(FP);
//	return 1;
//}
//
//
//int CGMFile::SaveGrd(vector<double> data, GridDataInfo datainfo)
//{
//	FILE* FP;
//	if ((FP=fopen((_bstr_t)m_FileName,"w"))==NULL)
//	{
//		AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
//		return 0;
//	}	
//	fprintf(FP,"DSAA\n%d\t%d\n%lf\t%lf\n%lf\t%lf\n%lf\t%lf\n",datainfo.columns,datainfo.rows,datainfo.gridXmin,datainfo.gridXmax,datainfo.gridYmin,datainfo.gridYmax,datainfo.gridZmin,datainfo.gridZmax);
//	for (int i = 0; i < datainfo.rows*datainfo.columns; i++)
//	{
//		fprintf(FP,"%lf\t",data[i]);
//	}fclose(FP);
//	return 1;
//}
//
//int CGMFile::SaveGrd(double* data, GridDataInfo datainfo)
//{
//	FILE* FP;
//	if ((FP=fopen((_bstr_t)m_FileName,"w"))==NULL)
//	{
//		AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
//		return 0;
//	}	
//	fprintf(FP,"DSAA\n%d\t%d\n%lf\t%lf\n%lf\t%lf\n%lf\t%lf\n",datainfo.columns,datainfo.rows,datainfo.gridXmin,datainfo.gridXmax,datainfo.gridYmin,datainfo.gridYmax,datainfo.gridZmin,datainfo.gridZmax);
//	for (int i = 0; i < datainfo.rows*datainfo.columns; i++)
//	{
//		fprintf(FP,"%lf\t",data[i]);
//	}fclose(FP);
//	return 1;
//}
int CGMFile::SaveGrd(double** data, CGM_GRDDATAINFO datainfo)
{
	FILE* FP;
	if ((FP = fopen((_bstr_t)m_FileName, "w")) == NULL)
	{
		MessageBox(NULL,_T("打开")+m_FileName+_T("失败!"),_T("保存Grd警告"),MB_OK);
		return 0;
	}
	
	fprintf(FP,"DSAA\n%d\t%d\n%lf\t%lf\n%lf\t%lf\n%.16f\t%.16f\n",datainfo.GetNumber_x(),datainfo.GetNumber_y(),datainfo.m_AxisBounds[0],datainfo.m_AxisBounds[1],
		datainfo.m_AxisBounds[2],datainfo.m_AxisBounds[3],datainfo.m_Ranges[0],datainfo.m_Ranges[1]);
	for (int i = 0; i < datainfo.GetNumber_y(); i++)
	{
		for (int j = 0; j < datainfo.GetNumber_x(); j++)
		{
			fprintf(FP,"%.16f\t",data[i][j]);
		}
		fprintf(FP,"\n");
	}fclose(FP);
	return 1;
}

//int CGMFile::SaveDat(vector<double> dataX,vector<double> dataY)
//{
//	//直接保存CGMFIPSView中的数据X和数据Y向量
//	FILE* FP;
//	if ((FP=fopen((_bstr_t)m_FileName,"w"))==NULL)
//	{
//		AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
//		return 0;
//	}
//	for (int i = 0; i<dataX.size(); i++)
//	{
//		fprintf(FP,"%lf\t%lf\n",dataX[i],dataY[i]);
//	}
//	fclose(FP);
//	return 1;
//}
//
//
//int CGMFile::SaveDatFromGrd(GridDataInfo datainfo,vector<double> data, int XorY , int index)
//{
//	FILE* FP;
//	if ((FP=fopen((_bstr_t)m_FileName,"w"))==NULL)
//	{
//		AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
//		return 0;
//	}
//	switch (XorY)
//	{
//	case 0:
//		{
//			double dx=(datainfo.gridXmax-datainfo.gridXmin)/(double)(datainfo.columns-1);
//			for (int i = 0; i < datainfo.columns; i++)
//			{
//				fprintf(FP,"%lf\t%lf\n",datainfo.gridXmin+i*dx,data[index*datainfo.columns+i]);
//			}
//		}
//		break;
//	case 1:
//		{
//			double dy=(datainfo.gridYmax-datainfo.gridYmin)/(double)(datainfo.rows-1);
//			for (int i = 0; i < datainfo.rows; i++)
//			{
//				fprintf(FP,"%lf\t%lf\n",datainfo.gridYmin+i*dy,data[index+i*datainfo.columns]);
//			}
//		}
//		break;
//	default:
//		break;
//	}fclose(FP);
//	return 1;
//}
//
//int CGMFile::SaveDatFromGrd(GridDataInfo datainfo,double** data, int XorY , int index)
//{
//	FILE* FP;
//	if ((FP=fopen((_bstr_t)m_FileName,"w"))==NULL)
//	{
//		AfxMessageBox(_T("打开")+m_FileName+_T("失败!"));
//		return 0;
//	}
//	switch (XorY)
//	{
//	case 0:
//		{
//			double dx=(datainfo.gridXmax-datainfo.gridXmin)/(double)(datainfo.columns-1);
//			for (int i = 0; i < datainfo.columns; i++)
//			{
//				fprintf(FP,"%lf\t%lf\n",datainfo.gridXmin+i*dx,data[index][i]);
//			}
//		}
//		break;
//	case 1:
//		{
//			double dy=(datainfo.gridYmax-datainfo.gridYmin)/(double)(datainfo.rows-1);
//			for (int i = 0; i < datainfo.rows; i++)
//			{
//				fprintf(FP,"%lf\t%lf\n",datainfo.gridYmin+i*dy,data[i][index]);
//			}
//		}
//		break;
//	default:
//		break;
//	}fclose(FP);
//	return 1;
//}

