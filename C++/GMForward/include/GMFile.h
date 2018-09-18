#pragma once

#ifndef GMFILE
#define GMFILE


class CGMFile
{
public:
	CGMFile(void);
	~CGMFile(void);
private:
	CString m_FileName;
public:
	int SetFileName(CString filename);
	int CGMFile::ReadGrd(CGM_GRDDATAINFO& datainfo,vector<double>& GMData);
	int ReadDat(vector<double>& dataX, vector<double>& dataY);
	int SaveGrd(vector<double> data, CGM_GRDDATAINFO datainfo);
	int SaveDat(vector<double> dataX,vector<double> dataY);
	int SaveDatFromGrd(CGM_GRDDATAINFO datainfo,vector<double> data, int XorY , int index);
	int SaveDatFromGrd(CGM_GRDDATAINFO datainfo,double** data, int XorY , int index);
	int SaveGrd(double* data, CGM_GRDDATAINFO datainfo);
	int SaveGrd(double** data, CGM_GRDDATAINFO datainfo);
	int ReadGrd(CGM_GRDDATA& GMData);
	int dsa(void);
};


#endif
