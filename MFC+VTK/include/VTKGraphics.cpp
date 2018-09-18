//#include "stdafx.h"
#include "VTKGraphics.h"
#include "vector"

using namespace std;
CVTKGraphics::CVTKGraphics(void)
	:m_FeatureAngle(60)
	,m_ContourNumber(40)
{
	RGBColorStruct rgb; 
	m_ColorArray.push_back({ RED });
	m_ColorArray.push_back({ GREEN });
	m_ColorArray.push_back({ BLUE });
	m_ColorArray.push_back({ YELLOW });
	m_ColorArray.push_back({ MAGENTA });
	m_ColorArray.push_back({ CYAN });
	m_ColorArray.push_back({ WHITE });
	m_ColorArray.push_back({ BLACK });
}


CVTKGraphics::~CVTKGraphics(void)
{
}

int CVTKGraphics::GrdToStructuredGridData(double** grddata, vtkStructuredGrid* sgrd,CGM_GRDDATAINFO grdinfo)
{
	int number_x=grdinfo.GetNumber_x(),number_y=grdinfo.GetNumber_y();
	int dims[3];dims[0]=number_x;dims[1]=number_y;dims[2]=1;
	double x[3];
	vtkPoints* points = vtkPoints::New();
	points->Allocate(number_x*number_y);

	double dx=fabs(grdinfo.m_AxisBounds[1]-grdinfo.m_AxisBounds[0])/(number_x-1);
	double dy=fabs(grdinfo.m_AxisBounds[2]-grdinfo.m_AxisBounds[3])/(number_y-1);
	//double dz=fabs(meshinfo.bounds[4]-meshinfo.bounds[5])/(meshinfo.NodeNum_z-1);
	int offset,jOffset;
	//for (int k = 0; k < meshinfo.NodeNum_z; k++)
	{
		x[2]=grdinfo.m_Height_data;
		for (int j = 0; j < number_y; j++)
		{
			x[1]=grdinfo.m_AxisBounds[2]+j*dy;
			jOffset = j * dims[0];
			for (int i = 0; i <number_x; i++)
			{
				x[0]=grdinfo.m_AxisBounds[0]+i*dx;
				offset = i + jOffset;
				points->InsertPoint(offset,x);
			}
		}
	}
	sgrd->SetDimensions(number_x,number_y,1);
	sgrd->SetPoints(points);
	points->Delete();

	//添加属性值
	vtkSmartPointer<vtkDoubleArray> scalardata=vtkSmartPointer<vtkDoubleArray>::New();
	int *dimension=sgrd->GetDimensions();
	//设置属性值
	scalardata->SetNumberOfComponents(1);
	scalardata->SetNumberOfTuples(dimension[0]*dimension[1]*dimension[2]);
	int nnn=0;
	
	for (int i = 0; i < dimension[1]; i++)
	{
		for (int j = 0; j < dimension[0]; j++)
		{
	
			scalardata->InsertTuple1(nnn,grddata[i][j]);
			nnn++;
		}
	}
	sgrd->GetPointData()->SetScalars(scalardata);

	//测试
	//3Dview
	//m_VTKGraphics.Surfer3D(sgrd,m_GM3DActor,m_GMColorTable);
	////contour
	//m_VTKGraphics.Contour(sgrd,m_GMFillContourActor,m_GMContourLineActor,m_GMContourLabelsActor,m_GMColorTable);

//输出到文件
	/*vtkSmartPointer<vtkStructuredGridWriter> sgrdwrite=vtkSmartPointer<vtkStructuredGridWriter>::New();
	sgrdwrite->SetFileName("测试目录\\结构化数据.dat");
	sgrdwrite->SetInputData(strgrd);
	sgrdwrite->Write();*/
	return 0;
}
int CVTKGraphics::GrdToStructuredGridData(CGM_GRDDATA gmdata, vtkStructuredGrid* sgrd)
{
	int number_x , number_y;
	gmdata.GetDimension(number_x, number_y);
	int dims[3]; dims[0] = number_x; dims[1] = number_y; dims[2] = 1;
	double x[3];
	vtkSmartPointer< vtkPoints> points = vtkSmartPointer< vtkPoints>::New();
	points->Allocate(number_x*number_y);
	double* axisbounds = gmdata.GetBounds();
	double dx = fabs(axisbounds[1] - axisbounds[0]) / (number_x - 1);
	double dy = fabs(axisbounds[2] - axisbounds[3]) / (number_y - 1);
	//double dz=fabs(meshinfo.bounds[4]-meshinfo.bounds[5])/(meshinfo.NodeNum_z-1);
	int offset, jOffset;
	//for (int k = 0; k < meshinfo.NodeNum_z; k++)
	{
		x[2] = gmdata.GetHeight();
		for (int j = 0; j < number_y; j++)
		{
			x[1] = axisbounds[2] + j*dy;
			jOffset = j * dims[0];
			for (int i = 0; i <number_x; i++)
			{
				x[0] = axisbounds[0] + i*dx;
				offset = i + jOffset;
				points->InsertPoint(offset, x);
			}
		}
	}
	sgrd->SetDimensions(number_x, number_y, 1);
	sgrd->SetPoints(points);

	//添加属性值
	vtkSmartPointer<vtkDoubleArray> scalardata = vtkSmartPointer<vtkDoubleArray>::New();
	int *dimension = sgrd->GetDimensions();
	//设置属性值
	scalardata->SetNumberOfComponents(1);
	scalardata->SetNumberOfTuples(dimension[0] * dimension[1] * dimension[2]);
	int nnn = 0;
	double* data = gmdata.GetData();
	for (int i = 0; i < dimension[1]; i++)
	{
		for (int j = 0; j < dimension[0]; j++)
		{

			scalardata->InsertTuple1(nnn, data[i*number_x+j]);
			nnn++;
		}
	}
	sgrd->GetPointData()->SetScalars(scalardata);

	//测试
	//3Dview
	//m_VTKGraphics.Surfer3D(sgrd,m_GM3DActor,m_GMColorTable);
	////contour
	//m_VTKGraphics.Contour(sgrd,m_GMFillContourActor,m_GMContourLineActor,m_GMContourLabelsActor,m_GMColorTable);

	//输出到文件
	/*vtkSmartPointer<vtkStructuredGridWriter> sgrdwrite=vtkSmartPointer<vtkStructuredGridWriter>::New();
	sgrdwrite->SetFileName("测试目录\\结构化数据.dat");
	sgrdwrite->SetInputData(strgrd);
	sgrdwrite->Write();*/
	return 0;
}
int CVTKGraphics::SaveGrd(vtkStructuredGrid* sgrd,CString filename)
{
	FILE* FP;
	if ((FP = fopen((_bstr_t)filename, "w")) == NULL)
	{
		MessageBox(NULL,_T("打开")+filename+_T("失败!"),_T(""),MB_OK);
		return 0;
	}
	double* bounds=sgrd->GetBounds();
	double* range=sgrd->GetScalarRange();
	int* dims=sgrd->GetDimensions();
	fprintf(FP,"DSAA\n%d\t%d\n%lf\t%lf\n%lf\t%lf\n%.16f\t%.16f\n",dims[0],dims[1],bounds[0],bounds[1],bounds[2],bounds[3],range[0],range[1]);
	vtkDoubleArray* scalardata=(vtkDoubleArray*)sgrd->GetPointData()->GetScalars();
	for (int i = 0; i < dims[1]; i++)
	{
		for (int j = 0; j < dims[0]; j++)
		{
			fprintf(FP,"%.16f\t",scalardata->GetValue(j+i*dims[0]));
		}
		fprintf(FP,"\n");
	}
	fclose(FP);
	return 1;
}

int CVTKGraphics::CreateColorTable(vtkLookupTable* colorTable, char* colorFile)
{
	//读取文件创建颜色表 
	int colornumber=255;
	double percentage;
	int r,g,b,a;
	vector<double>pervector;
	vector<int>rvector,gvector,bvector,aerfavector;
	char tempchr[20];
	int version1,version2;
	int count=0;

	FILE* colorfp;
	if ((colorfp=fopen(colorFile,"r"))==NULL)
	{
		MessageBox(NULL,_T("颜色表文件读取失败"),_T("出错提示"),MB_OK);
		colorTable->SetNumberOfColors(256);
		colorTable->SetHueRange(0.2,0.78);
		colorTable->Build();
	}else
	{
		fscanf(colorfp,"%s%d%d",tempchr,&version1,&version2);
		if (version1==1 && version2==1)
		{
			while(!feof(colorfp))
			{
				fscanf(colorfp,"%lf%d%d%d",&percentage,&r,&g,&b);
				pervector.push_back(percentage);
				rvector.push_back(r);
				gvector.push_back(g);
				bvector.push_back(b);
			}
		}else if(version1==2 && version2==1)
		{
			while(!feof(colorfp))
			{
				fscanf(colorfp,"%lf%d%d%d%d",&percentage,&r,&g,&b,&a);
				pervector.push_back(percentage);
				rvector.push_back(r);
				gvector.push_back(g);
				bvector.push_back(b);
				aerfavector.push_back(a);
			}
		}else
		{
			MessageBox(NULL,_T("颜色文件不是标准的Surfer颜色文件"),_T("出错提示"),MB_OK);
			colorTable->SetNumberOfColors(256);
			colorTable->SetHueRange(0.2,0.78);
			colorTable->Build();
			return 0;
		}
		if(pervector[pervector.size()-1]==pervector[pervector.size()-2])pervector.pop_back();//去掉文件末尾的回车
		colorTable->SetNumberOfColors(colornumber+1);

		//FILE* fp=fopen("测试colorscale.clr","w");fprintf(fp,"ColorMap 1 1\n");
		for (int i = 0; i < pervector.size()-1; i++)
		{
			int n1=(int)(colornumber*pervector[i]/100.0);
			int n2=(int)(colornumber*pervector[i+1]/100.0);
			int dn=n2-n1;
			double dr=(rvector[i+1]-rvector[i])/dn;
			double dg=(gvector[i+1]-gvector[i])/dn;
			double db=(bvector[i+1]-bvector[i])/dn;
			double dp=(pervector[i+1]-pervector[i])/dn;
			for (int j = n1; j <n2; j++)
			{
				colorTable->SetTableValue(j,(rvector[i]+(j-n1)*dr)/255.0,(gvector[i]+(j-n1)*dg)/255.0,(bvector[i]+(j-n1)*db)/255.0);
				//fprintf(fp,"%lf\t%d\t%d\t%d\n",pervector[i]+(j-n1)*dp,(int)(rvector[i]+(j-n1)*dr),(int)(gvector[i]+(j-n1)*dg),(int)(bvector[i]+(j-n1)*db));
			}
		}
		colorTable->SetTableValue(colornumber,(rvector[rvector.size()-1])/255.0,(gvector[gvector.size()-1])/255.0,(bvector[bvector.size()-1])/255.0);
		//fprintf(fp,"%lf\t%d\t%d\t%d\n",pervector[pervector.size()-1],(int)(rvector[rvector.size()-1]),(int)(gvector[gvector.size()-1]),(int)(bvector[bvector.size()-1]));
		//fclose(fp);
		colorTable->Build();
		fclose(colorfp);
	}

	return 0;
}
int CVTKGraphics::SetColorScale(vtkScalarBarActor* scalbar,char* title,double *pos,double width,double height,bool isHorizontal)
{
	vtkSmartPointer<vtkTextProperty>scalBarText=vtkSmartPointer<vtkTextProperty>::New();
	scalBarText->SetFontSize(10);
	scalBarText->SetColor(1,0,0);
	scalBarText->SetFontFamilyToTimes();
	scalbar->SetTitle(title);
	if (isHorizontal)
	{
		scalbar->SetOrientationToHorizontal();
	}else
	{
		scalbar->SetOrientationToVertical();
	}
	
	scalbar->SetTitleTextProperty(scalBarText);
	scalbar->SetLabelTextProperty(scalBarText);
	scalbar->SetHeight(height);
	scalbar->SetWidth(width);
	scalbar->SetPosition(pos);

	return 0;
}

int CVTKGraphics::SetCubeAxesProp(vtkCubeAxesActor* cubeaxes,vtkActor* outlineActor,vtkActor* BourderFaceActor,double* bounds)
{
	COLORREF Xcolor=RGB(255,0,0),Ycolor=RGB(0,255,0),Zcolor=RGB(0,0,255),m_OutlineColor=RGB(255-179,255-179,255-179);
	int fontsize=30;
	COLORREF BourderColor=RGB(100,0,100);
	double Bounds[6];
	for (int i = 0; i < 6; i++)
	{
		Bounds[i]=bounds[i];
	}
	
	//1 轮廓线
	 int nV = 2;      // No. of vertices
	int i;
	// Create points and cells for the spiral
	//1.1三维坐标赋值
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for(i = 0; i < nV; i++)
	{
		points->InsertPoint(i, Bounds[i],Bounds[i+2],Bounds[i+4]);
	}
	//1.2单元
	vtkSmartPointer<vtkCellArray> lines =
		vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(nV);
	for (i = 0; i < nV; i++)
	{
		lines->InsertCellPoint(i);
	}
	//1.3根据点和单元构建polydata
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);
	polyData->SetLines(lines);
	//轮廓
	vtkSmartPointer<vtkPolyDataNormals>normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputData(polyData);
	vtkSmartPointer<vtkOutlineFilter>outline = vtkSmartPointer<vtkOutlineFilter>::New();
	outline->SetInputConnection(normals->GetOutputPort());
	vtkSmartPointer<vtkPolyDataMapper>mapOutline = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapOutline->SetInputConnection(outline->GetOutputPort());
	outlineActor->SetMapper(mapOutline);
	outlineActor->GetProperty()->SetColor(GetRValue(m_OutlineColor)/255.0,GetGValue(m_OutlineColor)/255.0,GetBValue(m_OutlineColor)/255.0);
	outlineActor->GetProperty()->SetLineWidth(2);
	
	//2. 设置坐标轴
	cubeaxes->SetFlyModeToClosestTriad();
	cubeaxes->XAxisMinorTickVisibilityOff();		//关闭次级刻度
	cubeaxes->YAxisMinorTickVisibilityOff();
	cubeaxes->ZAxisMinorTickVisibilityOff();
	//if (Is3DImage)
	{
		cubeaxes->SetBounds(Bounds[0],Bounds[1],Bounds[2],Bounds[3],Bounds[4],Bounds[5]);
	}//else
	{
		//m_CubeAxes->SetBounds(Bounds[0],Bounds[1],Bounds[2],Bounds[3],0,0);
	}
	
	/*cubeaxes->SetXAxisRange(0,100);
	cubeaxes->SetYAxisRange(0,100);
	cubeaxes->SetZAxisRange(0,50);*/
	//m_CubeAxes->SetXLabelFormat("%6.4f");
	cubeaxes->SetXTitle("X");
	cubeaxes->SetYTitle("Y");
	cubeaxes->SetZTitle("Depth");
	//设置字体，颜色
	//if (Is3DImage)
	{
		for (int i = 0; i < 3; i++)
		{
			cubeaxes->GetTitleTextProperty(i)->SetFontSize(15);
			cubeaxes->GetTitleTextProperty(i)->SetFontFamilyToTimes();
			cubeaxes->GetLabelTextProperty(i)->SetFontFamilyToTimes();
			cubeaxes->GetLabelTextProperty(i)->SetFontSize(fontsize);
		}
		cubeaxes->GetXAxesLinesProperty()->SetColor(GetRValue(m_OutlineColor)/255.0,GetGValue(m_OutlineColor)/255.0,GetBValue(m_OutlineColor)/255.0);
		cubeaxes->GetYAxesLinesProperty()->SetColor(GetRValue(m_OutlineColor)/255.0,GetGValue(m_OutlineColor)/255.0,GetBValue(m_OutlineColor)/255.0);
		cubeaxes->GetZAxesLinesProperty()->SetColor(GetRValue(m_OutlineColor)/255.0,GetGValue(m_OutlineColor)/255.0,GetBValue(m_OutlineColor)/255.0);
		cubeaxes->GetTitleTextProperty(0)->SetColor(GetRValue(Xcolor)/255.0,GetGValue(Xcolor)/255.0,GetBValue(Xcolor)/255.0);
		cubeaxes->GetTitleTextProperty(1)->SetColor(GetRValue(Ycolor)/255.0,GetGValue(Ycolor)/255.0,GetBValue(Ycolor)/255.0);
		cubeaxes->GetTitleTextProperty(2)->SetColor(GetRValue(Zcolor)/255.0,GetGValue(Zcolor)/255.0,GetBValue(Zcolor)/255.0);
		cubeaxes->GetLabelTextProperty(0)->SetColor(GetRValue(Xcolor)/255.0,GetGValue(Xcolor)/255.0,GetBValue(Xcolor)/255.0);
		cubeaxes->GetLabelTextProperty(1)->SetColor(GetRValue(Ycolor)/255.0,GetGValue(Ycolor)/255.0,GetBValue(Ycolor)/255.0);
		cubeaxes->GetLabelTextProperty(2)->SetColor(GetRValue(Zcolor)/255.0,GetGValue(Zcolor)/255.0,GetBValue(Zcolor)/255.0);
	}//else
	{
		/*for (int i = 0; i < 3; i++)
		{
			m_CubeAxes->GetTitleTextProperty(i)->SetFontSize(25);
			m_CubeAxes->GetTitleTextProperty(i)->SetFontFamily(2);
			m_CubeAxes->GetTitleTextProperty(i)->SetColor(0,0,0);
			m_CubeAxes->GetLabelTextProperty(i)->SetColor(0,0,0);
			m_CubeAxes->GetLabelTextProperty(i)->SetFontFamilyToTimes();
			m_CubeAxes->GetLabelTextProperty(i)->SetFontSize(25);
		}
		m_CubeAxes->SetXAxisTickVisibility(true);
		m_CubeAxes->SetYAxisTickVisibility(true);
		m_CubeAxes->SetTickLocationToOutside();
		m_CubeAxes->GetXAxesLinesProperty()->SetColor(0,0,0);
		m_CubeAxes->GetYAxesLinesProperty()->SetColor(0,0,0);*/
	}
	
	cubeaxes->GetXAxesLinesProperty()->SetLineWidth(3);
	cubeaxes->GetYAxesLinesProperty()->SetLineWidth(3);
	cubeaxes->GetZAxesLinesProperty()->SetLineWidth(3);
	//3. 设置边框面
	vtkSmartPointer<vtkCubeSource>cube = vtkSmartPointer<vtkCubeSource>::New();
	cube->SetBounds(Bounds);

	vtkSmartPointer<vtkPolyDataMapper>cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	cubeMapper->SetInputConnection(cube->GetOutputPort());

	BourderFaceActor->SetMapper(cubeMapper);
	BourderFaceActor->GetProperty()->SetColor(GetRValue(BourderColor)/255.0,GetGValue(BourderColor)/255.0,GetBValue(BourderColor)/255.0);
	BourderFaceActor->GetProperty()->SetOpacity(0.1);
	return 0;
}
int CVTKGraphics::Image(vtkStructuredGrid* sgrd,vtkActor* imageactor,vtkLookupTable* colortable)
{
	//Image
	vtkSmartPointer<vtkPlane> cutplane=vtkSmartPointer<vtkPlane>::New();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> plane=vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	int *extent=sgrd->GetExtent();
	double* bounds=sgrd->GetBounds();
	plane->SetExtent(extent[0],extent[1],extent[2],extent[3],extent[4],extent[4]);
	plane->SetInputData(sgrd);
	vtkSmartPointer<vtkPolyDataMapper> cutMapper=vtkSmartPointer<vtkPolyDataMapper>::New();
	cutMapper->SetInputConnection(plane->GetOutputPort());
	cutMapper->SetScalarRange(sgrd->GetScalarRange());
	imageactor->SetMapper(cutMapper);
	//imageactor->SetVisibility(true);

	return 1;
}
int CVTKGraphics::Contour(vtkStructuredGrid* sgrd,vtkActor* contour,vtkActor* fillcontour,vtkActor2D* contourlabels,vtkLookupTable* colortable)
{
	//Contour
	//int numberofcontour=200;
	//vtkSmartPointer<vtkPlane> cutplane=vtkSmartPointer<vtkPlane>::New();
	//vtkSmartPointer<vtkStructuredGridGeometryFilter> plane=vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	//plane->SetExtent(grdinfo.m_AxisBounds[0],grdinfo.m_AxisBounds[1],grdinfo.m_AxisBounds[2],grdinfo.m_AxisBounds[3],grdinfo.m_Height_data,grdinfo.m_Height_data);
	//plane->SetInputData(sgrd);
	//vtkSmartPointer<vtkContourFilter> contours =
	//vtkSmartPointer<vtkContourFilter>::New();
	//contours->SetInputConnection(plane->GetOutputPort());
	//contours->GenerateValues(numberofcontour, sgrd->GetScalarRange());
//	vtkSmartPointer<vtkPolyDataMapper> contourLineMapper =
//	vtkSmartPointer<vtkPolyDataMapper>::New();
//	contourLineMapper->SetInputConnection(contours->GetOutputPort());
//	contourLineMapper->SetScalarRange(sgrd->GetScalarRange());
	//contourLineMapper->ScalarVisibilityOff();
//	m_GMContourLineActor->SetMapper(contourLineMapper);
//	m_GMContourLineActor->SetVisibility(false);

	//FillContour
	vtkSmartPointer<vtkPlane> cutplane=vtkSmartPointer<vtkPlane>::New();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> plane=vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	double* bounds=sgrd->GetBounds();
	plane->SetExtent(bounds[0],bounds[1],bounds[2],bounds[3],0,0);
	plane->SetInputData(sgrd);
	vtkSmartPointer<vtkBandedPolyDataContourFilter> bandedContours =
	vtkSmartPointer<vtkBandedPolyDataContourFilter>::New();
	bandedContours->SetInputConnection(plane->GetOutputPort());
	bandedContours->SetScalarModeToValue();
	bandedContours->GenerateContourEdgesOn();
	bandedContours->GenerateValues( m_ContourNumber, sgrd->GetScalarRange());
	//根据数据范围建立色标
	//填充
	vtkSmartPointer<vtkPolyDataMapper> contourMapper =
	vtkSmartPointer<vtkPolyDataMapper>::New();
	contourMapper->SetInputConnection(bandedContours->GetOutputPort());
	contourMapper->SetScalarRange(sgrd->GetScalarRange());
	contourMapper->SetScalarModeToUseCellData();
	contourMapper->SetLookupTable(colortable);
	fillcontour->SetMapper(contourMapper);
	fillcontour->GetProperty()->SetInterpolationToFlat();
	//fillcontour->SetVisibility(true);
	//等值线
	vtkSmartPointer<vtkPolyDataMapper> contourLineMapper =
	vtkSmartPointer<vtkPolyDataMapper>::New();
	contourLineMapper->SetInputData(bandedContours->GetContourEdgesOutput());
	contourLineMapper->SetScalarRange(sgrd->GetScalarRange());
	//contourLineMapper->ScalarVisibilityOff();//等值线着色关闭
	contour->SetMapper(contourLineMapper);
	//contour->SetVisibility(false);
	//contourlabel
	vtkSmartPointer<vtkMaskPoints> mask=vtkSmartPointer<vtkMaskPoints>::New();
	mask->SetInputConnection(bandedContours->GetOutputPort());
	mask->SetOnRatio(2);
	mask->SetMaximumNumberOfPoints(m_ContourNumber);
	mask->RandomModeOn();
	vtkSmartPointer<vtkSelectVisiblePoints> visPts=vtkSmartPointer<vtkSelectVisiblePoints>::New();
	visPts->SetInputConnection(mask->GetOutputPort());
	vtkSmartPointer<vtkLabeledDataMapper> ldm=vtkSmartPointer<vtkLabeledDataMapper>::New();
	ldm->SetInputConnection(mask->GetOutputPort());
	ldm->SetLabelFormat("%.2f");
	ldm->SetLabelModeToLabelScalars();
	ldm->GetLabelTextProperty()->SetFontFamilyToTimes();
	ldm->GetLabelTextProperty()->SetFontSize(10);
	ldm->GetLabelTextProperty()->SetColor(0,0,0);
	contourlabels->SetMapper(ldm);
	//contourlabels->SetVisibility(true);

	return 1;
}
int CVTKGraphics:: Surfer3D(vtkStructuredGrid* sgrd,vtkLODActor* surferactor,vtkLookupTable* colortable)
{
	vtkSmartPointer<vtkPlane> cutplane=vtkSmartPointer<vtkPlane>::New();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> plane=vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	int *extent=sgrd->GetExtent();
	double* bounds=sgrd->GetBounds();
	plane->SetExtent(extent[0],extent[1],extent[2],extent[3],0,0);
	plane->SetInputData(sgrd);
	double m_ZLength=pow((fabs(bounds[1]-bounds[0])*fabs(bounds[2]-bounds[3])),0.5)/4.0;
	double Zscale;
	Zscale=(sgrd->GetScalarRange()[1]-sgrd->GetScalarRange()[0])/m_ZLength;
	vtkWarpScalar *warp=vtkWarpScalar::New();
	warp->SetInputConnection(plane->GetOutputPort()/*pGrid*/);
	warp->SetNormal(0,0,1);
	warp->UseNormalOn();
	warp->SetScaleFactor(1.0/Zscale);//调节z幅度大小的
	warp->Update();
	
	vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
	normals->SetInputConnection(warp->GetOutputPort());
	normals->SplittingOff();
	normals->SetFeatureAngle(m_FeatureAngle);
	vtkPolyDataMapper *demMapper=vtkPolyDataMapper::New();
	demMapper->SetInputConnection(normals->GetOutputPort());
	normals->Update();
	demMapper->SetScalarRange(sgrd->GetScalarRange());
	demMapper->SetLookupTable(colortable);
	//_VTK_SetGM3DColorBar();//色度条
	//-------------------------------------------------
	surferactor->SetMapper(demMapper);
	//surferactor->SetVisibility(true);
	warp->Delete();
	normals->Delete();
	demMapper->Delete();
	return 1;
}
int CVTKGraphics::Surfer3D(vtkStructuredGrid* sgrd, vtkLODActor* surferactor, vtkLookupTable* colortable, vector<POINT3D>* BoundingPoint)
{
	vtkSmartPointer<vtkPlane> cutplane = vtkSmartPointer<vtkPlane>::New();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> plane = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	int *extent = sgrd->GetExtent();
	double* bounds = sgrd->GetBounds();
	plane->SetExtent(extent[0], extent[1], extent[2], extent[3], 0, 0);
	plane->SetInputData(sgrd);
	double m_ZLength = pow((fabs(bounds[1] - bounds[0])*fabs(bounds[2] - bounds[3])), 0.5) / 4.0;
	double Zscale;
	Zscale = (sgrd->GetScalarRange()[1] - sgrd->GetScalarRange()[0]) / m_ZLength * 4;
	vtkWarpScalar *warp = vtkWarpScalar::New();
	warp->SetInputConnection(plane->GetOutputPort()/*pGrid*/);
	warp->SetNormal(0, 0, 1);
	warp->UseNormalOn();
	warp->SetScaleFactor(1.0);//调节z幅度大小的
	warp->Update();
	//----提取三维图的点的实际坐标----------------------------------------------------------------
	vtkPointSet* pointset = warp->GetOutput();
	//CString str; str.Format("%d  %d  %d  %d", pointset->GetNumberOfPoints(),extent[1],extent[3],extent[5]); AfxMessageBox(str);
	double* firstpoint = pointset->GetPoint(0);
	BoundingPoint[0].clear();
	BoundingPoint[2].clear();
	BoundingPoint[1].clear();
	BoundingPoint[3].clear();
	POINT3D point3d;
	int pointnumber = pointset->GetNumberOfPoints();
	int pointnumX = extent[1] * 2;//因为这个经过算法的数据点数与原始数据点数不同了，貌似比我们给定的网格节点要多，所以从我们需要找的这个列开始在增加一列进行判断搜索，确保正确和取完
	int num1 = pointnumber - pointnumX;
	//x轴平行的两个面
	for (size_t i = 0; i < pointnumX; i++)
	{
		firstpoint = pointset->GetPoint(i);
		point3d.x = firstpoint[0];
		point3d.y = firstpoint[1];
		point3d.z = firstpoint[2];
		if (point3d.y == bounds[2])
		{
			BoundingPoint[0].push_back(point3d);
		}
		firstpoint = pointset->GetPoint(num1+i);
		point3d.x = firstpoint[0];
		point3d.y = firstpoint[1];
		point3d.z = firstpoint[2];
		if (point3d.y == bounds[3])
		{
			BoundingPoint[2].push_back(point3d);
		}
	}
	//y轴平行的两个面
	for (size_t i = 0; i < pointnumber; i++)
	{
		firstpoint = pointset->GetPoint(i);
		point3d.x = firstpoint[0];
		point3d.y = firstpoint[1];
		point3d.z = firstpoint[2];
		if (point3d.x == bounds[1])
		{
			BoundingPoint[1].push_back(point3d);
		}
		if (point3d.x == bounds[0])
		{
			BoundingPoint[3].push_back(point3d);
		}
	}
	//--------------------------------------------------------------------------------------
	vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
	normals->SetInputConnection(warp->GetOutputPort());
	normals->SplittingOff();
	normals->SetFeatureAngle(m_FeatureAngle);
	vtkPolyDataMapper *demMapper = vtkPolyDataMapper::New();
	demMapper->SetInputConnection(normals->GetOutputPort());
	normals->Update();
	demMapper->SetScalarRange(sgrd->GetScalarRange());
	demMapper->SetLookupTable(colortable);
	//_VTK_SetGM3DColorBar();//色度条
	//-------------------------------------------------
	surferactor->SetMapper(demMapper);
	//surferactor->SetVisibility(true);
	warp->Delete();
	normals->Delete();
	demMapper->Delete();
	return 1;
}
int CVTKGraphics::SurfaceTexture(vtkStructuredGrid* sgrd, vtkLODActor* surferactor, vtkLookupTable* colortable, vector<POINT3D>* BoundingPoint,double* pos)
{
	vtkSmartPointer<vtkPlane> cutplane = vtkSmartPointer<vtkPlane>::New();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> plane = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	int *extent = sgrd->GetExtent();
	double* bounds = sgrd->GetBounds();
	plane->SetExtent(extent[0], extent[1], extent[2], extent[3], 0, 0);
	plane->SetInputData(sgrd);
	double m_ZLength = pow((fabs(bounds[1] - bounds[0])*fabs(bounds[2] - bounds[3])), 0.5) / 4.0;
	vtkSmartPointer<vtkWarpScalar> warp = vtkSmartPointer<vtkWarpScalar>::New();
	warp->SetInputConnection(plane->GetOutputPort()/*pGrid*/);
	warp->SetNormal(0, 0, 1);
	warp->UseNormalOn();
	warp->SetScaleFactor(1.0 );//调节z幅度大小的
	warp->Update();
	//2. 提取三维图的点的实际坐标----------------------------------------------------------------
	//CString str; str.Format("%d %d %d %d %d %d", extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]); AfxMessageBox(str);
	vtkPointSet* pointset = warp->GetOutput();
	double* firstpoint = pointset->GetPoint(0);
	/*FILE* fp = fopen("test3.txt", "w");
	if (fp == NULL)
	{
		AfxMessageBox("失败");
	}*/
	BoundingPoint[4].clear();
	BoundingPoint[6].clear();
	BoundingPoint[5].clear();
	BoundingPoint[7].clear();
	POINT3D point3d;
	int pointnumber = pointset->GetNumberOfPoints();
	int pointnum1 = extent[1] * 2;//因为这个经过算法的数据点数与原始数据点数不同了，貌似比我们给定的网格节点要多，所以从我们需要找的这个列开始在增加一列进行判断搜索，确保正确和取完
	int num1 = pointnumber - pointnum1;
	for (size_t i = 0; i <pointnum1; i++)
	{
		firstpoint = pointset->GetPoint(i);
		point3d.x = firstpoint[0];
		point3d.y = firstpoint[1];
		point3d.z = firstpoint[2];
		if (point3d.y == bounds[2])
		{
			BoundingPoint[4].push_back(point3d);
		}
		firstpoint = pointset->GetPoint(num1+i);
		point3d.x = firstpoint[0];
		point3d.y = firstpoint[1];
		point3d.z = firstpoint[2];
		if (point3d.y == bounds[3])
		{
			BoundingPoint[6].push_back(point3d);
		}
	}
	//y轴平行的两个面
	for (size_t i = 0; i < pointnumber; i++)
	{
			firstpoint = pointset->GetPoint(i);
			point3d.x = firstpoint[0];
			point3d.y = firstpoint[1];
			point3d.z = firstpoint[2];
			if (point3d.x == bounds[1])
			{
				BoundingPoint[5].push_back(point3d);
			}
			if (point3d.x == bounds[0])
			{
				BoundingPoint[7].push_back(point3d);
			}
	}
	//fclose(fp);
	//3. 根据船的位置和波浪的实际坐标大小更新船坐标的z值
		//计算与给定坐标最接近的网格节点(x0,y0)
	double* findpos = pointset->GetPoint(pointset->FindPoint(pos[0], pos[1], bounds[5]));
	pos[2] = findpos[2];
	//--------------------------------------------------------------------------------------
	vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputConnection(warp->GetOutputPort());
	//normals->SplittingOff();
	//normals->SetFeatureAngle(m_FeatureAngle);
	//--------------------------------------------------------------------------------------
	//vtkSmartPointer<vtkTextureMapToCylinder> tmapper = vtkSmartPointer<vtkTextureMapToCylinder>::New();
	//tmapper->SetInputConnection(normals->GetOutputPort());
	//vtkSmartPointer<vtkPNGReader> texture_reader = vtkSmartPointer<vtkPNGReader>::New();//载入图像文件数据
	//texture_reader->SetFileName("resources\\textures\\sea_foam.png");
	//vtkSmartPointer<vtkTexture> texture = vtkSmartPointer<vtkTexture>::New();
	//texture->SetInputConnection(texture_reader->GetOutputPort());
	//texture->InterpolateOn();
	//vtkSmartPointer<vtkTransformTextureCoords> xform = vtkSmartPointer<vtkTransformTextureCoords>::New();
	//xform->SetInputConnection(tmapper->GetOutputPort());
	//xform->SetScale(1, 1, 4);
	vtkSmartPointer<vtkPolyDataMapper> weightedTransMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//weightedTransMapper->SetInputConnection(xform->GetOutputPort()); 
	//------------------------------------------------------------------
	weightedTransMapper->SetInputConnection(normals->GetOutputPort());
	weightedTransMapper->SetScalarRange(sgrd->GetScalarRange());
	//----------------------海水颜色表----------------------------------------
	//CreateOceanWaterColorTable(colortable);
	//-------------------------------------------------------------------------
	weightedTransMapper->SetLookupTable(colortable);
	//-------------------------------------------------
	surferactor->SetMapper(weightedTransMapper);
	//surferactor->SetTexture(texture);
	surferactor->GetProperty()->SetOpacity(0.7);
	//surferactor->SetVisibility(true);
	//surferactor->GetProperty()->SetColor(27 / 255.0, 57 / 255.0, 109 / 255.0);
	//surferactor->GetProperty()->SetAmbient(100);
	surferactor->GetProperty()->SetDiffuse(0.0);
	//surferactor->GetProperty()->LightingOn();

	return 1;
}
int CVTKGraphics::CreateOceanWaterColorTable(vtkLookupTable* colorTable, vtkLookupTable* boundingcolor)
{
	colorTable->SetNumberOfColors(1);
	colorTable->SetTableValue(0, 27 / 255.0, 57 / 255.0, 109 / 255.0);
	//colorTable->SetTableValue(1, 44 / 255.0, 69 / 255.0, 106 / 255.0);
	//colorTable->SetTableValue(2, 84 / 255.0, 135 / 255.0, 172 / 255.0);
	colorTable->Build();

	//-----
	boundingcolor->SetNumberOfColors(255);
	boundingcolor->SetTableValue(0, 27 / 255.0, 57 / 255.0, 109 / 255.0);
	boundingcolor->SetTableValue(0, 27 / 255.0, 57 / 255.0, 255 / 255.0);
	boundingcolor->Build();
	return 1;
}

void CVTKGraphics::Plot3DLine(vector<POINT3D> xyzvector, vtkSmartPointer<vtkActor> lineactor, double linewidth)
{
	unsigned int nV = xyzvector.size();      // No. of vertices
	unsigned int nTv = 10;       // No. of surface elements for each tube vertex

	unsigned int i;
	//1.1三维坐标赋值
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (i = 0; i < nV; i++)
	{
		points->InsertPoint(i, xyzvector[i].x,xyzvector[i].y, xyzvector[i].z);
	}

	//1.2单元
	vtkSmartPointer<vtkCellArray> lines =
		vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(nV);
	for (i = 0; i < nV; i++)
	{
		lines->InsertCellPoint(i);
	}
	//1.3根据点和单元构建polydata
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);
	polyData->SetLines(lines);

	//粗细
	vtkSmartPointer<vtkDoubleArray> tubeRadius =
		vtkSmartPointer<vtkDoubleArray>::New();
	tubeRadius->SetName("Cable");
	tubeRadius->SetNumberOfTuples(nV);
	for (i = 0; i<nV; i++)
	{
		tubeRadius->SetTuple1(i, linewidth);
	}
	polyData->GetPointData()->AddArray(tubeRadius);
	polyData->GetPointData()->SetActiveScalars("Cable");

	//1.5过滤
	vtkSmartPointer<vtkTubeFilter> tube
		= vtkSmartPointer<vtkTubeFilter>::New();
	tube->SetNumberOfSides(nTv);
	tube->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
	
	//纹理
	vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputData(polyData);
	vtkSmartPointer<vtkTextureMapToCylinder> tmapper = vtkSmartPointer<vtkTextureMapToCylinder>::New();
	tmapper->SetInputConnection(normals->GetOutputPort());
	vtkSmartPointer<vtkPNGReader> texture_reader = vtkSmartPointer<vtkPNGReader>::New();//载入图像文件数据
	//vtkSmartPointer<vtkJPEGReader> texture_reader = vtkSmartPointer<vtkJPEGReader>::New();//载入图像文件数据
	//vtkSmartPointer<vtkBMPReader> texture_reader = vtkSmartPointer<vtkBMPReader>::New();//载入图像文件数据
	texture_reader->SetFileName("resources\\textures\\Sand.png");
	texture_reader->Update();
	vtkSmartPointer<vtkTexture> texture = vtkSmartPointer<vtkTexture>::New();
	texture->SetInputConnection(texture_reader->GetOutputPort());
	texture->InterpolateOn();
	vtkSmartPointer<vtkTransformTextureCoords> xform = vtkSmartPointer<vtkTransformTextureCoords>::New();
	xform->SetInputConnection(tmapper->GetOutputPort());
	//xform->SetScale(1, 1, 4);
	tube->SetInputConnection(xform->GetOutputPort());

	//1.6映射
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tube->GetOutputPort());
	mapper->ScalarVisibilityOn();
	mapper->SetScalarModeToUsePointFieldData();

	//1.7演员
	lineactor->SetMapper(mapper);
	lineactor->SetTexture(texture);
	
}
int CVTKGraphics::Export2PS(vtkSmartPointer<vtkWin32OpenGLRenderWindow> renwin, CString filename)
{
	vtkSmartPointer<vtkGL2PSExporter> exp = vtkSmartPointer<vtkGL2PSExporter>::New();
	exp->SetRenderWindow(renwin);
	exp->SetFileFormatToPS();
	//exp->UsePainterSettings();
	exp->CompressOff();
	exp->SetSortToSimple();
	//exp->DrawBackgroundOn();
	//exp->SetLineWidthFactor(1.0);
	//exp->SetPointSizeFactor(1.0);
	exp->SetFilePrefix((_bstr_t)filename);
	exp->Write3DPropsAsRasterImageOn();
	exp->Write();

	return 0;
}
int CVTKGraphics::Export2PS(vtkSmartPointer<vtkRenderWindow> renwin, CString filename)
{
	vtkSmartPointer<vtkGL2PSExporter> exp = vtkSmartPointer<vtkGL2PSExporter>::New();
	exp->SetRenderWindow(renwin);
	exp->SetFileFormatToPS();
	//exp->UsePainterSettings();
	exp->CompressOff();
	exp->SetSortToSimple();
	//exp->DrawBackgroundOn();
	//exp->SetLineWidthFactor(1.0);
	//exp->SetPointSizeFactor(1.0);
	exp->SetFilePrefix((_bstr_t)filename);
	exp->Write3DPropsAsRasterImageOn();
	exp->Write();

	return 0;
}

int CVTKGraphics::ChartXY(vtkContextView* pContextView, vtkTable* pTable, char* Title, int xcols, int type, int markstyle, float linewidth, unsigned char aerfa)
{
	//首先清空contextview的item
	pContextView->GetScene()->ClearItems();

	vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
	pContextView->GetScene()->AddItem(chart.GetPointer());
	chart->SetShowLegend(true);
	chart->SetTitle(Title);
	chart->GetTitleProperties()->SetFontFamilyToTimes();
	chart->GetTitleProperties()->SetFontSize(24);
	

	// Add multiple line plots, setting the colors etc
	vtkPlotLine *line;
	int colorindex = 0;
	for (int i = 1; i < pTable->GetNumberOfColumns(); i++)
	{
		line = vtkPlotLine::SafeDownCast(chart->AddPlot(type));
		line->SetInputData(pTable, xcols, i);
		line->SetColor(m_ColorArray[colorindex].r, m_ColorArray[colorindex].g, m_ColorArray[colorindex].b, aerfa);
		line->SetWidth(linewidth);
		line->SetMarkerStyle(i);
		colorindex++;
		if (colorindex >= m_ColorArray.size())
		{
			colorindex = 0;
		}
		//chart->SetPlotCorner(line, 0);//多条曲线同时现实，两个坐标轴
	}
	//坐标轴
	vtkAxis *left = chart->GetAxis(vtkAxis::LEFT);
	vtkAxis *bottom = chart->GetAxis(vtkAxis::BOTTOM);
	vtkChartLegend* legend = chart->GetLegend();
	left->SetTitle(pTable->GetColumn(1)->GetName());
	left->GetTitleProperties()->SetFontSize(20);
	left->GetTitleProperties()->SetColor(0, 0.0, 0.0);
	left->GetTitleProperties()->SetFontFamilyToTimes();
	left->GetLabelProperties()->SetFontFamilyToTimes();
	left->GetLabelProperties()->SetFontSize(20);
	left->GetLabelProperties()->SetColor(0, 0, 0);

	bottom->SetTitle(pTable->GetColumn(0)->GetName());
	bottom->GetTitleProperties()->SetFontSize(20);
	bottom->GetTitleProperties()->SetColor(0.0, 0.0, 0);
	bottom->GetTitleProperties()->SetFontFamilyToTimes();
	bottom->GetLabelProperties()->SetFontFamilyToTimes();
	bottom->GetLabelProperties()->SetFontSize(20);
	bottom->GetLabelProperties()->SetColor(0, 0, 0);
	
	legend->GetLabelProperties()->SetFontFamilyToTimes();
	legend->GetLabelProperties()->SetColor(0, 0, 0);
	legend->GetLabelProperties()->SetFontSize(15);

	
	// Set up a scientific plot...
	//chart->SetDrawAxesAtOrigin(true);
	//chart->SetShowLegend(true);
	//chart->GetAxis(vtkAxis::LEFT)->SetRange(1.0, -1.5);
	//chart->GetAxis(vtkAxis::LEFT)->SetNotation(2);
	//chart->GetAxis(vtkAxis::LEFT)->SetPrecision(1);
	//chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
	//chart->GetAxis(vtkAxis::BOTTOM)->SetRange(-1.0, 1.5);
	//chart->GetAxis(vtkAxis::BOTTOM)->SetNotation(2);
	//chart->GetAxis(vtkAxis::BOTTOM)->SetPrecision(1);
	//chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::FIXED);

	pContextView->Render();
	return 0;
}