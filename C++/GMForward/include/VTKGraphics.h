#pragma once

#define VTK_CREATE(type,name) \
	vtkSmartPointer<type> name = vtkSmartPointer<type>::New();
#include "GMDPS_proj.h"
//测试
#include "vtkTecplotReader.h"

//#include "vtkRenderingCore_AUTOINIT_vtkInteractionStyle_vtkRenderingFreeType_vtkRenderingFreeTypeOpenGL_vtkRenderingOpenGL.h"
#include "vtkLODActor.h"
#include "vtkLookupTable.h"
#include "vtkCubeAxesActor.h"
#include "vtkTextProperty.h"
#include "vtkActor.h"
#include "vtkPointSource.h"
#include "vtkScalarBarActor.h"
#include "vtkPolyDataNormals.h"
#include "vtkWin32OpenGLRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkAxesActor.h"
#include "vtkProperty.h"
#include "vtkInteractorStyleImage.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkTextMapper.h"
#include "vtkTextActor.h"
#include "vtkTextActor3D.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkOutlineFilter.h"
#include "vtkCellArray.h"
#include "vtkPolygon.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkCubeSource.h"
#include "vtkSphereSource.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkAppendPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkContourFilter.h"
#include "vtkMaskPoints.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkLabeledDataMapper.h"
#include "vtkWarpScalar.h"
#include "vtkThreshold.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridWriter.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkBandedPolyDataContourFilter.h"
#include "vtkDoubleArray.h"
#include "vtkSphereSource.h"
#include "vtkTextureMapToCylinder.h"
#include "vtkTransformTextureCoords.h"
#include "vtkPNGReader.h"
#include "vtkBMPReader.h"
#include "vtkPlaneSource.h"
#include "vtkTextureMapToPlane.h"
#include "vtkOBJReader.h"
#include "vtkTexturedSphereSource.h"
#include "vtkJPEGReader.h"
#include "vtkLight.h"
#include "vtkBoundingBox.h"
#include "vtk3DSImporter.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkCaptionActor2D.h"
#include "vtkAnnotatedCubeActor.h"
#include "vtkTubeFilter.h"
//#include "vtkGL2PSExporter.h"
#include "vtkContextView.h"
#include "vtkChartXY.h"
#include "vtkContextScene.h"
#include "vtkContextView.h"
#include "vtkFloatArray.h"
#include "vtkGL2PSExporter.h"
#include "vtkNew.h"
#include "vtkPlot.h"
#include "vtkPlotLine.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkTable.h"
#include "vtkContextActor.h"
#include "vtkRendererCollection.h"
#include "vtkAxis.h"
#include "vtkChartLegend.h"

#include "comutil.h "
#ifndef VTKGRAPHICS
#define VTKGRAPHICS

//颜色枚举
#define RED 255,0,0
#define GREEN 0,255,255
#define BLUE 0,0,255
#define YELLOW 255,255,0
#define MAGENTA 255,0,255
#define CYAN 0,255,255
#define WHITE 255,255,255
#define BLACK 0,0,0

//三维空间的点
struct POINT3D
{
	double x, y, z;
};
struct RGBColorStruct
{
	unsigned char r,g,b;
};
class CVTKGraphics
{
public:
	CVTKGraphics(void);
	~CVTKGraphics(void);
	int GrdToStructuredGridData(CGM_GRDDATA gmdata, vtkStructuredGrid* sgrd);
	int GrdToStructuredGridData(double** grddata, vtkStructuredGrid* sgrd,CGM_GRDDATAINFO grdinfo);
	int SaveGrd(vtkStructuredGrid* sgrd,CString filename);//保存为grd
	int Surfer3D(vtkStructuredGrid* sgrd,vtkLODActor* surferactor,vtkLookupTable* colortable);
	int Contour(vtkStructuredGrid* sgrd,vtkActor* contour,vtkActor* fillcontour,vtkActor2D* contourlabels,vtkLookupTable* colortable);
	int Image(vtkStructuredGrid* sgrd,vtkActor* imageactor,vtkLookupTable* colortable);
	int SetColorScale(vtkScalarBarActor* scalbar,char* title,double *pos,double width,double height,bool isHorizontal=true);
	int CreateColorTable(vtkLookupTable* colorTable, char* colorFile);
	int SetCubeAxesProp(vtkCubeAxesActor* cubeaxes,vtkActor* outlineActor,vtkActor* BourderFaceActor,double* bounds);
protected:
	int m_FeatureAngle;//三维显示数据的精细度
	int m_ContourNumber;//等值线个数
public:
	vector<RGBColorStruct> m_ColorArray;
	int SurfaceTexture(vtkStructuredGrid* sgrd, vtkLODActor* surferactor, vtkLookupTable* colortable, vector<POINT3D>* BoundingPoint,double* pos);
	int CreateOceanWaterColorTable(vtkLookupTable* colorTable, vtkLookupTable* boundingcolor);
	int Surfer3D(vtkStructuredGrid* sgrd, vtkLODActor* surferactor, vtkLookupTable* colortable, vector<POINT3D>* BoundingPoint);
	void Plot3DLine(vector<POINT3D> xyzvector, vtkSmartPointer<vtkActor> lineactor, double linewidth=0.4);
	int Export2PS(vtkSmartPointer<vtkWin32OpenGLRenderWindow>renwin, CString filename);
	int Export2PS(vtkSmartPointer<vtkRenderWindow> renwin, CString filename);
	int ChartXY(vtkContextView* pContextView, vtkTable* pTable,char* Title, int xcols = 0, int type = vtkChart::LINE, int markstyle = vtkPlotLine::CIRCLE, float linewidth = 2.0, unsigned char aerfa=255);
};

#endif
