clc;
clear;
%倾斜的有限长水平圆柱体重力位
%
%*****************参数*******************
%----------磁性参数----------------------
K=0.2;                      %磁化率
Mag=50000;                  %地球磁场
I=45;                       %磁化倾角（度）
A=30;                       %磁化偏角（度）
I=I/180*pi;A=A/180*pi;
MagPar=[K,Mag,I,A];

%--------------------------------------------
Radius=0.4;
Length=200;
Density=1;
dx=5;dy=5;
xmin=0;xmax=300;
ymin=0;ymax=400;
x0=(xmax-xmin)/2;
y0=(ymax-ymin)/2;
z0=10;
AxisRange=[xmin,xmax,ymin,ymax];
CenterPosition=[x0,y0,z0];
DXDY=[dx,dy];
%----------------------------------------
x1=(xmax-xmin)/2;          %倾斜圆柱体两端点坐标
x2=(xmax-xmin)/2;
y1=(ymax-ymin)/2;
y2=(ymax-ymin)/2;
z1=10;
z2=60;
CylinderPos=[x1,x2,y1,y2,z1,z2];%圆柱体两端点位置
%--------------------------------------------------------------
[Ta,x,y]=Hax_Cylinder3D_Slop(Radius,Length,MagPar,CylinderPos,AxisRange,DXDY);
