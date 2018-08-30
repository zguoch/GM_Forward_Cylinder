clc;
clear;
%多个有限长水平圆柱体磁异常叠加正演
%
%*****************参数*******************
%---------------磁性参数-----------------
K=0.2;                          %磁化率
Mag=50000*10^-9;                %地球磁场(T)
I=45;%磁化倾角（度）
A=30;%磁化偏角（度）
I=I/180*pi;A=A/180*pi;
MagPar=[K,Mag,I,A];
%------------几何参数---------------------
dx=5;dy=5;
xmin=0;xmax=300;
ymin=0;ymax=400;

Radius=0.4;
Length=200;
x0=(xmax-xmin)/2;
y0=(ymax-ymin)/2;
z0=10;
AxisRange=[xmin,xmax,ymin,ymax];
CenterPosition=[x0,y0,z0];
DXDY=[dx,dy];
%****************************************
[Ta,x,y]=Ta_Cylinder3D(Radius,Length,MagPar,CenterPosition,AxisRange,DXDY);
figure(1);pcolor(Ta);shading interp;colorbar;
figure(2);contour(Ta);title('Ta等值线');
figure(3);surf(Ta);colormap(hsv),shading interp;colorbar;