function [ Ta,x,y ] = Ta_Cylinder3D( Radius,Length,MagPar,CenterPosition,AxisRange,DXDY )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   Detailed explanation goes here
%   有限长水平圆柱体Hax正演函数
%   输入参数：半径（m），长度（m），磁性参数(4×1）,中心位置 （x0,y0,D）,坐标范围（xmin,xmax,ymin,ymax）,点距线距(dx,dy)
%   输出参数：重力位——二维数组（）
%   郭志馗，中国地质大学（武汉），2014.12.31
%   zhikuiguo@live.cn


I=MagPar(3);%磁化倾角（度）
A=MagPar(4);%磁化偏角（度）

[Hax,x,y]=Hax_Cylinder3D(Radius,Length,MagPar,CenterPosition,AxisRange,DXDY);

[Hay,x,y]=Hay_Cylinder3D(Radius,Length,MagPar,CenterPosition,AxisRange,DXDY);

[Za,x,y]=Za_Cylinder3D(Radius,Length,MagPar,CenterPosition,AxisRange,DXDY);

Ta=Hax*cos(I)*cos(A)+Hay*cos(I)*sin(I)+Za*sin(I);

end

