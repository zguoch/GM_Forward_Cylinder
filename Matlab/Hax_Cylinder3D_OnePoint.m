function [ Hax] = Hax_Cylinder3D_OnePoint( Radius,Length,MagPar,CenterPosition,ForwardPos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   有限长水平圆柱体Hax正演函数(计算某一个单点）
%   输入参数：半径（m），长度（m），磁性参数(4×1）,中心位置 （x0,y0,D）,计算位置（x,y,z)
%   输出参数：重力位——二维数组（）
%   郭志馗，中国地质大学（武汉），2014.12.31
%   zhikuiguo@live.cn

u0=4*pi*10^-7;%真空中的磁导率
MagSusceptibility=MagPar(1);%磁化率(SI)，
Mag_Earth=MagPar(2);%地球磁场（nT），
I=MagPar(3);%磁化倾角（度）
A=MagPar(4);%磁化偏角（度）
R=Radius;%半径
S=pi*R^2;
M=S*MagSusceptibility*Mag_Earth;
Mx=M*cos(I)*sin(A);
My=M*cos(I)*cos(A);
Mz=M*sin(I);
L=Length/2;
D=CenterPosition(3);%中心埋深
y0=CenterPosition(2);%中心水平位置y
x0=CenterPosition(1);%中心水平位置x

xx=ForwardPos(1)-x0;
yy=ForwardPos(2)-y0;
 L1=(L-yy);L2=(L+yy);
temp1=sqrt(xx^2+D^2+L1^2);temp2=sqrt(xx^2+D^2+L2^2);
haxx=1/(xx^2+D^2)*(sign(L1)*abs(L1)/temp1*((xx^2-D^2)/(xx^2+D^2)-xx^2/temp1^2)+...
    sign(L2)*abs(L2)/temp2*((xx^2-D^2)/(xx^2+D^2)-xx^2/temp2^2));
haxy=xx*(1/temp1^3-1/temp2^3);
haxz=-D*xx/(xx^2+D^2)*(L1/temp1*(2/(xx^2+D^2)+1/temp1^2)+...
    L2/temp2*(2/(xx^2+D^2)+1/temp2^2));
Hax=u0/4/pi*(Mx*haxx+My*haxy+Mz*haxz);

end

