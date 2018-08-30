clc;
clear;
%有限长水平圆柱体重力异常正演，延伸方向与y方向一致
G=6.67*10^-11;
R=0.4;%半径
den=1;%g/cm3
S=pi*R^2;
LinDen=den*S;%线密度
u0=4*pi*10^-7;
I=45;I=I/180*pi;
A=30;A=A/180*pi;
J=1/pi*R^2;
M=10;chl=M*u0/50000*10^9
Mx=M*cos(I)*sin(A);
My=M*cos(I)*cos(A);
Mz=M*sin(I);
L=100;%半长度
D=10;%中心埋深
y0=200;%中心水平位置y
x0=150;%中心水平位置x

x=100:1:200;
y=0:1:400;

for i=1:length(y)
    yy=y(i)-y0;
    for j=1:length(x)
        xx=x(j)-x0;
        %重力异常
        vz(i,j)=10^8*G*LinDen*D*((L-yy)/sqrt(xx^2+D^2+(L-yy)^2)+(L+yy)/sqrt(xx^2+D^2+(L+yy)^2))/(xx^2+D^2);
       L1=(L-yy);L2=(L+yy);
       b1=L1^2;b2=L2^2;
       temp1=sqrt(xx^2+D^2+L1^2);temp2=sqrt(xx^2+D^2+L2^2);
       v1=sign(L1)*log((temp1-sqrt(b1))/(temp1+sqrt(b1)));
       v2=sign(L2)*log((temp2-sqrt(b2))/(temp2+sqrt(b2)));
        v(i,j)=-10^8*G*LinDen*(v1+v2)/2;%重力位
        %水平方向导数vx
        vx1=sign(L1)*abs(L1)/temp1;
        vx2=sign(L2)*abs(L2)/temp2;
        vx(i,j)=-10^8*G*LinDen*xx/(xx^2+D^2)*(vx1+vx2);
        %水平方向导数vy
        vy(i,j)=-10^8*G*LinDen*(1/temp1-1/temp2);
        %二阶导数vxx
        vxx1=(L1)/temp1^3;
        vxx2=(L2)/temp2^3;
        vxx(i,j)=10^8*G*LinDen*((xx^2-D^2)/(xx^2+D^2)^2*(vx1+vx2)-xx^2/(xx^2+D^2)*(vxx1+vxx2));
        %二次导数vxy
        vxy(i,j)=10^8*G*LinDen*(xx/temp1^3-xx/temp2^3);
        %二次导数vxz
        vxz1=(L1/(temp1)+L2/(temp2))*2*xx/(xx^2+D^2)^2;
        vxz2=(L1/temp1^3+L2/temp2^3)*xx/(xx^2+D^2);
        vxz(i,j)=-10^8*G*LinDen*D*(vxz1+vxz2);
        %二次导数vyy
        vyy(i,j)=10^8*G*LinDen*((yy-L)/temp1^3-(L+yy)/temp2^3);
        %二次导数vyz
        vyz(i,j)=-10^8*G*LinDen*D*(1/temp1^3-1/temp2^3);
        %二次导数vzz
        vzz1=(xx^2-D^2)/(xx^2+D^2)^2*(L1/temp1+L2/temp2); 
        vzz2=D^2/(xx^2+D^2)*(L1/temp1^3+L2/temp2^3);
        vzz(i,j)=-10^8*G*LinDen*(vzz1-vzz2);
        %Hax
        haxx=1/(xx^2+D^2)*(sign(L1)*abs(L1)/temp1*((xx^2-D^2)/(xx^2+D^2)-xx^2/temp1^2)+...
            sign(L2)*abs(L2)/temp2*((xx^2-D^2)/(xx^2+D^2)-xx^2/temp2^2));
        haxy=xx*(1/temp1^3-1/temp2^3);
        haxz=-D*xx/(xx^2+D^2)*(L1/temp1*(2/(xx^2+D^2)+1/temp1^2)+...
            L2/temp2*(2/(xx^2+D^2)+1/temp2^2));
        Hax0(i,j)=u0/4/pi*(Mx*haxx+My*haxy+Mz*haxz);
        %Hay
        hayy=-(L1/temp1^3+L2/temp2^3);
        hayz=-D*(1/temp1^3-1/temp2^3);
        Hay0(i,j)=u0/4/pi*(Mx*haxy+My*hayy+Mz*hayz);
        %Za
        hazz=1/(xx^2+D^2)*(L1/temp1*((D^2-xx^2)/(xx^2+D^2)+D^2/temp1^2)+...
            L2/temp2*((D^2-xx^2)/(xx^2+D^2)+D^2/temp2^2));
        Za0(i,j)=u0/4/pi*(Mx*haxz+My*hayz+Mz*hazz);
    end
end

%计算磁异常

Hax=u0/4/pi/G/den*(Mx*vxx+My*vxy+Mz*vxz);
Hay=u0/4/pi/G/den*(Mx*vxy+My*vyy+Mz*vyz);
Za=u0/4/pi/G/den*(Mx*vxz+My*vyz+Mz*vzz);
Ta=Hax*cos(I)*cos(A)+Hay*cos(I)*sin(I)+Za*sin(I);
% figure(1);mesh(vz);
% figure(2);contour(v);
% figure(1);contour(vx);title('水平一次导数Vx');
% figure(1);contour(vy);title('水平一次导数Vy');
% figure(1);contour(vxx);title('二次导数Vxx');
% figure(1);contour(vxy);title('二次导数Vxy');
% figure(1);contour(vxz);title('二次导数Vxz');
% figure(1);contour(vyy);title('二次导数Vyy');
% figure(1);contour(vyz);title('二次导数Vyz');
% figure(1);contour(vzz);title('二次导数Vzz');
figure(1);mesh(vz);title('vz');
figure(2);contour(Hay);title('Hay');
figure(3);contour(Za);title('Za');
figure(4);mesh(Ta);title('Ta');
save_grd('有限长水平圆柱体vz.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vz)),max(max(vz)),vz);
save_grd('有限长水平圆柱体v.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(v)),max(max(v)),v);
save_grd('有限长水平圆柱体vx.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vx)),max(max(vx)),vx);
save_grd('有限长水平圆柱体vy.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vy)),max(max(vy)),vy);
save_grd('有限长水平圆柱体vxx.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vxx)),max(max(vxx)),vxx);
save_grd('有限长水平圆柱体vxy.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vxy)),max(max(vxy)),vxy);
save_grd('有限长水平圆柱体vxz.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vxz)),max(max(vxz)),vxz);
save_grd('有限长水平圆柱体vyy.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vyy)),max(max(vyy)),vyy);
save_grd('有限长水平圆柱体vyz.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vyz)),max(max(vyz)),vyz);
save_grd('有限长水平圆柱体vzz.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(vzz)),max(max(vzz)),vzz);
save_grd('有限长水平圆柱体Hax.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(Hax)),max(max(Hax)),Hax);
save_grd('有限长水平圆柱体Hay.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(Hay)),max(max(Hay)),Hay);
save_grd('有限长水平圆柱体Za.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(Za)),max(max(Za)),Za);
save_grd('有限长水平圆柱体Ta.grd',length(x),length(y),min(x),max(x),min(y),max(y),min(min(Za)),max(max(Za)),Za);







