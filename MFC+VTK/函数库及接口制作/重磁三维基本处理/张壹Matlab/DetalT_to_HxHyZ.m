clear;
clc;

filename='立方体磁异常';
[h,l,ymin,ymax,xmin,xmax,zmin,zmax,z,dx,dy]=opengrd(filename);%打开DeltaT文件

I0=45;%给出地磁场参数
A0=0;

M0=cosd(I0)*cosd(A0);%计算地磁场方向余弦
L0=cosd(I0)*sind(A0);
N0=sind(I0);

L=2^(floor(log2(l))+2);%扩边到2^n个数据
H=2^(floor(log2(h))+2);

[T,k1,k2,k3,k4]=kuobian_duoxiangshi(z,H,L,2,2,2);%对文件按需扩边
T=z;
H=h;
L=l;
fftT=fft2(T);%转换频谱

[U,V]=cal_UV(H,L,dx,dy);%计算波数
%{
fftHx=(U./((L0*U+M0*V)-1i*N0*sqrt(U.^2+V.^2)-1i*eps)).*fftT;
fftHy=(V./((L0*U+M0*V)-1i*N0*sqrt(U.^2+V.^2)-1i*eps)).*fftT;
fftZ =(sqrt(U.^2+V.^2)./(1i*(L0*U+M0*V)+N0*sqrt(U.^2+V.^2)+eps)).*fftT;
%}

% fftHx=(U./((L0*U+M0*V)-1i*N0*sqrt(U.^2+V.^2)-1i*eps)).*fftT;
fftHx=((N0*sqrt(U.^2+V.^2)*1i+(L0*U+M0*V)).*U./((N0^2*(U.^2+V.^2))+(L0*U+M0*V).^2+eps)).*fftT;
fftHy=(V./((L0*U+M0*V)-1i*N0*sqrt(U.^2+V.^2)-1i*eps)).*fftT;
fftZ =(sqrt(U.^2+V.^2)./(1i*(L0*U+M0*V)+N0*sqrt(U.^2+V.^2)+eps)).*fftT;

Hx=real(ifft2(fftHx));
Hy=real(ifft2(fftHy));
Z =real(ifft2(fftZ));

hx=Hx;
hy=Hy;
z=Z;
% hx=zeros(h,l);
% hy=zeros(h,l);
% z =zeros(h,l);
% for i=1:h
%     for j=1:l
%         hx(i,j)=Hx(i+k4,j+k1);
%         hy(i,j)=Hy(i+k4,j+k1);
%         z(i,j) =Z(i+k4,j+k1);
%     end
% end

filenamex=strcat(filename,'_Hx');
filenamey=strcat(filename,'_Hy');
filenamez=strcat(filename,'_Z');

savegrd(xmin,xmax,ymin,ymax,hx,filenamex);
savegrd(xmin,xmax,ymin,ymax,hy,filenamey);
savegrd(xmin,xmax,ymin,ymax,z,filenamez);
contour(hx);







