%%%=======绘制三维重力异常图===============================================
%{
M=800;
D=500;
G=6.67;
 x=-2000:10:2000;
y=-2000:10:2000;
x1=2000:10:6000;
y1=2000:10:6000;
N=401;
g=zeros(N);
g1=zeros(N);
for i=1:N
for j=1:N
g(i,j)=G*M*D./(x(i).^2+y(j).^2+D^2).^1.5;
g1(i,j)=G*M*D./(x(i).^2+y(j).^2+D^2).^1.5;
end
end
x0=zeros(N);
y0=zeros(N);
for i=1:N
    x0(i,:)=x;
    y0(:,i)=y';
end
figure(1);
mesh(x0,y0,g);
figure(2)
mesh(g1);
%}
%=======绘制叠加重力异常=========================================================

N=400;
D=80;%大圆柱体的埋深
d=25;%小圆柱体的埋深
R=30;%大圆柱体的半径
r=15;%小圆柱体的半径
X=2;%剩余密度
z1=X*pi*(R.^2);%单位长度大圆柱体的剩余质量
z2=X*pi*(r.^2);%单位长度小圆柱体的剩余质量
G=6.67;%万有引力常量
x=zeros(0,2*N+1);
for i=1:(2*N+1)
   x(i)=i-1-N;
end
M=zeros(0,2*N+1);
F=zeros(0,2*N+1);
K=zeros(0,2*N+1);

for j=1:(2*N+1)
    g0=2*G*z2*d/((x(j)+300).^2+d.^2);%从左向右依次表示
    g1=2*G*z1*D/((x(j)+200).^2+D.^2);
    g2=2*G*z2*d/((x(j)+150).^2+d.^2);
       g6=2*G*z2*d/((x(j)).^2+d.^2);
    g3=2*G*z2*d/((x(j)-150).^2+d.^2);
    g4=2*G*z1*D/((x(j)-200).^2+D.^2);
    g5=2*G*z2*d/((x(j)-300).^2+d.^2);
%     g7=2*G*z1*D/((x(j)-150).^2+D.^2);
%     g8=2*G*z1*D/((x(j)+150).^2+D.^2);
    M(j)=g0+g1+g2+g3+g4+g5+g6;  %总异常
    F(j)=g0+g2+g3+g5+g6;    %局部异常
    K(j)=g1+g4;        %区域异常
  
end
curve=plot(x,M,'r',x,F,'g',x,K,'k-');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
title('模型重力异常');
%legend('模型叠加异常','模型局部异常');
%set(lsline,'LineWidth',10);
 set(curve(1),'linewidth',3);%设置曲线1的粗细
 set(curve(2),'linewidth',1.5);%设置曲线1的粗细
 set(curve(3),'linewidth',1.5);%设置曲线1的粗细
 %legend('模型叠加异常','模型局部异常','模型区域异常');
 %%%%=======模型EMD分解===================================================== 
 
%  q=zeros(0,2*N+1);
%  for i=1:(2*N+1)
%   q(i)=q(0)+M(i);
%  end
[imf,ort,nbits]=emd(M); 
figure;
imfsize=size(imf);
i=imfsize(1);
for j=1:i
    subplot(i/2+1,2,j);plot(x,imf(j,:));
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end   
%}    
%================================等密度差圆柱体模型============================== 
%{
N=1500;
D=600;%大圆柱体的埋深
%D1=1000;%更深大圆柱体的埋深
d=150;%小圆柱体的埋深
%R1=80;%更深大圆柱体的半径
R=30;%大圆柱体的半径
r=15;%小圆柱体的半径
X=1000;%剩余密度
z1=X*pi*(R.^2);%单位长度大圆柱体的剩余质量
z2=X*pi*(r.^2);%单位长度小圆柱体的剩余质量
%z3=X*pi*(R1.^2);%单位长度更深大圆柱体的剩余质量
G=6.667*power(10,-11);%万有引力常量
x=zeros(1,2*N);
for i=1:(2*N)
   x(i)=i-1-N;
end
M=zeros(1,2*N);
F=zeros(1,2*N);
K1=zeros(1,2*N);
%Q=zeros(1,2*N);
for j=1:(2*N)
    g0=2*G*z2*r/((x(j)+1200).^2+d.^2);%从左向右依次表示
    g1=2*G*z1*R/((x(j)+800).^2+D.^2);
    g2=2*G*z2*r/((x(j)+500).^2+d.^2);
    g3=2*G*z2*r/((x(j)-500).^2+d.^2);
    g4=2*G*z1*R/((x(j)-800).^2+D.^2);
    g5=2*G*z2*r/((x(j)-1200).^2+d.^2);
   %g6=2*G*z3*R1/((x(j)).^2+D1.^2);
    M(j)=g0+g1+g2+g3+g4+g5;  %总异常
    F(j)=g0+g2+g3+g5;%局部异常
    K1(j)=g1+g4;%区域异常
    %Q(j)=g6;
end
figure
curve=plot(x,M*power(10,6),'r',x,F*power(10,6),'g',x,K1*power(10,6),'k');
%curve=plot(x,M,'r',x,F,'g',x,K,'k',x,Q,'c');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
title('模型重力异常');
%legend('模型叠加异常','模型局部异常');
%set(lsline,'LineWidth',10);
set(curve(1),'linewidth',3);%设置曲线1的粗细
set(curve(2),'linewidth',1.5);%设置曲线1的粗细
set(curve(3),'linewidth',1.5);%设置曲线1的粗细
%set(curve(4),'linewidth',1.5);%设置曲线1的粗细
 %legend('模型叠加异常','模型局部异常','模型区域异常');
%legend('模型叠加异常','模型局部异常','模型区域异常');
% rdxn=zeros(1,2*N);%获取行向量
% for i=1:(2*N)
%     rdxn(i)=M(i);
% end
[imf,ort,nbits]=emd(M);
%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);

for j=1:i
    
    subplot(i,1,j);curve=plot(x,imf(j,:)*power(10,6));set(curve(1),'linewidth',1.5);%设置曲线1的粗细
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end

%对IMF进行重构，观察重构误差========================
% sum=0;
% for j=1:i
%    sum= sum+imf(j,:);
% end
% figure
% diffence=sum-M;
% plot(x,diffence);title('重力重构误差');
% ylabel('重力重构误差值  Δg/g.u');xlabel('坐标值  X/m');   
%figure
%===================================计算|Ck|和|g|的幅值======================================================
MAX=M(find(diff(sign(diff(M)))==-2)+1)
%求极大值
MIN=M(find(diff(sign(diff(M)))==2)+1)
%求极小值
DOC=min(MIN);%求取最小值
%比较两个边界值与最小的大小
if(M(1)<DOC)
    DOC=M(1);
end
if(M(3000)<DOC)
    DOC=M(3000)
end
g=max(MAX)-DOC
%g=max(MAX)-min(MIN)
%%%%%==========叠加异常的幅度值======================
rexdn=zeros(1,3000);
rexdn=imf(1,:);
MAX=rexdn(find(diff(sign(diff(rexdn)))==-2)+1)
%求极大值
MIN=rexdn(find(diff(sign(diff(rexdn)))==2)+1)
%求极小值
doc=min(MIN);%求取最小值
%比较两个边界值与最小的大小
if(rexdn(1)<doc)
    doc=rexdn(1);
end
if(rexdn(3000)<doc)
    doc=rexdn(3000)
end

Ck=max(MAX)-doc
%Ck=max(MAX)-min(MIN)
%%%%%==========IMF分量的幅度值===============
figure
d=150;%球体模型，可根据对数功率谱计算要分离异常体的埋深，这里已知
%p=3*x^.2/(x.^2+d^.2);
%p=3*power(x,2)/(power(x,2)+power(d,2));%球体模型系数
p=2*power(x,2)/(power(x,2)+power(d,2));%圆柱体模型系数
t=Ck/g;
K=1.5/p;%水平方向的叠加，取系数1~2之间的数据
greg=imf(1,:)+K*t*imf(3,:);
subplot(211)
plot(x,F*power(10,6),'k',x,greg*power(10,6),'r-');
%axis([-1500,1500,1,10]);
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
legend('模型理论局部异常','EMD分解重构局部异常');
subplot(212)
diffence=F-greg;
plot(x,diffence*power(10,6));
title('重构误差');xlabel('坐标值 X/m');ylabel('误差值 Δg/g.u ');

figure
subplot(211)
en=imf(1,:)+imf(2,:)+imf(3,:)-greg;%重构区域异常
plot(x,K1*power(10,6),'k',x,en*power(10,6),'m');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
%axis([-1500,1500,1,6]);
legend('模型理论区域异常','EMD分解重构区域异常');
subplot(212)
Diffence=K1-en;
plot(x,Diffence*power(10,6));
title('重构误差');xlabel('坐标值 X/m');ylabel('误差值 Δg/g.u ');
figure
subplot(211)
curve=plot(x,greg*power(10,6));title('EMD分解重构局部重力异常');
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
subplot(212)
curve=plot(x,en*power(10,6));title('EMD分解重构区域重力异常');
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');

%==========径向对数功率谱========================================
figure
fs=1;%采样频率
N=3000; n=0:N-1; q=n*1/N;
%yy=greg;%局部重力异常
%yy=en;%区域重力异常
yy=M;%叠加重力异常
y=fft(yy,N);
R=real(y);%实部
I=imag(y);%虚部
E=power(R,2)+power(I,2);
Xn=log(E);
curve=plot(q,Xn,'b');axis([0,0.05,-50,-15]);
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
title('径向对数功率谱（Radial logarithmic power spectrum）');
ylabel('lnE(ω)');xlabel('ω');
%=======================================================================

%}
%====================正负密度差球体模型===================================
%{
N=2000;
D=500;%大球体的埋深
%D1=1000;%更深大圆柱体的埋深
d=50;%小球体的埋深  单位m
%R1=80;%更深大圆柱体的半径
R=50;%大球体的半径
r=35;%小球体的半径
X=1000;%大球剩余密度 单位g/cm3
X0=4000;%小球剩余密度   单位g/cm3
MM=(4/3)*pi*power(R,3)*X;%大球体的剩余质量 单位为kg
mm=(4/3)*pi*power(r,3)*X0;%小球体的剩余质量  单位为kg
% z1=X*pi*(R.^2);%单位长度大圆柱体的剩余质量
% z2=X*pi*(r.^2);%单位长度小圆柱体的剩余质量
%z3=X*pi*(R1.^2);%单位长度更深大圆柱体的剩余质量
G=6.667*power(10,-11);%万有引力常量
X1=-1000;%大球剩余密度
X10=-4000;%小球剩余密度  单位g/cm3
MM1=(4/3)*pi*power(R,3)*X1;%大球体的剩余质量
mm1=(4/3)*pi*power(r,3)*X10;%小球体的剩余质量
x=zeros(1,2*N);
for i=1:(2*N)
   x(i)=i-1-N;
end
M=zeros(1,2*N);
F=zeros(1,2*N);
K1=zeros(1,2*N);
%Q=zeros(1,2*N);
for j=1:(2*N)
    g0=G*mm*d/(((x(j)+1200).^2+d.^2)^.3/2);
    g1=G*MM*D/(((x(j)+800).^2+D.^2)^.3/2);
    g2=G*mm*d/(((x(j)+400).^2+d.^2)^.3/2);
    
%     g0=2*G*z2*r/((x(j)+1200).^2+d.^2);%从左向右依次表示
%     g1=2*G*z1*R/((x(j)+800).^2+D.^2);
%     g2=2*G*z2*r/((x(j)+500).^2+d.^2);
    %%%=================================================
           
    %%%=================================================
     g3=G*mm1*d/(((x(j)-400).^2+d.^2)^.3/2);
     g4=G*MM1*D/(((x(j)-800).^2+D.^2)^.3/2);
     g5=G*mm1*d/(((x(j)-1200).^2+d.^2)^.3/2);
%     g3=2*G*z2*r/((x(j)-500).^2+d.^2);
%     g4=2*G*z1*R/((x(j)-800).^2+D.^2);
%     g5=2*G*z2*r/((x(j)-1200).^2+d.^2);
   %g6=2*G*z3*R1/((x(j)).^2+D1.^2);
    M(j)=g0+g1+g2+g3+g4+g5;  %总异常
    F(j)=g0+g2+g3+g5;%局部异常
    K1(j)=g1+g4;%区域异常
    %Q(j)=g6;
end
figure
curve=plot(x,M,'r',x,F,'g',x,K1,'k');
%curve=plot(x,F,'r');
%curve=plot(x,M,'r',x,F,'g',x,K,'k',x,Q,'c');
ylabel('重力值Δg  m/s2');xlabel('坐标值  X/m');
title('模型重力异常');
legend('模型叠加异常','模型局部异常','模型局部异常');
%set(lsline,'LineWidth',10);
set(curve(1),'linewidth',3);%设置曲线1的粗细
%set(curve(2),'linewidth',1.5);%设置曲线1的粗细
%set(curve(3),'linewidth',1.5);%设置曲线1的粗细
%set(curve(4),'linewidth',1.5);%设置曲线1的粗细
 %legend('模型叠加异常','模型局部异常','模型区域异常');
%legend('模型叠加异常','模型局部异常','模型区域异常');
% rdxn=zeros(1,2*N);%获取行向量
% for i=1:(2*N)
%     rdxn(i)=M(i);
% end

[imf,ort,nbits]=emd(M);
%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);

for j=1:i
    
    subplot(i,1,j);curve=plot(x,imf(j,:));set(curve(1),'linewidth',1.5);%设置曲线1的粗细
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end

%对IMF进行重构，观察重构误差========================
% sum=0;
% for j=1:i
%    sum= sum+imf(j,:);
% end
% figure
% diffence=sum-M;
% plot(x,diffence);title('重力重构误差');
% ylabel('重力重构误差值  Δg/g.u');xlabel('坐标值  X/m');   
%figure
%===================================计算|Ck|和|g|的幅值======================================================

MAX=M(find(diff(sign(diff(M)))==-2)+1)
%求极大值
MIN=M(find(diff(sign(diff(M)))==2)+1)
%求极小值
DOC=min(MIN);%求取最小值
EXCEL=max(MAX);%求取最大值
if(M(1)>EXCEL)
    EXCEL=M(1);
end
if(M(4000)>EXCEL)
    EXCEL=M(4000)
end
%比较两个边界值与最小值的大小
if(M(1)<DOC)
    DOC=M(1);
end
if(M(4000)<DOC)
    DOC=M(4000)
end
%g=EXCEL-DOC
g=max(MAX)-min(MIN)
%上式为极大值为最大值，极小值为最小值的情况
%%%%%==========叠加异常的幅度值======================
rexdn=zeros(1,4000);
rexdn=imf(1,:);
MAX=rexdn(find(diff(sign(diff(rexdn)))==-2)+1)
%求极大值
MIN=rexdn(find(diff(sign(diff(rexdn)))==2)+1)
%求极小值
doc=min(MIN);%求取最小值
excel=max(MAX);%求取最大值
if(rexdn(1)>excel)
    excel=rexdn(1);
end
if(M(4000)>excel)
    excel=rexdn(4000)
end
%比较两个边界值与最小的大小
if(rexdn(1)<doc)
    doc=rexdn(1);
end
if(rexdn(4000)<doc)
    doc=rexdn(4000)
end

%Ck=excel-doc
Ck=max(MAX)-min(MIN)
%上式为极大值为最大值，极小值为最小值的情况
%%%%%==========IMF分量的幅度值===============
figure
d=50;%球体模型，可根据对数功率谱计算要分离异常体的埋深，这里已知
%p=3*x^.2/(x.^2+d^.2);
p=3*power(x,2)/(power(x,2)+power(d,2));%球体模型系数
%p=2*power(x,2)/(power(x,2)+power(d,2));%圆柱体模型系数
t=Ck/g;
K=4/p;%水平方向的叠加，取系数1~2之间的数据
greg=imf(1,:)+K*t*imf(2,:);
subplot(211)
plot(x,F,'k',x,greg,'r-');
%axis([-1500,1500,1,10]);
ylabel('重力值Δg  m/s2');xlabel('坐标值  X/m');
legend('模型理论局部异常','EMD分解重构局部异常');
subplot(212)
diffence=F-greg;
plot(x,diffence);
title('重构误差');xlabel('坐标值 X/m');ylabel('误差值 Δg/mGal ');

figure
subplot(211)
%en=imf(1,:)+imf(2,:)+imf(3,:)-greg;%重构区域异常
en=M-greg;
plot(x,K1,'k',x,en,'m');
ylabel('重力值Δg  m/s2');xlabel('坐标值  X/m');
%axis([-1500,1500,1,6]);
legend('模型理论区域异常','EMD分解重构区域异常');
subplot(212)
Diffence=K1-en;
plot(x,Diffence);
title('重构误差');xlabel('坐标值 X/m');ylabel('误差值Δg  m/s2 ');
figure
subplot(211)
curve=plot(x,greg);title('EMD分解重构局部重力异常');
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
ylabel('重力值Δg  m/s2');xlabel('坐标值  X/m');
subplot(212)
curve=plot(x,en);title('EMD分解重构区域重力异常');
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
ylabel('重力值Δg  m/s2');xlabel('坐标值  X/m');


%==========径向对数功率谱========================================
figure
fs=1;%采样频率
N=4000; n=0:N-1; q=n*1/N;
%yy=greg;%局部重力异常
%yy=en;%区域重力异常
yy=M;%叠加重力异常
y=fft(yy,N);
R=real(y);%实部
I=imag(y);%虚部
E=power(R,2)+power(I,2);
Xn=log(E);
curve=plot(q,Xn,'b');axis([0,0.05,-10,20]);
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
title('径向对数功率谱（Radial logarithmic power spectrum）');
ylabel('lnE(ω)');xlabel('ω');
%=======================================================================
%}

% for i=1:100
%     for j=1:100
%             E=power(R(i)(j),2)+power(I(i)(j),2);
%     end
% end
% E=power(R,2)+power(I,2)
%imshow(abs(y),[-1 5],'notruesize');
%

% Xn=log(E);
% curve=plot(Xn,'b');axis([0,0.05,-10,20]);
% set(curve(1),'linewidth',1.5);%设置曲线1的粗细
% title('径向对数功率谱（Radial logarithmic power spectrum）');
% ylabel('lnE(ω)');xlabel('ω');
%============================圆柱体和球体==============================================================
%{
N=1500;
D=600;%大圆柱体的埋深
%D1=1000;%更深大圆柱体的埋深
d=150;%小圆柱体的埋深
%R1=80;%更深大圆柱体的半径
R=30;%大圆柱体的半径
r=15;%小圆柱体的半径
X=1000;%剩余密度
z1=X*pi*(R.^2);%单位长度大圆柱体的剩余质量
z2=X*pi*(r.^2);%单位长度小圆柱体的剩余质量
%z3=X*pi*(R1.^2);%单位长度更深大圆柱体的剩余质量
G=6.667*power(10,-11);%万有引力常量
%==============================================================

d0=50;%小球体的埋深
r0=35;%小球体的半径
X0=4000;%小球剩余密度   单位kg/m3
mm0=(4/3)*pi*power(r,3)*X0;%小球体的剩余质量  单位为kg
%g2=G*mm0*d0/(((x(j)+400).^2+d0.^2)^.3/2);

%================================================================
x=zeros(1,2*N);
for i=1:(2*N)
   x(i)=i-1-N;
end
M=zeros(1,2*N);
F=zeros(1,2*N);
K1=zeros(1,2*N);
%Q=zeros(1,2*N);
for j=1:(2*N)
    g0=2*G*z2*r/((x(j)+1200).^2+d.^2);%从左向右依次表示
    g1=2*G*z1*R/((x(j)+800).^2+D.^2);
    g2=2*G*z2*r/((x(j)+500).^2+d.^2);
    %g3=2*G*z2*r/((x(j)-500).^2+d.^2);
    %g3=G*mm0*d0/(((x(j)-500).^2+d0.^2)^.3/2);
    g3
    g4=2*G*z1*R/((x(j)-800).^2+D.^2);
    %g5=2*G*z2*r/((x(j)-1200).^2+d.^2);
    %g5=G*mm0*d0/(((x(j)-1200).^2+d0.^2)^.3/2);
   %g6=2*G*z3*R1/((x(j)).^2+D1.^2);
    M(j)=g0+g1+g2+g3+g4+g5;  %总异常
    F(j)=g0+g2+g3+g5;%局部异常
    K1(j)=g1+g4;%区域异常
    %Q(j)=g6;
end
figure
curve=plot(x,M*power(10,6),'r',x,F*power(10,6),'g',x,K1*power(10,6),'k');
%curve=plot(x,M,'r',x,F,'g',x,K,'k',x,Q,'c');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
title('模型重力异常');
%legend('模型叠加异常','模型局部异常');
%set(lsline,'LineWidth',10);
set(curve(1),'linewidth',3);%设置曲线1的粗细
set(curve(2),'linewidth',1.5);%设置曲线1的粗细
set(curve(3),'linewidth',1.5);%设置曲线1的粗细
%set(curve(4),'linewidth',1.5);%设置曲线1的粗细
 %legend('模型叠加异常','模型局部异常','模型区域异常');
%legend('模型叠加异常','模型局部异常','模型区域异常');
% rdxn=zeros(1,2*N);%获取行向量
% for i=1:(2*N)
%     rdxn(i)=M(i);
% end
[imf,ort,nbits]=emd(M);
%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);

for j=1:i
    
    subplot(i,1,j);curve=plot(x,imf(j,:)*power(10,6));set(curve(1),'linewidth',1.5);%设置曲线1的粗细
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end

%对IMF进行重构，观察重构误差========================
% sum=0;
% for j=1:i
%    sum= sum+imf(j,:);
% end
% figure
% diffence=sum-M;
% plot(x,diffence);title('重力重构误差');
% ylabel('重力重构误差值  Δg/g.u');xlabel('坐标值  X/m');   
%figure
%===================================计算|Ck|和|g|的幅值======================================================
MAX=M(find(diff(sign(diff(M)))==-2)+1)
%求极大值
MIN=M(find(diff(sign(diff(M)))==2)+1)
%求极小值
DOC=min(MIN);%求取最小值
%比较两个边界值与最小的大小
if(M(1)<DOC)
    DOC=M(1);
end
if(M(3000)<DOC)
    DOC=M(3000)
end
g=max(MAX)-DOC
%g=max(MAX)-min(MIN)
%%%%%==========叠加异常的幅度值======================
rexdn=zeros(1,3000);
rexdn=imf(1,:);
MAX=rexdn(find(diff(sign(diff(rexdn)))==-2)+1)
%求极大值
MIN=rexdn(find(diff(sign(diff(rexdn)))==2)+1)
%求极小值
doc=min(MIN);%求取最小值
%比较两个边界值与最小的大小
if(rexdn(1)<doc)
    doc=rexdn(1);
end
if(rexdn(3000)<doc)
    doc=rexdn(3000)
end

Ck=max(MAX)-doc
%Ck=max(MAX)-min(MIN)
%%%%%==========IMF分量的幅度值===============
figure
d=150;%球体模型，可根据对数功率谱计算要分离异常体的埋深，这里已知
%p=3*x^.2/(x.^2+d^.2);
%p=3*power(x,2)/(power(x,2)+power(d,2));%球体模型系数
p=2*power(x,2)/(power(x,2)+power(d,2));%圆柱体模型系数
t=Ck/g;
K=1.5/p;%水平方向的叠加，取系数1~2之间的数据
greg=imf(1,:)+K*t*imf(3,:);
subplot(211)
plot(x,F*power(10,6),'k',x,greg*power(10,6),'r-');
%axis([-1500,1500,1,10]);
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
legend('模型理论局部异常','EMD分解重构局部异常');
subplot(212)
diffence=F-greg;
plot(x,diffence*power(10,6));
title('重构误差');xlabel('坐标值 X/m');ylabel('误差值 Δg/g.u ');

figure
subplot(211)
en=imf(1,:)+imf(2,:)+imf(3,:)-greg;%重构区域异常
plot(x,K1*power(10,6),'k',x,en*power(10,6),'m');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
%axis([-1500,1500,1,6]);
legend('模型理论区域异常','EMD分解重构区域异常');
subplot(212)
Diffence=K1-en;
plot(x,Diffence*power(10,6));
title('重构误差');xlabel('坐标值 X/m');ylabel('误差值 Δg/g.u ');
figure
subplot(211)
curve=plot(x,greg*power(10,6));title('EMD分解重构局部重力异常');
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
subplot(212)
curve=plot(x,en*power(10,6));title('EMD分解重构区域重力异常');
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');

%==========径向对数功率谱========================================
figure
fs=1;%采样频率
N=3000; n=0:N-1; q=n*1/N;
%yy=greg;%局部重力异常
%yy=en;%区域重力异常
yy=M;%叠加重力异常
y=fft(yy,N);
R=real(y);%实部
I=imag(y);%虚部
E=power(R,2)+power(I,2);
Xn=log(E);
curve=plot(q,Xn,'b');axis([0,0.05,-50,-15]);
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
title('径向对数功率谱（Radial logarithmic power spectrum）');
ylabel('lnE(ω)');xlabel('ω');
%}
%=========一维径向对数功率谱==============
data=load('modle1Za.dat');
ts=data(:,1);%读取第一列数据
rdxn=data(:,2);
fs=1.0/10;%采样频率
N=length(ts); n=0:N-1; q=n*fs/N;
%yy=greg;%局部重力异常
%yy=en;%区域重力异常
yy=rdxn;%叠加重力异常
y=fft(yy,N);
R=real(y);%实部
I=imag(y);%虚部
E=power(R,2)+power(I,2);
Xn=log(E);
figure
curve=plot(q,Xn,'b');axis([0,fs/2,4,24]);
set(curve(1),'linewidth',1.5);%设置曲线1的粗细
title('径向对数功率谱（Radial logarithmic power spectrum）');
ylabel('lnE(ω)');xlabel('ω');

%
%==========EMD用于二维重力异常的分离======================================
%{
N=1500;
D=600;%大圆柱体的埋深
%D1=1000;%更深大圆柱体的埋深
d=150;%小圆柱体的埋深
%R1=80;%更深大圆柱体的半径
R=30;%大圆柱体的半径
r=15;%小圆柱体的半径
X=1;%剩余密度
z1=X*pi*(R.^2);%单位长度大圆柱体的剩余质量
z2=X*pi*(r.^2);%单位长度小圆柱体的剩余质量
%z3=X*pi*(R1.^2);%单位长度更深大圆柱体的剩余质量
G=6.67;%万有引力常量
x=zeros(1,2*N);
for i=1:(2*N)
   x(i)=i-1-N;
end
M=zeros(1,2*N);
F=zeros(1,2*N);
K=zeros(1,2*N);
%Q=zeros(1,2*N);
for j=1:(2*N)
    g0=2*G*z2*r/((x(j)+1200).^2+d.^2);%从左向右依次表示
    g1=2*G*z1*R/((x(j)+800).^2+D.^2);
    g2=2*G*z2*r/((x(j)+500).^2+d.^2);
    
    g3=2*G*z2*r/((x(j)-500).^2+d.^2);
    g4=2*G*z1*R/((x(j)-800).^2+D.^2);
    g5=2*G*z2*r/((x(j)-1200).^2+d.^2);
   %g6=2*G*z3*R1/((x(j)).^2+D1.^2);
    M(j)=g0+g1+g2+g3+g4+g5;  %总异常
    F(j)=g0+g2+g3+g5;
    K(j)=g1+g4;
    %Q(j)=g6;
end
figure
curve=plot(x,M,'r',x,F,'g',x,K,'k');
%curve=plot(x,M,'r',x,F,'g',x,K,'k',x,Q,'c');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
title('模型重力异常');
%legend('模型叠加异常','模型局部异常');
%set(lsline,'LineWidth',10);
set(curve(1),'linewidth',3);%设置曲线1的粗细
set(curve(2),'linewidth',1.5);%设置曲线1的粗细
set(curve(3),'linewidth',1.5);%设置曲线1的粗细
%set(curve(4),'linewidth',1.5);%设置曲线1的粗细
 %legend('模型叠加异常','模型局部异常','模型区域异常');
%legend('模型叠加异常','模型局部异常','模型区域异常');
% rdxn=zeros(1,2*N);%获取行向量
% for i=1:(2*N)
%     rdxn(i)=M(i);
% end
[imf,ort,nbits]=emd(M);
%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);
for j=1:i
    subplot(i,1,j);plot(x,imf(j,:));
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end
%对IMF进行重构，观察重构误差========================
w=zeros(1,2*N);
w(1)=zeros(1,2*N);
for j=1:i
   w(j+1)= w(j)+imf(j,:)
end
diffence=w-M;
plot(x,diffence);title('重力重构误差');
ylabel('重力重构误差值  Δg/g.u');xlabel('坐标值  X/m');
%}



































