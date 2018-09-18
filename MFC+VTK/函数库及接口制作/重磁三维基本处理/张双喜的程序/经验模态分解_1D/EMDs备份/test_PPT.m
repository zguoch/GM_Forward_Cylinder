
clc
clear all
close all

N=1024/2;
n=0:N-1;
t=[0:1:1000];
%输入信号
%y=5*cos(2*pi*20*t)+4*sin(2*pi*40*t)+3*cos(pi*t);
y=sin(2*pi*30*t)+2*sin(2*pi*60*t);
figure;
%whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
plot(t,y);title('原始信号');

%FFT显示频谱
fs=10;%采样频率
X=fft(y,N);

% figure;
% plot(abs(X));
figure;
%whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
plot(n*fs/N,abs(X));title('原始信号的频谱图');

%EMD变换
[imf,ort,nbits]=emd(y);
imfsize=size(imf);
i=imfsize(1);

%各IMF分量及残余量
newrdxn=y;
figure;
for j=1:i
    subplot(i/2+1,2,j);plot(imf(j,:));%Ylabel('IMF');
     Ylabel(sprintf('IMF%d', j));
    newrdxn=newrdxn-imf(j,:);
end
% subplot(i/2+1,2,j+1);
% plot(newrdxn);Ylabel('残余量');

%每分解一次后的残余量
figure;
newrdxn=y;
for k=1:i
    newrdxn=newrdxn-imf(k,:);
    subplot(i/2+1,2,k);plot(newrdxn);title('残余量');
end

%HHT变换
[A,f,tt] = hhspectrum(imf(1:end-1,:));
[im,tt,ff] = toimage(A,f,tt,512);
disp_hhs(im,tt);
hgsave('HH_Spectrum');
%HHT边际谱
for kkk=1:512
    aaa(kkk)=sum(im(kkk,:));
end
 figure;
plot(ff(1,:),aaa);title('边际谱');
%figure;
%plot(n*fs/2/N,aaa);

