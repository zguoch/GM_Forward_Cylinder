
%spo2_0p6.txt;%ecg1410.txt;%spo2_3p.txt;%spo2_0p2.txt ;% spo2_3p.txt;spo2_0p4.txt spo2_3p.txt
%result = spo2_0p2;%pecg1410.txt;
% len=500;
% x = [1:len];
clear
clc
close all

load eeg0410.txt;
rdxn=eeg0410;
% load eeg0410.mat;
% rdxn=eeg0410;
N=512;
n=0:N-1;
fs=100;
X=fft(rdxn,N);
figure;
plot(n*fs/N,abs(X));

%第二，第三个字节是红光脉搏波数据
%rdxn = result([(x-1)*5+4])*256 +result([(x-1)*5+5]);

figure;
plot(rdxn);
title('脑电图—EEG');

[imf,ort,nbits]=emd(rdxn);
imfsize=size(imf);
i=imfsize(1);

%各IMF分量
figure;
for j=1:i
    subplot(i/2+1,2,j);plot(imf(j,:));
end
%title('EMD分解');
%每分解一次后的残余量
% figure;
% newrdxn=rdxn;
% for k=1:i
%     newrdxn=newrdxn-imf(k,:);
%     subplot(i/2+1,2,k);plot(newrdxn);
% end
% title('残余量');


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
plot(ff(1,:),aaa);
%plot(aaa);
% figure;

%plot(n*100/N,aaa);