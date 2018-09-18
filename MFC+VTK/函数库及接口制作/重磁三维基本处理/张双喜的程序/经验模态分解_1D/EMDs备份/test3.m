clear
close all
clc
load spo2_0p2.txt;%spo2_3p.txt;%spo2_0p2.txt ;% spo2_3p.txt;spo2_0p4.txt spo2_3p.txt
result = spo2_0p2;
len=500;
x = [1:len];

%第二，第三个字节是红光脉搏波数据
rdxn = result([(x-1)*5+4])*256 +result([(x-1)*5+5]);

% figure;
% plot(abs(rdxn));


% 
[imf,ort,nbits]=emd(rdxn);
imfsize=size(imf);
i=imfsize(1);

%各IMF分量
figure;
for j=1:i
    subplot(i/2+1,2,j);plot(imf(j,:));
end
%每分解一次后的残余量
figure;
newrdxn=rdxn;
for k=1:i
    newrdxn=newrdxn-imf(k,:);
    subplot(i/2+1,2,k);plot(newrdxn);
    end
%HHT变换

[A,f,tt] = hhspectrum(imf(1:end-1,:));
% emd_visu(rdxn,x,imf)
% hgsave('all_imfs')
% 
%
im=toimage(A,f,tt);
disp_hhs(im,tt);
hgsave('HH_Spectrum')