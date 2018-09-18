
clear;
clc;
% 读取文本数据文件

filename='model_two_TAHG.dat'; % 读取文件名
fname=filename(1:size(filename,2)-4);
data=load(filename);
ts=data(:,1); % 读取文本读取第1列
rdxn=data(:,2);

N=length(rdxn); % 信号的长度
D1=zeros(N,2); % 对第一阶高频开辟空间N行2列，多出的第一列用于存储时间数据
D2=zeros(N,2); % 第二阶高频
D3=zeros(N,2); % 第三阶高频
D4=zeros(N,2); % 第四阶高频
A=zeros(N,2); % 第四阶低频

% EMD分解
[imf,ort,nbits]=emd(rdxn);
imfsize=size(imf);
m=imfsize(1);% 分解阶次

figure(1);
plot(ts,rdxn,'k');title('原始异常');

figure(2);
newrdxn=rdxn';
for j=1:m
    subplot(m/2+1,2,j);plot(imf(j,:));% Ylabel('IMF');
     ylabel(sprintf('IMF%d', j));
    newrdxn=newrdxn-imf(j,:);
end

%{
% 重构异常
figure(3);
subplot(311)
plot(ts,rdxn,'b');
subplot(312)
plot(ts,imf(1,:)+imf(2,:)+imf(3,:),'b');
subplot(313)
plot(ts,imf(4,:)+imf(5,:),'b');
%}

for k=1:m
    t=imf(k,:);  
    for j=1:N
        NUM=k;
        switch NUM
            case 1
            D1(j,1)=ts(j);% 文件的第一列输出时间序列
            D1(j,2)=t(j);
            case 2
            D2(j,1)=ts(j);
            D2(j,2)=t(j);
            case 3
            D3(j,1)=ts(j);
            D3(j,2)=t(j);   
            case 4
            D4(j,1)=ts(j);
            D4(j,2)=t(j);  
            case 5
            A(j,1)=ts(j);
            A(j,2)=t(j);            
        end
    end
end


filename1=strcat(fname,'_EMD_1','.txt');
filename2=strcat(fname,'_EMD_2','.txt');
filename3=strcat(fname,'_EMD_3','.txt');
filename4=strcat(fname,'_EMD_4','.txt');
filename5=strcat(fname,'_EMD_5','.txt');

% 打开数据文件用于写入文件 
fidD1=fopen(filename1,'wt');
fidD2=fopen(filename2,'wt');
fidD3=fopen(filename3,'wt');
fidD4=fopen(filename4,'wt');
fidA=fopen(filename5,'wt');

% 向文本中存入数据
disp('文本存入进程：');
for j0=1:N
   disprog(j0,N,10); % 函数调用 10为显示步长
   for i0=1:2
       fprintf(fidD1,'%.4f   ',D1(j0,i0)); % 保留四位有效数字
       fprintf(fidD2,'%.4f   ',D2(j0,i0)); 
       fprintf(fidD3,'%.4f   ',D3(j0,i0));
       fprintf(fidD4,'%.4f   ',D4(j0,i0));
       fprintf(fidA,'%.4f   ',A(j0,i0)); 
   end
       % 每行数据输出后进行换行
       fprintf(fidD1,'\n');
       fprintf(fidD2,'\n');
       fprintf(fidD3,'\n');
       fprintf(fidD4,'\n');
       fprintf(fidA,'\n');
end
fclose(fidD1);
fclose(fidD2);
fclose(fidD3);
fclose(fidD4);
fclose(fidA);


% 图形显示
%{
figure ;
subplot(611)
plot(ts,rdxn);
subplot(612);
plot(ts,D1(:,2));%1阶高频
subplot(613);
plot(ts,D2(:,2));%2阶高频
subplot(614);
plot(ts,D3(:,2));%3阶高频
subplot(615);
plot(ts,D4(:,2));%4阶高频
subplot(616);
plot(ts,A(:,2));%4阶低频
%}

% 重构异常
%{
rdxn_new=imf(3,:)+imf(4,:)+imf(5,:);
%保存数据
filename7=strcat(fname,'_经验模态分解(EMD)_(3+4+5)','.txt');
fidD7=fopen(filename7,'wt');
 
for m=1:N
       fprintf(fidD7,'%.4f   %.6f',ts(m),rdxn_new(m));%保留四位有效数字
       fprintf(fidD7,'\n');
end      
fclose(fidD7);
%}

figure ;
subplot(211)
plot(ts,rdxn);title('原始信号');
subplot(212)
plot(ts,imf(4,:)+imf(5,:));title('去噪后的信号');

%{
figure ;
plot(ts,rdxn,'k',ts,rdxn_new,'r');title('原始信号');
%}











