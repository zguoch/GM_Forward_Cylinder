
function [Data_cut]=Signal_Cut_1D(Data_extension,k1,k2)
% Data_extension 表示扩边后的数据； k1,k2表示左、扩边的点数
% Data_cut 表示缩编后的数据

Length=length(Data_extension); 
NUM=Length-k1-k2; % 缩边后信号长度
Data_cut=zeros(1,NUM);
for i=1:Length
    if(i>=k1+1 && i<=Length-k2)
       Data_cut(i-k1)=Data_extension(i);
    end
end





