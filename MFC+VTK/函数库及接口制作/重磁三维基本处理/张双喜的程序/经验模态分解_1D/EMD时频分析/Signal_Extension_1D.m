
function [Data_extension,k1,k2]=Signal_Extension_1D(Data,signal_length,K)
% data 表示原始信号；K 表示扩边后信号的长度
% Data_extension 表示返回的扩边信号
% k1,k2 分别表示返回的信号两边扩边的点数

Length=K-signal_length; % 需要延拓后信号长度
Data_extension=zeros(1,K);
if mod(Length,2)==0    % 若为偶数，左、右分配的点数相等
    k1=Length/2;
    k2=k1;
      for i=1:K
          if i<=(Length/2)
             Data_extension(i)=Data(Length/2-i+2);
          elseif (i>(Length/2) && i<=(signal_length+Length/2))
             Data_extension(i)=Data(i-Length/2);
          elseif (i>(signal_length+Length/2))
             Data_extension(i)=Data(signal_length-(i-(signal_length+Length/2)));
          end
      end 
else  % 若为奇数，左侧数值大，右侧数值小  
    k1=(Length+1)/2;
    k2=(Length-1)/2;
      for i=1:K
          if i<=((Length+1)/2)
             Data_extension(i)=Data((Length+1)/2-i+2);
          elseif (i>((Length+1)/2) && i<=(signal_length+(Length+1)/2))
             Data_extension(i)=Data(i-(Length+1)/2);
          elseif (i>(signal_length+(Length+1)/2))
             Data_extension(i)=Data(signal_length-(i-(signal_length+(Length+1)/2)));
          end
      end         
end