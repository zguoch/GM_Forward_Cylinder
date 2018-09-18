function value=SNR(sig1,sig2,m)
% calculate SNR
%
% value=SNR(sig1,sig2,m);
%
% sig1 - original signal
% sig2 - reconstructed signal
%

error(nargchk(2,3,nargin));

if (nargin<3) m='null'; end

sig1=reshape(sig1,1,prod(size(sig1)));
sig2=reshape(sig2,1,prod(size(sig2)));

switch m
 case 'rm'; sig1=remmean(sig1); sig2=remmean(sig2);
 case 'rv'; sig1=remmean(sig1); sig2=remmean(sig2); sig1=sig1/std(sig1); sig2=sig2/std(sig2);
end

value=20*log10(norm(sig1)/norm(sig1-sig2));
return

