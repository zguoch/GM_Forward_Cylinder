
% HSPEC: Hilbert Amplitude Spectrum
%
% [S,freq]=hspec(imf,N);
%
% S   - Time-frequency-amplitude matrix
%       Columns are indexed in time, rows in frequency, values are amplitudes
%
% freq- instantaneous frequencies of each component
%
% imf - Matrix of intrinsic mode functions (each as a row)
%
% N   - Number of frequency cells
%
% See:  Huang et al, Royal Society Proceedings on Math, Physical, 
%       and Engineering Sciences, vol. 454, no. 1971, pp. 903-995, 
%       8 March 1998
%
% Remark: the graphical representation is the Hilbert Energy Spectrum
% that is: 20 * log ( S * S )
%
% Author: Ivan Magrin-Chagnolleau  <ivan@ieee.org>
% 

function [S,freq] = hspec(imf,N);


L = size(imf,1); % Number of components in the decomposition

%-------------------------------------------------------------------------
% loop for on each component

S = []; % Matrix which will contain the time-frequency-amplitude representation

      
clear x z m p freq
   
x = imf'; % now each column is a component
z = hilbert(x); % analytic signal
m = abs(z); % module of z
p = angle(z); % phase of z

for i = 1:L
   
   freq(:,i) = instfreq(z(:,i)); % instantaneous frequency
   
   % if the function instfreq is not available...
   % p(:,i) = unwrap(p(:,i)); % unwrap phase
   % freq(:,i) = abs(diff(p(:,i))); % derivative of the phase and absolute value
   % to have always positive frequencies
   
   ceilfreq(:,i) = ceil(freq(:,i)*N); % to have integer values - also do a smoothing given the number
   % of frequency cells
   
   for j = 1:length(x)-2
      
      S(ceilfreq(j,i),j+1) = m(j+1,i);
      
   end
      
end

eps = 0.00001; % to avoid zero values before the log
E = 20 * log ( S.^2 + eps ); % Hilbert energy spectrum

% plot S
figure;
t=1:length(x); % time sample
f=t/length(x)*0.5; % normalized frequency
imagesc(t,f,E); % !!! I am not sure it is the best way to visualize it !!!
colorbar;
set(gca,'YDir','normal');
xlabel('Time Sample');
ylabel('Normalized Frequency');

return
