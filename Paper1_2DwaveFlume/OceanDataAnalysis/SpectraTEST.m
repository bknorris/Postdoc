%%%TEST TEST TEST TEST TEST%%%%
clear

% load('MKK_HRaqd_Burst151_158_P.mat')
% load('MKK_RBRC15_Burst151_158_P.mat')
rng default

fs = 1000;
t = 0:1/fs:5-1/fs;
P = cos(2*pi*100*t) + randn(size(t));
%Begin spectra code to debug
nfft = 256;
dt = 4;

m = nfft;                                    %number of windows to average spectra, max(T) output will be m/2
noverlap = round(nfft/4,0);

%% The following lines from calcspectrum.m %%%
P = P(:);
n = max(size(P));	% Number of data points
index = 1:m;
w = hanning(m);		% Window specification; change this if you want:
k = fix((n-noverlap)/(m-noverlap));	% Number of windows
KMU = k*norm(w)^2;	% Normalizing scale factor
Pxx = zeros(m,1);
for i=1:k
    xw = w.*detrend(P(index));
    index = index + (m - noverlap);
    Xx = (abs(fft(xw,nfft)).^2)./KMU;
    Pxx = Pxx + Xx;
end
% Select first half and let 1st point (DC value) be replaced
% with 2nd point
select = [2 2:m/2];
Pxx = Pxx(select);
Pxx = Pxx*2/dt;
Cpp_i = Pxx;

f_i = (select - 1)'*dt/m;
figure(1)
plot(f_i,Cpp_i),hold on

[Cpp_m,f_m] = pwelch(P,hanning(nfft),round(nfft/4,0),nfft,dt);
plot(f_m,Cpp_m)
