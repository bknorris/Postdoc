% Calculate shallow water wave parameters (c and L)
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%  Wave period (s)
T = 30:30:300;
% Water depths vector (m)
h = 1;
% Gravitational acceleration constant (ms^-2)
g = 9.81;
% Shallow water wave phase speed $ c = sqrt{gh} $
% (ms^-1)
c = sqrt(g.*h);
%Shallow water wavelength
L = wavedispersion(T,h,'pade');
fprintf('Wave celerity: %5.1f m/s \n',c);
fprintf('Wavelength: %5.1f m \r\n',L);