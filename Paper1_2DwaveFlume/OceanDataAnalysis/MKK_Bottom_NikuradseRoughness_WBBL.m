%% Calculate Nikuradse Roughness length for roughWallFunctions
% The model requires an estimate the Nikuradse roughness lengthscale (Ks)
% to estimate the height of the Wave Bottom Boundary Layer (WBBL)
%
%
%
% This is Version 2 of this script
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;
dataPath = 'd:\USGS\Data\InstrumentData\WavesData\';
figPath = 'd:\USGS\Documents\Models\Figures\'; %where to save figs
load([dataPath 'MKK18HR101aqdHR_waves.mat']); %'MKK18LR101aqdHR_waves.mat'
saveFigures = 0;

%% Copy Bursts from MKK_Inputs_for_irregularWaveModel.m
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
if length(AQDbursts) == 1
    burstTxt = sprintf('%d',AQDbursts(1));
else
    burstTxt = sprintf('%d-%d',AQDbursts(1),AQDbursts(2));
end
%% Calculate Ks
idx = find(wave.burst >= AQDbursts(1) & wave.burst <= AQDbursts(2));

%Extract rotated wave velocities from waveStats. Note, discard lower 3 bins (these
%are obscured by the bed)
Umean = mean(wave.Urot(1:end-3,idx),2);
zuv = wave.zuv(1:end-3);
bins = 8:14;

ff1 = figure(1);
xs = Umean./Umean(15);ys = log(zuv);
plot(xs,ys,'.k','linewidth',1.5),hold on
pf = polyfit(xs(bins),ys(bins),1);
pv = polyval(pf,xs(bins));
plot(xs(bins),ys(bins),'.r','linewidth',1.5),hold on
plot(xs(bins),pv,'r-','linewidth',1.5)
ylabel('\bf\itln(z) (m)')
xlabel('\bf\itU_z/U_{z = 0.13} (m/s)')
% eqn = sprintf('u(z) = %0.1dln(z)-%0.1f',pf(1),pf(2));
% tt(1) = text(0.6,0.45,eqn,'units','normalized');
z0 = exp(-pf(2)/pf(1));
Kn = 30*z0;
X = ['Kn= ',num2str(Kn)];
display(X);

%Calculate WBBL height
time = wave.time;
Aw = nanmean(wave.Uw./wave.omegar,1); %wave orbital excursion

Wdel = Kn*(0.09*(Aw./Kn).^0.82); %eqn. from Fredsoe & Deigaard, 1992
C4 = [0.072 0.896 0.111];
C5 = [-0.25 -0.469 -0.246];

Wdel2 = zeros(3,length(Aw));
for i = 1:3
    Wdel2(i,:) = Aw.*C4(i).*(Aw./Kn).^(C5(i));
end

ff2 = figure(2);
c = lines(4);
p = zeros(4,1);
p(1) = plot(time,Wdel.*1000,'color',c(1,:));hold on
for i = 1:3
    p(i+1) = plot(time,Wdel2(i,:).*1000,'color',c(i+1,:));
end
text = {'F & D';'Jonsson';'Sleath';'JSF'};
legend(p,text)
datetickzoom('x','dd HH:MM:SS')
ylabel('\delta_w (mm)')
title('Wave Boundary Layer Estimate')
