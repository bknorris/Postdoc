%% Calculate Nikuradse Roughness length for roughWallFunctions
% The model requires an estimate the Nikuradse roughness lengthscale (Ks)
% to employ roughWallFunctions for U or K.
%
% This is Version 1 of this script
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;
dataPath = 'C:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
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

%Plot mean profile with minimum cell size (level 3 refinement)
ff2 = figure(2);
plot(Umean./Umean(15),zuv.*100,'-.k','linewidth',1.5)
set(gca,'ytick',0:1.5:max(zuv)*100,'ygrid','on',...
    'ylim',[0 30])
ylabel('\bf\itz (cm)')
xlabel('\bf\itU_z/U_{z = 0.13} (m/s)')

prettyfigures('text',11,'labels',12,'box',1)

if saveFigures == 1
    export_fig(ff1,[figPath  'z0_fits_burst' burstTxt],'-png','-r600','-nocrop')
    export_fig(ff2,[figPath  'Modelbins_velProfile_burst' burstTxt],'-png','-r600','-nocrop')
end