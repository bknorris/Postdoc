%Plot pressure spectra for selected instruments from the 2018 Molokai
%deployment
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'C:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Data\InstrumentData\Figures\Spectra\'; %where to save figs
insts = {'MKK18C1501rbr_waves.mat'}; %C14 is offshore of C15
p = zeros(1,length(insts));
name = cell(1,length(insts));
ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[100 80   550   450]);
for i = 1:length(insts)
    load([dataPath insts{i}])
    hold on;
    %     DOF = str2double(wave.atts.DOF);
    %     alfa = 1 - 0.95;
    %     c = chi2inv([1-alfa/2 alfa/2],DOF);
    %     c = DOF./c; %spectra confidence intervals
    %
    %     %Plot confidence intervals first
        F = mean(nanmean(wave.f(:,:,:),2),3);
        Spp = mean(nanmean(wave.Spp(:,:,:),3),2);
    %     area(F,Spp.*c(2),'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
    %     area(F,Spp.*c(1),'FaceColor',[1 1 1],'LineStyle','none')
    %     plot(F,Spp,'LineWidth',1.5,'color','k')
    %     name{i} = regexp(insts{i},'(?<=MKK18)(.*?)([0-9]{2})','match');
    %     title(name{i})
    p(i) = plot(F,Spp,'linewidth',1.5);    
end
ylabel('\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
xlabel('\bf\itf [Hz]')
set(gca,'xlim',[0 0.5])
leg = legend(p,{'C14';'C15'});
set(leg,'box','off')
prettyfigures('text',12,'labels',13,'box',1,'tickdir','out','tlength',[0.005 0.005])
export_fig(ff,[figPath 'SurfaceElevSpectra_C15'],'-png')
    
