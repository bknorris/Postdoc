%Plot wave heights in swell and IG band for selected instruments from 
%the 2018 Molokai deployment
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'C:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\bknorris\Documents\Data\InstrumentData\Figures\Spectra\'; %where to save figs
insts = {'MKK18C1401rbr_waves.mat';'MKK18C1601rbr_waves.mat'}; %C14 is offshore of C16
sp = zeros(1,length(insts));
p = zeros(length(insts),1);
name = cell(1,length(insts));
cc = brewermap(2,'Paired');
ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   450]);
for i = 1:length(insts)
    load([dataPath insts{i}])
    sp(1) = subplot(2,1,1);
    
    hold on;
    p(i) = plot(wave.time,wave.HrmsSS.*sqrt(2),'LineWidth',1.5,'color',cc(i,:));
    
    sp(2) = subplot(2,1,2);
    hold on;
    plot(wave.time,wave.HrmsIG.*sqrt(2),'LineWidth',1.5,'color',cc(i,:))
    name{i} = regexp(insts{i},'(?<=MKK18)(.*?)([0-9]{2})','match');
end
leg = legend(p,[name{:}]);
set(leg,'position',[0.85 0.49 0.05 0.05])
ylabel(sp(1),'\bf\itH_{s, SS} [m]')
ylabel(sp(2),'\bf\itH_{s, IG} [m]')
datetick(sp(2),'x','dd-mm-yy HH:MM:SS')
set(sp(1),'ylim',[0 0.5],'ytick',0:0.15:0.5)
set(sp(2),'ylim',[0 0.15],'ytick',0:0.05:0.15)
set(sp(1),'xticklabel',[])
prettyfigures('text',12,'labels',13,'box',1,'tickdir','out','tlength',[0.005 0.005])
% export_fig(ff,[figPath 'Hs_SS_IG_C14_C16'],'-png')

