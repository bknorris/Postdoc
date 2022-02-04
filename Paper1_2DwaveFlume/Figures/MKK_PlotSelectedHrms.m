%Plot rms wave heights in swell and IG band for selected instruments from 
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
name = cell(1,length(insts));
cc = brewermap(2,'Accent');
ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   450]);
for i = 1:length(insts)
    load([dataPath insts{i}])
    sp(i) = subplot(length(insts),1,i);
    p = zeros(length(insts),1);
    hold on;
    p(1) = plot(wave.time,wave.HrmsSS,'LineWidth',1.5,'color',cc(1,:));
    p(2) = plot(wave.time,wave.HrmsIG,'LineWidth',1.5,'color',cc(2,:));
    name{i} = regexp(insts{i},'(?<=MKK18)(.*?)([0-9]{2})','match');
    title(name{i})
    legend(p,{'SS Band';'IG Band'});
end
ylabel(sp(1),'\bf\itH_{rms, SS, IG}')
ylabel(sp(2),'\bf\itH_{rms, SS, IG}')
datetick(sp(2),'x','dd-mm-yy HH:MM:SS')
set(sp,'ylim',[0 0.4],'ytick',0:0.1:0.4)
set(sp(1),'xticklabel',[])
prettyfigures('text',12,'labels',13,'box',1,'tickdir','out','tlength',[0.005 0.005])
export_fig(ff,[figPath 'HrmsSS_IG_C14_C16'],'-png')

