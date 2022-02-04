%Script to plot wave statistics from the C-line of the MKK deployment
% These data will be used for decision-making regarding the model domain
% development
%
% Updates:
% 01/14/21: added the ability to select a burst so the code will display
% the relevant statistics (instead of needing to copy these from the
% figure)
%
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;
dataPath = 'C:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\bknorris\Documents\Data\InstrumentData\Figures\Spectra\'; %where to save figs
insts = {'MKK18C1501rbr_waves.mat';'MKK18HR101aqdHR_waves.mat '};
startTime = datenum('24-Jun-2018 20:17:04');
endTime = datenum('25-Jun-2018 20:17:04');
whichBurst = [144 54]; %determine with MKK_EvolutionarySpectrum_V2.m

%Fields to display
burstNo = zeros(size(insts));
depths = zeros(size(insts));
Hs = zeros(size(insts));
Tm = zeros(size(insts));

%Plot Routine
sp = zeros(1,length(insts));
p = zeros(length(insts),1);
name = cell(1,length(insts));
cc = brewermap(2,'Paired');
ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   500]);
for i = 1:length(insts)
    load([dataPath insts{i}])
    idx = find(wave.time>=startTime & wave.time<=endTime);
    burstID = find(wave.burst == whichBurst(i));
    
    sp(1) = subplot(4,1,1);
    hold on;
    plot(wave.time(idx),wave.burst(idx),'LineWidth',1.5,'color',cc(i,:))
    burstNo(i) = wave.burst(burstID);
    
    sp(2) = subplot(4,1,2);
    hold on;
    plot(wave.time(idx),wave.depth(idx),'LineWidth',1.5,'color',cc(i,:))
    depths(i) = wave.depth(burstID);
    
    sp(3) = subplot(4,1,3);
    hold on;
    p(i) = plot(wave.time(idx),wave.Hs(idx),'LineWidth',1.5,'color',cc(i,:));
    Hs(i) = wave.Hs(burstID);
    
    sp(4) = subplot(4,1,4);
    hold on;
    plot(wave.time(idx),wave.Tm(idx),'LineWidth',1.5,'color',cc(i,:))
    name{i} = regexp(insts{i},'(?<=MKK18)(.*?)([0-9]{2})','match');
    Tm(i) = wave.Tm(burstID);
end
leg = legend(p,[name{:}]);
ylabel(sp(1),'\bf\itBurst No.')
ylabel(sp(2),'\bf\itDepth [m]')
ylabel(sp(3),'\bf\itH_{s} [m]')
ylabel(sp(4),'\bf\itT_{m} [s]')
set([sp(1) sp(2) sp(3)],'xticklabel',[])
linkaxes(sp,'x')
datetickzoom(sp(4),'x','dd/mm HH:MM:SS')
xlabel(sp(4),'Day in 2018')
set(leg,'position',[0.88 0.6 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1,'tickdir','out','tlength',[0.005 0.005])

%% Display model settings here!
for i = 1:length(insts)
fprintf('Instrument: %s\n',char(name{i}))
fprintf('Burst: %0.0f \n Avg. Depth: %0.2f m \n Hs: %0.1f m \n Tm: %0.1f s\n\n',burstNo(i),depths(i),Hs(i),Tm(i))
end