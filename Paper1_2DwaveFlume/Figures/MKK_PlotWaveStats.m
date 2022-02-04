%Script to plot wave statistics from the C-line of the MKK deployment
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;

path = 'C:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
files = dir([path '*_waves.mat']);
fileList = {files.name};
%Load all files, store in temporary structure
for i = 1:length(fileList)
    inst = regexp(fileList{i},'(?<=MKK18)(.*?)([0-9]{2})','match');
    load([path fileList{i}])
    dat.(inst{:}) = wave;clear wave
end

f1 = figure(1);
set(f1,'PaperOrientation','landscape',...
    'position',[100 80   1050   550]);
fn = fieldnames(dat);
p = zeros(length(fn),1);
r = zeros(length(fn),1);
s = zeros(length(fn),1);
cc = brewermap(length(fn),'Spectral');
for i = 1:length(fn)
    sp(1) = subplot(3,1,1);
    hold on
    for j = 1:length(dat.(fn{i}).burst)
        p(i) = plot(dat.(fn{i}).time(:,j),dat.(fn{i}).depth(:,j),'color',cc(i,:));
    end
    hold off
    
    sp(2) = subplot(3,1,2);
    hold on
    for j = 1:length(dat.(fn{i}).burst)
        r(i) = plot(dat.(fn{i}).time(:,j),dat.(fn{i}).Tm(:,j),'color',cc(i,:));
    end
    hold off
    
    sp(3) = subplot(3,1,3);
    hold on
    for j = 1:length(dat.(fn{i}).burst)
        s(i) = plot(dat.(fn{i}).time(:,j),dat.(fn{i}).Hs(:,j),'color',cc(i,:));
    end
    hold off
end
leg = legend(p,fn);
set(leg,'position',[0.86 0.465 0.1 0.1])
linkaxes(sp,'x')
set([sp(1) sp(2)],'xticklabel',[])
set(sp(1),'position',[0.1 0.7 0.75 0.23])
set(sp(2),'position',[0.1 0.4 0.75 0.23])
set(sp(3),'position',[0.1 0.1 0.75 0.23])
ylabel(sp(1),'[m]'),title(sp(1),'Depth')
ylabel(sp(2),'[s]'),title(sp(2),'T_m')
ylabel(sp(3),'[m]'),title(sp(3),'H_s')
datetickzoom(sp(3),'x','mm-dd HH:MM:SS')


f2 = figure(2);
set(f2,'PaperOrientation','landscape',...
    'position',[100 80   650   450]);
fn = fieldnames(dat);
p = zeros(length(fn),1);
cc = brewermap(length(fn),'Spectral');
hold on;
for i = 1:length(fn)
    p(i) = plot(dat.(fn{i}).f(:,1,1),mean(dat.(fn{i}).Spp(:,1,:),3),'color',cc(i,:));
end
hold off;
set(gca,'xlim',[0 0.5])
ylabel('S_{\eta}(f) [m^2/{\Deltaf}]')
xlabel('f [Hz]')
leg = legend(p,fn);
set(leg,'box','off')
prettyfigures('text',12,'labels',13,'box',1,'tickdir','in','tlength',[0.005 0.005])


