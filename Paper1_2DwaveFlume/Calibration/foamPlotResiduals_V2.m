% Plot residuals from an OpenFOAM run
%
% Use residuals calculated during runtime in controlDict
%
% This is Version 2 of this script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close
fdir = 'c:\Users\user\Documents\Models\Paper2_OptimizingRestoration\ModelRuns\Scenarios\Scenario_50\Model\postProcessing\residuals\0\';
figPath = 'c:\Users\user\Documents\Models\Figures\';
modelName = 'Scenario_42';
file = dir(fdir);file = file(~ismember({file.name},{'.','..'}));
fid = fopen([fdir file.name]);
line = fgetl(fid);header = fgetl(fid);
data = textscan(fid,'%n%s%n%n%n%n%n%n%s%s%n%n%n%s%s%n%n%n%s%s%n%n%n%s');
header = split(header);header = {header{2:25}};
residuals = cell2struct(data,header,2);

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[200 400   1200   450]);
hold on
fn = fieldnames(residuals);
toPlot = [3 4 6 7 11 12 16 17 21 22];
cc = brewermap(length(toPlot),'Paired');
legendTxt = cell(size(toPlot));
p = zeros(size(toPlot));
for i = 1:length(toPlot)
    p(i) = plot(residuals.Time,residuals.(fn{toPlot(i)}),'-','linewidth',1.5,'color',cc(i,:));
    legendTxt{i} = fn{toPlot(i)};
end
legendTxt = regexprep(legendTxt,'_',' ');
leg = legend(p,legendTxt);
set(leg,'position',[0.9 0.35 0.05 0.2])
set(gca,'xlim',[min(residuals.Time) max(residuals.Time)],...
    'ylim',[10^-14 10^0],...
    'yscale','log',...
    'position',[0.1 0.15 0.75 0.8])
xlabel('Model Time (s)')
ylabel('Residuals')
prettyfigures('text',12,'labels',13,'tickdir','out','tlength',[0.005 0.005])
hold off
% export_fig(ff,[figPath modelName '_residuals'],'-png','-r600','-nocrop')
