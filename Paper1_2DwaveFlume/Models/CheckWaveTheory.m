clear, close all

modelInfo = 'C:\Users\bknorris\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Scenarios\ModelScenarios.csv';
fid = fopen(modelInfo);
header = textscan(fgetl(fid),repmat('%s',1,4),'delimiter',',');
data = textscan(fid,repmat('%f',1,4),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

for i = 1:length(modelInfo.ScenarioNumber)
    h = modelInfo.h;
    Hs = modelInfo.Hs;
    Tp = modelInfo.Tp;
    SN = modelInfo.ScenarioNumber;
    
    %Show wave theory for each scenario
    fprintf('Scenario %d\n',SN(i))
    fprintf('Hs = %0.1f\n',Hs(i))
    fprintf('Tp = %0.1f\n',Tp(i))
    fprintf('h = %0.1f\n\n',h(i))
    
    f1 = NonLinearWaveSolver(Hs(i),Tp(i),h(i));
    ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[2200 -500   1650   1250]);
    pause()
    close(ff)
end