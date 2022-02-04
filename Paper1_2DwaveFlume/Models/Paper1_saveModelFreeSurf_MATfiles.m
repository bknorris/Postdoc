% Load and save Paper1_2DwaveFlume model free surface data from the
% postProcessing folders in each model and save the data in MAT-format to
% disk.
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
versionNO = 'V1'; %set to same version as this script
saveData = 1; %set to 0 to turn off!

%% HR MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\';
modelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\';
modelInfo = 'HR_Domain_model_toRun.csv';

%HR MODEL Default settings:
waveGauges = [31 66 68 70 72 74 76];

%% LR MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\DataAnalysis\';
% workingDir = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\DataAnalysis\';
% % modelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Scenarios\';
% modelDir = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Scenarios\';
% modelInfo = 'LR_Domain_models_toRun.csv';
% 
% %LR MODEL Default settings:
% waveGauges = [31 63 65 67 69 71 73];

%% BEGIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
modelScenarios = dir(modelDir);
modelScenarios = modelScenarios(~ismember({modelScenarios.name},{'.','..'}));
dirFlags = [modelScenarios.isdir];
modelScenarios = modelScenarios(dirFlags); %process only directories!
modelScenarios = {modelScenarios.name};
count = 1;

%% Load and save model data!
for i = 1:length(modelScenarios)
    %Find model folders
    modelDataPath = [modelDir modelScenarios{i} '\MKKmodel\postProcessing\freeSurface\'];
    modelFolders = dir(modelDataPath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    modelFolders = {modelFolders.name};
%     scenarioNumber = regexp(modelScenarios{i},'\_(.*)','tokens');
    scenarioNumber = regexp(modelScenarios{i},'\_V(.*)','tokens');
    modelID = find(strcmp(string(modelInfo.scenarioNumber),scenarioNumber{:}));
    if isempty(modelID) %process only the folders listed in the CSV file
        continue
    else
        %Output to command window to track progress
        fprintf('Processing folder %0.0f of %0.0f: %s\n',count,length(modelInfo.scenarioNumber),modelScenarios{i})
        tic
        rampTime = modelInfo.rampTime(modelID);
        modelDuration = modelInfo.runTime(modelID);
        modelDepth = modelInfo.h(modelID);
        
        %Sort modelFolders
        sortModels = str2double(modelFolders); %windows sorts 1 then 10 then 100 then 2 then 20, etc.
        
        %Skip 0 time and rampTime
        [~,sortID] = sort(sortModels);
        skipRampTime = find(sortModels(sortID) > rampTime,1,'first');
        sortID = sortID(skipRampTime:end);
        
        %Model timesteps
        time = sortModels(sortID);
        
        %Preallocate data structure
        data = struct();
        data.time = time;
        data.waveGauges = waveGauges;
        
        %Load the data
        for j = 1:length(sortID)
            modelFiles = dir([modelDataPath modelFolders{sortID(j)} '\']);
            modelFiles = {modelFiles.name};
            
            %% LOAD FREE SURFACE DATA
            isFile = contains(modelFiles,'alpha.water_freeSurface');
            fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
            file = textscan(fid,'%n%n%n%n','headerlines',2);
            x = file{1};y = file{2};z = file{3};water = file{4};
            %Loop through wave gauges and save data to structure
            for k = 1:length(waveGauges)
                fieldName = sprintf('WG%0.0f',k);
                %Check to see if there are null values in the model data file
                if any((isempty(x) || isempty(y) || isempty(z) || isempty(water)))
                    fprintf('Timestep %s contains null values!\n',modelFolders{sortID(j)})
                    if j-1 ~= 0
                        data.(fieldName).eta(:,j) = NaN(length(data.(fieldName).eta(:,j-1)),1);
                    else
                        continue
                    end
                else
                    [ud,ix,iy]=uniquetol(x,1e-4);
                    output = [ud, accumarray(iy,z,[],@mean)];
                    [~,InstxID] = min(abs(output(:,1)-waveGauges(k)));
                    %Allocate cropped data to structure
                    data.(fieldName).eta(:,j) = output(InstxID,2);
                end
            end
            %Clear data to free memory
            clear x y z water file
            fclose(fid);
        end
        disp('Saving data to file')
        %Optionally save out data to disk
        if saveData == 1
            save([workingDir modelScenarios{i} '_freeSurf_' versionNO '.mat'],'-struct','data','-mat','-v7.3')
        end
        clear data
        disp(['Model data processing completed in ' num2str(toc/60) ' minutes'])
        count = count+1;
    end
end
