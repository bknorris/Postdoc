% Load and save Paper1_2DwaveFlume model data from the postProcessing
% folder and save to disk.
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
saveData = 1; %set to 0 to turn off!

% HR MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\DataAnalysis\';
% modelInfo = 'HR_Domain_models.csv';

%HR MODEL Default settings:
% EDIT THESE!!

% LR MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\DataAnalysis\';
modelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Scenarios\';
modelInfo = 'LR_Domain_models.csv';

%LR MODEL Default settings:
loX = 62.75;hiX = 73.75;
loZ = -5;hiZ = [0 -0.25 -1 -1.5 -2]; %0.5, 1, 2, 3, 4 m models in that order
modelTimeStep = 0.125; %model timestep in seconds

%Load modelInfo and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,9),'delimiter',',');
data = textscan(fid,repmat('%f',1,9),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
modelScenarios = dir(modelDir);
modelScenarios = modelScenarios(~ismember({modelScenarios.name},{'.','..'}));
modelScenarios = {modelScenarios.name};

%% Load and save model data!
for i = 1:length(modelScenarios)
    %Output to command window to track progress
    fprintf('Processing folder %0.0f of %0.0f: %s\n',i,length(modelScenarios),modelScenarios{i})
    tic
    %Find model folders
    modelDataPath = [modelDir modelScenarios{i} '\MKKmodel\postProcessing\bathySample\surface\'];
    modelFolders = dir(modelDataPath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    modelFolders = {modelFolders.name};
    scenarioNumber = regexp(modelScenarios{i},'\_(.*)','tokens');
    modelID = find(strcmp(string(modelInfo.scenarioNumber),scenarioNumber{:}));
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
    
    %Find hiZ from modelInfo
    if modelDepth == 0.5
        zMax = hiZ(1);
    elseif modelDepth == 1
        zMax = hiZ(2);
    elseif modelDepth == 2
        zMax = hiZ(3);
    elseif modelDepth == 3
        zMax = hiZ(4);
    elseif modelDepth == 4
        zMax = hiZ(5);
    end
    
    %Preallocate data structure
    data = struct();
    
    %Load the data
    for j = 1:length(sortID)
        modelFiles = dir([modelDataPath modelFolders{sortID(j)} '\']);
        modelFiles = {modelFiles.name};
        
        %% LOAD U DATA
        isFile = contains(modelFiles,'U_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};Ux = file{4};Uy = file{5};Uz = file{6};
        %Crop the data file by the area defined in user settings
        cropX = find(x>=loX & x<=hiX);
        cropZ = find(z>=loZ & z<=zMax);
        cropXZ = intersect(cropX,cropZ);
        
        %Allocate cropped data to structure
        data.U.x(:,j) = x(cropXZ);
        data.U.y(:,j) = y(cropXZ);
        data.U.z(:,j) = z(cropXZ);
        data.U.Ux(:,j) = Ux(cropXZ);
        data.U.Uy(:,j) = Uy(cropXZ);
        data.U.Uz(:,j) = Uz(cropXZ);
        
        %Clear data to free memory
        clear x y z Ux Uy Uz file
        fclose(fid);
        
        %% LOAD p_rgh DATA
        isFile = contains(modelFiles,'p_rgh_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};p_rgh = file{4};
        
        %Allocate cropped data to structure
        data.p_rgh.x(:,j) = x(cropXZ);
        data.p_rgh.y(:,j) = y(cropXZ);
        data.p_rgh.z(:,j) = z(cropXZ);
        data.p_rgh.p_rgh(:,j) = p_rgh(cropXZ);
        
        %Clear data to free memory
        clear x y z p_rgh file
        fclose(fid);
        
        %% LOAD k DATA
        isFile = contains(modelFiles,'k_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};k = file{4};
        
        %Allocate cropped data to structure
        data.k.x(:,j) = x(cropXZ);
        data.k.y(:,j) = y(cropXZ);
        data.k.z(:,j) = z(cropXZ);
        data.k.k(:,j) = k(cropXZ);
        
        %Clear data to free memory
        clear x y z p_rgh file
        fclose(fid);
        
        %% LOAD k DATA
        isFile = contains(modelFiles,'k_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};k = file{4};
        
        %Allocate cropped data to structure
        data.k.x(:,j) = x(cropXZ);
        data.k.y(:,j) = y(cropXZ);
        data.k.z(:,j) = z(cropXZ);
        data.k.k(:,j) = k(cropXZ);
        
        %Clear data to free memory
        clear x y z p_rgh file
        fclose(fid);
    end
    disp('Saving data to file')
    %Optionally save out data to disk
    if saveData == 1
        save([workingDir modelScenarios{i} '_rawData.mat'],'data','-mat','-v7.3')
    end
    clear data
    disp(['Model data processing completed in ' num2str(toc/60) ' minutes'])
end
