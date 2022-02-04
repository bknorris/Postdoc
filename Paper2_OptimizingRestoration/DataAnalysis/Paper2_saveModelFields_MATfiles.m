% Load and save Paper2_OptimizingRestoration model field data from the postProcessing
% folders in each model and save the data in MAT-format to disk.
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

%% MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%
workingDir = 'c:\Users\user\Documents\Models\Paper2_OptimizingRestoration\ModelRuns\CPUTests\DataAnalysis\';
modelDir = 'c:\Users\user\Documents\Models\Paper2_OptimizingRestoration\ModelRuns\CPUTests\';
modelInfo = 'CPUTest_models_toRun.csv';

%HR MODEL Default settings:
loX = 0;hiX = 2.6;
loZ = -5;hiZ = [-0.005 -0.005 -0.005]; %Levels 1, 2, 3 in that order
modelTimeStep = 0.5; %model timestep in seconds

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
count = 0;

%% Load and save model data!
for i = 1:length(modelScenarios)
    %Find model folders
    modelDataPath = [modelDir modelScenarios{i} '\Model\postProcessing\bathySample\surface\'];
    modelFolders = dir(modelDataPath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    modelFolders = {modelFolders.name};
    scenarioNumber = regexp(modelScenarios{i},'(\d+)','tokens');
    try
        modelID = find(strcmp(string(modelInfo.scenarioNumber),scenarioNumber{2}));
    catch
        continue % Skip folders not identified in the csv file
    end
    count = count+1;
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
    
    %Find hiZ from modelInfo
    %         if modelDepth == 0.5
    %             zMax = hiZ(1);
    %         elseif modelDepth == 1
    %             zMax = hiZ(2);
    %         elseif modelDepth == 2
    %             zMax = hiZ(3);
    %         elseif modelDepth == 3
    %             zMax = hiZ(4);
    %         elseif modelDepth == 4
    %             zMax = hiZ(5);
    %         end
    zMax = hiZ(count);
    
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
        if any((isempty(x) || isempty(y) || isempty(z) || isempty(Ux) || isempty(Uy) || isempty(Uz)))
            fprintf('Timestep %s contains null values!\n',modelFolders{sortID(j)})
            if j-1 ~= 0
                data.U.x(:,j) = NaN(length(data.U.x(:,j-1)),1);
                data.U.y(:,j) = NaN(length(data.U.y(:,j-1)),1);
                data.U.z(:,j) = NaN(length(data.U.z(:,j-1)),1);
                data.U.Ux(:,j) = NaN(length(data.U.Ux(:,j-1)),1);
                data.U.Uy(:,j) = NaN(length(data.U.Uy(:,j-1)),1);
                data.U.Uz(:,j) = NaN(length(data.U.Uz(:,j-1)),1);
            else
                continue
            end
        else
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
        end
        %Clear data to free memory
        clear x y z Ux Uy Uz file
        fclose(fid);
        
        %% LOAD p_rgh DATA
        isFile = contains(modelFiles,'p_rgh_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};p_rgh = file{4};
        if any((isempty(x) || isempty(y) || isempty(z) || isempty(p_rgh)))
            fprintf('Timestep %s contains null values!\n',modelFolders{sortID(j)})
            if j-1 ~= 0
                data.p_rgh.x(:,j) = NaN(length(data.p_rgh.x(:,j-1)),1);
                data.p_rgh.y(:,j) = NaN(length(data.p_rgh.y(:,j-1)),1);
                data.p_rgh.z(:,j) = NaN(length(data.p_rgh.z(:,j-1)),1);
                data.p_rgh.p_rgh(:,j) = NaN(length(data.p_rgh.p_rgh(:,j-1)),1);
            else
                continue
            end
        else
            %Allocate cropped data to structure
            data.p_rgh.x(:,j) = x(cropXZ);
            data.p_rgh.y(:,j) = y(cropXZ);
            data.p_rgh.z(:,j) = z(cropXZ);
            data.p_rgh.p_rgh(:,j) = p_rgh(cropXZ);
        end
        %Clear data to free memory
        clear x y z p_rgh file
        fclose(fid);
        
        %% LOAD k DATA
        isFile = contains(modelFiles,'k_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};k = file{4};
        if any((isempty(x) || isempty(y) || isempty(z) || isempty(k)))
            fprintf('Timestep %s contains null values!\n',modelFolders{sortID(j)})
            if j-1 ~= 0
                data.k.x(:,j) = NaN(length(data.k.x(:,j-1)),1);
                data.k.y(:,j) = NaN(length(data.k.y(:,j-1)),1);
                data.k.z(:,j) = NaN(length(data.k.z(:,j-1)),1);
                data.k.k(:,j) = NaN(length(data.k.k(:,j-1)),1);
            else
                continue
            end
        else
            %Allocate cropped data to structure
            data.k.x(:,j) = x(cropXZ);
            data.k.y(:,j) = y(cropXZ);
            data.k.z(:,j) = z(cropXZ);
            data.k.k(:,j) = k(cropXZ);
        end
        %Clear data to free memory
        clear x y z p_rgh file
        fclose(fid);
        
        %% LOAD epsilon DATA
        isFile = contains(modelFiles,'epsilon_sampledSurface');
        fid = fopen([modelDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
        file = textscan(fid,'%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};epsilon = file{4};
        if any((isempty(x) || isempty(y) || isempty(z) || isempty(epsilon)))
            fprintf('Timestep %s contains null values!\n',modelFolders{sortID(j)})
            if j-1 ~= 0
                data.epsilon.x(:,j) = NaN(length(data.epsilon.x(:,j-1)),1);
                data.epsilon.y(:,j) = NaN(length(data.epsilon.y(:,j-1)),1);
                data.epsilon.z(:,j) = NaN(length(data.epsilon.z(:,j-1)),1);
                data.epsilon.epsilon(:,j) = NaN(length(data.epsilon.epsilon(:,j-1)),1);
            else
                continue
            end
        else
            %Allocate cropped data to structure
            data.epsilon.x(:,j) = x(cropXZ);
            data.epsilon.y(:,j) = y(cropXZ);
            data.epsilon.z(:,j) = z(cropXZ);
            data.epsilon.epsilon(:,j) = epsilon(cropXZ);
        end
        %Clear data to free memory
        clear x y z epsilon file
        fclose(fid);
        
    end
    disp('Saving data to file')
    %Optionally save out data to disk
    if saveData == 1
        save([workingDir modelScenarios{i} '_rawData_' versionNO '.mat'],'-struct','data','-mat','-v7.3') %%%CHECK THIS!!!
    end
    clear data
    disp(['Model data processing completed in ' num2str(toc/60) ' minutes'])
end
