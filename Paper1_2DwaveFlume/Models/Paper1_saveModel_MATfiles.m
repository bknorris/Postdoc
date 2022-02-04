% Load and save Paper1_2DwaveFlume model data from the postProcessing
% folders in each model and save the data in MAT-format to disk. This
% script combines the "Paper1_saveModelFields_MATfiles" and
% "Paper1_saveModelFreeSurf_MATfiles" into one and incorporates folder
% unzipping. WARNING: Using this script with many model files will result
% in long processing times!
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
versionNO = 'V2';
saveData = 1; %set to 0 to turn off!

%% HR MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\HR_Domain\';
modelDir = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Scenarios\';
modelInfo = 'HR_Domain_model_toRun_082821.csv';

%HR MODEL Default settings:
waveGauges = [31 65 66 67 68 69 70 71 72 73 74 75 76];
loX = 65;hiX = 76;
loZ = -5;hiZ = [0 -0.25 -1 -1.5 -2]; %0.5, 1, 2, 3, 4 m models in that order
modelTimeStep = 0.125; %model timestep in seconds

%% LR MODEL SETTINGS - COMMENT OUT IF NOT USING %%%%%%%%%%%%%%%%%%%%%%%%%%%
% workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\LR_Domain\';
% modelDir = 'h:\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Scenarios\';
% modelInfo = 'LR_Domain_model_toRun_082121.csv';
% 
% %LR MODEL Default settings:
% waveGauges = [31 63 64 65 66 67 68 69 70 71 72 73 74];
% loX = 63;hiX = 74;
% loZ = -5;hiZ = [0 -0.25 -1 -1.5 -2]; %0.5, 1, 2, 3, 4 m models in that order
% modelTimeStep = 0.125; %model timestep in seconds

%% BEGIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
modelScenarios = dir(modelDir);
modelScenarios = modelScenarios(~ismember({modelScenarios.name},{'.','..'}));
modelScenarios = {modelScenarios.name};
count = 1;
for i = 1:length(modelScenarios)
    %Find model folders
    scenarioNumber = regexp(modelScenarios{i},'\_(.*).zip*','tokens');
    modelID = find(strcmp(string(modelInfo.scenarioNumber),scenarioNumber{:}));
    if isempty(modelID) %process only the folders listed in the CSV file
        continue
    else
        %Output to command window to track progress
        fprintf('Processing folder %0.0f of %0.0f: %s\n',count,length(modelInfo.scenarioNumber),modelScenarios{i})
        t1 = tic;
        %% Unzip file into Temp folder for post-processing
        %Check to see if file is already unzipped
        if isfolder([workingDir 'TempFiles\' modelScenarios{i}(1:end-4)])
            disp('Folder is unzipped in TempFiles/ directory already...')
        else
            t2 = tic;
            disp('Unzipping file...')
            %unzip([modelDir modelScenarios{i}],[workingDir 'TempFiles\'])
            
            %Use 7-zip to expedite unzipping process
            [status,result] = system(['"C:\Program Files\7-Zip\7z.exe" -y x ' '"' [modelDir modelScenarios{i}] '"' ' -o' '"' [workingDir 'TempFiles\'] '"']);
            disp(['File unzipped in ' num2str(toc(t2)/60) ' minutes'])
        end
        %Get folder name
        folderName = regexp(modelScenarios{i},'(.*).zip','tokens');
        %% Process Free Surface Data
        disp('Processing freeSurf data...')
        t3 = tic;
        rampTime = modelInfo.rampTime(modelID);
        
        %Sort modelFolders
        freeSurfDataPath = [workingDir 'TempFiles\' folderName{:}{1} '\MKKmodel\postProcessing\freeSurface\'];
        modelFolders = dir(freeSurfDataPath);
        modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
        modelFolders = {modelFolders.name};
        
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
            modelFiles = dir([freeSurfDataPath modelFolders{sortID(j)} '\']);
            modelFiles = {modelFiles.name};
            
            %% LOAD FREE SURFACE DATA
            isFile = contains(modelFiles,'alpha.water_freeSurface');
            fid = fopen([freeSurfDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
            file = textscan(fid,'%n%n%n%n','headerlines',2);
            x = file{1};y = file{2};z = file{3};water = file{4};
            %Loop through wave gauges and save data to structure
            for k = 1:length(waveGauges)
                fieldName = sprintf('WG%0.0f',k);
                %Check to see if there are null values in the model data file
                if any((isempty(x) || isempty(y) || isempty(z) || isempty(water)))
                    warning('Timestep %s contains null values!\n',modelFolders{sortID(j)})
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
            save([workingDir folderName{:}{1} '_freeSurf_' versionNO '.mat'],'-struct','data','-mat','-v7.3')
        end
        clear data
        disp(['freeSurf data processing completed in ' num2str(toc(t3)/60) ' minutes'])
        
        %% Process Field Data
        disp('Processing field data...')
        t4 = tic;
        rampTime = modelInfo.rampTime(modelID);
        modelDepth = modelInfo.h(modelID);
        
        %Sort modelFolders
        fieldDataPath = [workingDir 'TempFiles\' folderName{:}{1} '\MKKmodel\postProcessing\bathySample\surface\'];
        modelFolders = dir(fieldDataPath);
        modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
        modelFolders = {modelFolders.name};
        
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
        data.time = time;
        data.xLim = [loX hiX];
        data.zLim = [loZ zMax];
        
        %Load the data
        for j = 1:length(sortID)
            modelFiles = dir([fieldDataPath modelFolders{sortID(j)} '\']);
            modelFiles = {modelFiles.name};
            
            %% LOAD U DATA
            isFile = contains(modelFiles,'U_sampledSurface');
            fid = fopen([fieldDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
            file = textscan(fid,'%n%n%n%n%n%n','headerlines',2);
            x = file{1};y = file{2};z = file{3};Ux = file{4};Uy = file{5};Uz = file{6};
            if any((isempty(x) || isempty(y) || isempty(z) || isempty(Ux) || isempty(Uy) || isempty(Uz)))
                warning('Timestep %s contains null values!\n',modelFolders{sortID(j)})
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
                try
                    %Allocate cropped data to structure
                    data.U.x(:,j) = x(cropXZ);
                    data.U.y(:,j) = y(cropXZ);
                    data.U.z(:,j) = z(cropXZ);
                    data.U.Ux(:,j) = Ux(cropXZ);
                    data.U.Uy(:,j) = Uy(cropXZ);
                    data.U.Uz(:,j) = Uz(cropXZ);
                catch
                    warning('U Timestep %s is not the same length as the rest of the data!\n',modelFolders{sortID(j)})
                end
            end
            %Clear data to free memory
            clear x y z Ux Uy Uz file
            fclose(fid);
            
            %% LOAD p_rgh DATA
            isFile = contains(modelFiles,'p_rgh_sampledSurface');
            fid = fopen([fieldDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
            file = textscan(fid,'%n%n%n%n','headerlines',2);
            x = file{1};y = file{2};z = file{3};p_rgh = file{4};
            if any((isempty(x) || isempty(y) || isempty(z) || isempty(p_rgh)))
                warning('Timestep %s contains null values!\n',modelFolders{sortID(j)})
                if j-1 ~= 0
                    data.p_rgh.x(:,j) = NaN(length(data.p_rgh.x(:,j-1)),1);
                    data.p_rgh.y(:,j) = NaN(length(data.p_rgh.y(:,j-1)),1);
                    data.p_rgh.z(:,j) = NaN(length(data.p_rgh.z(:,j-1)),1);
                    data.p_rgh.p_rgh(:,j) = NaN(length(data.p_rgh.p_rgh(:,j-1)),1);
                else
                    continue
                end
            else
                try
                    %Allocate cropped data to structure
                    data.p_rgh.x(:,j) = x(cropXZ);
                    data.p_rgh.y(:,j) = y(cropXZ);
                    data.p_rgh.z(:,j) = z(cropXZ);
                    data.p_rgh.p_rgh(:,j) = p_rgh(cropXZ);
                catch
                    warning('p_rgh Timestep %s is not the same length as the rest of the data!\n',modelFolders{sortID(j)})
                end
            end
            %Clear data to free memory
            clear x y z p_rgh file
            fclose(fid);
            
            %% LOAD k DATA
            isFile = contains(modelFiles,'k_sampledSurface');
            fid = fopen([fieldDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
            file = textscan(fid,'%n%n%n%n','headerlines',2);
            x = file{1};y = file{2};z = file{3};k = file{4};
            if any((isempty(x) || isempty(y) || isempty(z) || isempty(k)))
                warning('Timestep %s contains null values!\n',modelFolders{sortID(j)})
                if j-1 ~= 0
                    data.k.x(:,j) = NaN(length(data.k.x(:,j-1)),1);
                    data.k.y(:,j) = NaN(length(data.k.y(:,j-1)),1);
                    data.k.z(:,j) = NaN(length(data.k.z(:,j-1)),1);
                    data.k.k(:,j) = NaN(length(data.k.k(:,j-1)),1);
                else
                    continue
                end
            else
                try
                    %Allocate cropped data to structure
                    data.k.x(:,j) = x(cropXZ);
                    data.k.y(:,j) = y(cropXZ);
                    data.k.z(:,j) = z(cropXZ);
                    data.k.k(:,j) = k(cropXZ);
                catch
                    warning('k Timestep %s is not the same length as the rest of the data!\n',modelFolders{sortID(j)})
                end
            end
            %Clear data to free memory
            clear x y z p_rgh file
            fclose(fid);
            
            %% LOAD epsilon DATA
            isFile = contains(modelFiles,'epsilon_sampledSurface');
            fid = fopen([fieldDataPath modelFolders{sortID(j)} '\' modelFiles{isFile}]);
            file = textscan(fid,'%n%n%n%n','headerlines',2);
            x = file{1};y = file{2};z = file{3};epsilon = file{4};
            if any((isempty(x) || isempty(y) || isempty(z) || isempty(epsilon)))
                warning('Timestep %s contains null values!\n',modelFolders{sortID(j)})
                if j-1 ~= 0
                    data.epsilon.x(:,j) = NaN(length(data.epsilon.x(:,j-1)),1);
                    data.epsilon.y(:,j) = NaN(length(data.epsilon.y(:,j-1)),1);
                    data.epsilon.z(:,j) = NaN(length(data.epsilon.z(:,j-1)),1);
                    data.epsilon.epsilon(:,j) = NaN(length(data.epsilon.epsilon(:,j-1)),1);
                else
                    continue
                end
            else
                try
                    %Allocate cropped data to structure
                    data.epsilon.x(:,j) = x(cropXZ);
                    data.epsilon.y(:,j) = y(cropXZ);
                    data.epsilon.z(:,j) = z(cropXZ);
                    data.epsilon.epsilon(:,j) = epsilon(cropXZ);
                catch
                    warning('Epsilon Timestep %s is not the same length as the rest of the data!\n',modelFolders{sortID(j)})
                end
            end
            %Clear data to free memory
            clear x y z epsilon file
            fclose(fid);
        end
        disp('Saving data to file')
        %Optionally save out data to disk
        if saveData == 1
            save([workingDir folderName{:}{1} '_rawData_' versionNO '.mat'],'-struct','data','-mat','-v7.3') %%%CHECK THIS!!!
        end
        clear data
        disp(['Field data processing completed in ' num2str(toc(t4)/60) ' minutes'])
        count = count+1;
    end
    disp('Data processing complete, removing tempFile...')
    fclose('all');
    rmdir([workingDir 'TempFiles\' folderName{:}{1}],'s')
    disp(['Model data processing completed in ' num2str(toc(t1)/60) ' minutes'])
end