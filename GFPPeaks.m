%% GFP PEAK IDENTIFICATION
% This script automatically identifies the group-averaged GFP peak inside a
% period of interest (POI)

% Corentin Wicht
% Laboratory for Neurorehabilitation Science, UNIFR, CH
% 02.06.2020

% TO DO:
% - Do not let possibility to load different "Chanlocs", the files should
% always be avg. referenced!! Thus, need to load an avg. ref. mat file of
% EEGLAB chanlocs

% Method:
% 1. identifier la latence des individual GFP peak (une valeur en ms par
% sujet) visuellement (based on latency + topography of expected component) !!!
% 2. extraire le SD de ce groupe de valeur
% 3. extraire le peak de la moyenne des GFP ou du GFP de la moyenne (à rediscuter)
% 4. définir la POI comme GFPmean (3.) +-SD (2.).

clear variables
%% SETUP

% getting path of the script location
p = matlab.desktop.editor.getActiveFilename;
I_p = strfind(p,'\');
p2 = p(1:I_p(end)-1);

% Path of all needed functions
addpath(genpath([p2 '\Functions']));

% Format of files to load and method to apply
PromptIndivGFP = inputdlg({['Extension of the files to be loaded:' newline ...
    '1) ep' newline '2) eph'],['Method to use to identify individual GFP peaks' newline ...
    '1) semi-automatic' newline '2) manual' newline '3) Load previous data']},'Settings',1,{'1','1'});
PromptIndivGFP = str2double(PromptIndivGFP);

% Extension
if PromptIndivGFP(1)==1; Extension='.ep'; elseif PromptIndivGFP(1)==2; Extension='.eph'; end

% Method
if PromptIndivGFP(2)==1; Method='semi-automatic'; elseif PromptIndivGFP(2)==2; Method='manual'; 
elseif PromptIndivGFP(2)==3; Method='load'; end

% Path of your upper folder containing your .EP data
root_folder = uigetdir(p2,'Choose the path of your most upper folder containing your files.');
cd(root_folder)
FileList = dir(['**/*' Extension]);
% OutputPath = strsplit(root_folder,'\');
% OutputPath = strcat(OutputPath(1:end-1),'\');
% OutputPath = strcat(OutputPath{:});
OutputPath = root_folder;

% Sort FileList according to natural order
[~, Idx] = natsort({FileList.name});
FileList = FileList(Idx);

% Define design
% Default values for the GUI below 
to_display = [{'B';'Group';'OH_OH';'PBO_OH';''} cell(5,3)];

% Table
ScreenSize=get(0,'ScreenSize');
f = figure('Position', [ScreenSize(3)/2-500/2 ScreenSize(4)/2-500/2 500 600]);
p = uitable('Parent', f,'Data',to_display,'ColumnEdit',true(1,size(to_display,2)),'ColumnName',...
    {'Factor1','Factor2','IgnoreFolders','IgnoreFiles'},'CellEditCallBack','DesignList = get(gco,''Data'');');
p.Position = [50 -100 350 400];
uicontrol('Style', 'text', 'Position', [20 350 450 220], 'String',...
        {['DESIGN DEFINITION' newline ''],['The definition of the design needs to be structured in the following way:'...
        newline 'FOR THE FIRST TWO COLUMNS (i.e. analysis factors)'...
        newline '1st line = Within (W) or Between-subject (B)'...
        newline '2nd line = Name of the factor'...
        newline '3rd+ lines = Name of the levels'...
        newline ...
        newline 'FOR THE LAST TWO COLUMNS (i.e. data removal)'...
        newline 'The last two columns enables to remove undesired data (e.g. at the level of i) folders or ii) file names'...
        newline '' newline '! THE NAMES MUST BE THE SAME AS YOUR FILES/FOLDERS !']});
    
% Wait for t to close until running the rest of the script
waitfor(p)

% If no modifications to the example in the figure
if ~exist('DesignList','var')
    DesignList = to_display;
end

% Removing files/folders if specified
if nnz(~cellfun('isempty',(DesignList(:,3)))) % folders
    for k=1:sum(~cellfun('isempty',(DesignList(:,3))))
        FileList = FileList(~contains({FileList.folder},DesignList{k,3}));
    end
elseif nnz(~cellfun('isempty',(DesignList(:,4)))) % files
    for k=1:sum(~cellfun('isempty',(DesignList(:,4))))
        FileList = FileList(~contains({FileList.name},DesignList{k,4}));
    end
end

% Reshaping the design structure
Design=[];
Factors = [0 0];
for k=1:2 % never more than 2 factors
    if strcmpi(DesignList{1,k},'W')
        Factors(1) = Factors(1) + 1;
        Design.Within.(sprintf('Factor%d',Factors(1))){k}=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    elseif strcmpi(DesignList{1,k},'B')
        Factors(2) = Factors(2) + 1;
        Design.Between.(sprintf('Factor%d',Factors(2))){k}=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    end
    
    % Making sure that each levels are unique !!!
    %i.e., not self-contained such as "Pain" vs "NoPain"
    % It would mess with the "contains" function (line 244 XXX)
    TempLevels = DesignList(3:end,k);
    TempLevels = TempLevels(~cellfun('isempty',TempLevels));
    AllComb = nchoosek(TempLevels,2); % Number of unique combinations
    SelfCont = [false false]; % comparing them 2-by-2
    for m=1:size(AllComb,1)
        SelfCont(1) = contains(AllComb{m,1},AllComb{m,2}); % 2-in-1
        SelfCont(2) = contains(AllComb{m,2},AllComb{m,1}); % 1-in-2
        
        % In case of self-containing levels
        if any(SelfCont) 
           Direction = fastif(SelfCont(1)==true, [AllComb{m,2} ' in ' AllComb{m,1}], ...
           [AllComb{m,1} ' in ' AllComb{m,2}]);
           error(['WARNING: the following levels are self-contained for factor "' DesignList{2,k} '":' ...
               newline Direction, '.', newline ...
               'You should make them UNIQUE identifiers otherwise the program will fail to load the files separately for each level.']) 
        end
    end
end

% Rename conditions if double-blind
PromptBlind = questdlg('Would you like to rename the levels (e.g. in case of double-blind design) ?', ...
	'Unblinding','Yes','No','Yes');

% Type of unblinding
if strcmpi(PromptBlind,'Yes')
    PromptBlind = questdlg('Should the levels be renamed for all files or individually ?', ...
        'Unblinding type','All','Individual','Individual');
end

% Prompt blinding
% This is the case where for files contain the name of the blinded condition 
if strcmpi(PromptBlind,'All')
    DesignFields = fieldnames(Design);
    for k=1:numel(DesignFields) % for each within/between

        FactorFields = fieldnames(Design.(DesignFields{k}));
        for j=1:numel(FactorFields) % for each factor
            
            % Isolate the levels
            LevelsList = Design.(DesignFields{k}).(FactorFields{j}){:}';
            
            % Run prompt
            PromptInstructBlind = [['Enter the new names of each level for' newline ...
                [LevelsList{1} '-Subject factor "' LevelsList{2} '":'] newline ...
                LevelsList{3}] LevelsList(4:end)];
            PromptValBlind = repmat({''},1,length(LevelsList)-2);
            PromptInpBlind = inputdlg(PromptInstructBlind,'Unblinding',1,PromptValBlind);
            
            % Saving new levels' name along the old ones
            Design.(DesignFields{k}).(FactorFields{j}){:} = horzcat ...
                (Design.(DesignFields{k}).(FactorFields{j}){:},[LevelsList(1:2)';PromptInpBlind]);
        end
    end
    
% This is the case where for files contain the name of the session but 
% the sessions order was randomized
elseif strcmpi(PromptBlind,'Individual')
    
    % Save a temporary excel file (first column is filenames)
    Data = horzcat(['FILENAMES';{FileList.name}'],['NEW LEVEL NAMES';cell(length({FileList.name}),1)]);
    xlswrite('IndividualUnblinding.xlsx',Data)
    
    % Open excel
    e=actxserver('excel.application');
    eW=e.Workbooks;
    eF=eW.Open([pwd '\IndividualUnblinding.xlsx']); % open OutputTest.xls
    eS=eF.ActiveSheet;
    e.visible = 1; % If you want Excel visible.
    
%     % Message box stoping code execution
    MessageBox('The code will continue once you press OK','Wait for user input',30,250,70)
    
    try % to avoid errors if the excel file was saved/closed before the message box
        % edit sheet
        eF.Save;
        eF.Close; % close the file
        e.Quit; % close Excel entirely
    catch
    end
    
    % Load the data and then delete the excel file
    Unblinding = readcell([pwd '\IndividualUnblinding.xlsx']);
    Unblinding = Unblinding(2:end,:);
%     delete([pwd '\IndividualUnblinding.xlsx']);
end

% Analysis on all data or specific files
PromptAnalyses = questdlg('Would you like to perform the analysis on all your data or on specific files ?', ...
	'Selection of analysis','All data','Specific files','All data');

% Retrieve names from the FileList structure
AllNames=natsort(unique({FileList.name}));

% If user decides to restrict analysis to specified folders
if strcmp(PromptAnalyses,'Specific files')

    % Matrix to integrate in the following uitable
    to_display = [AllNames', repmat({false},[size(AllNames,2) 1])];

    % Select folders on which to apply analyses
    f = figure('Position', [ScreenSize(3)/2-500/2 ScreenSize(4)/2-500/2 400 400]);
    p=uitable('Parent', f,'Data',to_display,'ColumnEdit',[false true],'ColumnName',...
        {'Files', 'Perform analysis ?'},'CellEditCallBack','SbjList = get(gco,''Data'');',...
        'ColumnWidth', {140,120});
    uicontrol('Style', 'text', 'Position', [20 325 200 50], 'String',...
            {'Folder selection for GFP peak analysis','Click on the box of the participants files you want to perform analyses on'});
    % Wait for t to close until running the rest of the script
    waitfor(p)

    % Stores the files on which to apply IC decomposition
    ToAnalyse=find(cell2mat(SbjList(:,2)));

    % Recreates a new FileList structure based on selected folders
    FileList = FileList(ToAnalyse);
    AllNames = AllNames(ToAnalyse);
end

% Load each .ep file and store it in structure according to design
EEGData = []; Names = []; Fields = fieldnames(Design); 
if strcmpi(PromptBlind,'No'); Unblinding = AllNames'; end
for k=1:length(Fields) % For BS and/or WS category
    Factors = fieldnames(Design.(Fields{k}));
    for j=1:length(Factors) % For each factor
        
        % NO UNBLINDING
        if strcmpi(PromptBlind,'No')
            Levels = Design.(Fields{k}).(Factors{j}){:}(3:end);
            for m=1:length(Levels) % For each level of the factor
                LevelsFilesList = FileList(contains({FileList.name},Levels{m}));
                for n=1:length(LevelsFilesList) % For each file
                    % Load .ep(h) files
                    % Taken from "Matlabroutinesforcartoolers" functions 
                    % https://sites.google.com/site/cartoolcommunity/files
                    Path = [LevelsFilesList(n).folder '\' LevelsFilesList(n).name];
                    [numchannels,numtimeframes,samplingrate,thedata]=openeph(Path);
                    EEGData.(Levels{m})(:,:,n) = thedata; % Store in main database
                    Names.(Levels{m})(n,1) = {LevelsFilesList(n).name};
                    disp(['File ' LevelsFilesList(n).name ' successfully loaded'])
                    
                    % For excel output
                    Pos = find(contains(Unblinding(:,1),LevelsFilesList(n).name));
                    Unblinding{Pos,j+1} = Levels{m}; 
                end
            end
            
        % UNBLINDING    
        else
            if strcmpi(PromptBlind,'All')
                OldLevels = Design.(Fields{k}).(Factors{j}){:}(3:end,1);
                NewLevels = Design.(Fields{k}).(Factors{j}){:}(3:end,2);
                for m=1:length(OldLevels) % For each level of the factor
                    LevelsFilesList = FileList(contains({FileList.name},OldLevels{m}));
                    for n=1:length(LevelsFilesList) % For each file
                        % Load .ep(h) files
                        Path = [LevelsFilesList(n).folder '\' LevelsFilesList(n).name];
                        [numchannels,numtimeframes,samplingrate,thedata]=openeph(Path);
                        EEGData.(NewLevels{m})(:,:,n) = thedata; % Store in main database (new name!)
                        Names.(NewLevels{m})(n,1) = {LevelsFilesList(n).name};
                        disp(['File ' LevelsFilesList(n).name ' successfully loaded'])
                    end
                end
            elseif strcmpi(PromptBlind,'Individual')
                UniqueUnblinding = unique(Unblinding(:,2));
                Pos = ones(1,length(UniqueUnblinding));
                for n=1:size(Unblinding,1) % For each file
                    % Load .ep(h) files
                    Path = [FileList(n).folder '\' FileList(n).name];
                    [numchannels,numtimeframes,samplingrate,thedata]=openeph(Path);
                    WhichLev = find(contains(UniqueUnblinding,Unblinding(n,2))==1);
                    EEGData.(Unblinding{n,2})(:,:,Pos(WhichLev)) = thedata; % Store in main database (new name!)
                    Names.(Unblinding{n,2})(Pos(WhichLev),1) = Unblinding(n,1);
                    Pos(WhichLev) = Pos(WhichLev) + 1;
                    disp(['File ' FileList(n).name ' successfully loaded'])
                end
            end
        end
    end
end

% Prompt for EEG parameters
PromptInstructions = {'Enter the epoching interval (in ms)','Enter the sampling rate:','Enter the output folder name'};
PromptValues = {'-100 700','1024','Output'};
PromptInputs = inputdlg(PromptInstructions,'EEG parameters',1,PromptValues);
Epoch = str2num(PromptInputs{1});
SamplingRate = str2double(PromptInputs{2});
OutputFolder = PromptInputs{3}; mkdir([OutputPath '\' OutputFolder]);
if numchannels>63
    Chanlocs = load([p2 '\Functions\ChanLocs_CMSDRL.mat']); Chanlocs = Chanlocs.TEMP;
else
    Chanlocs = load([p2 '\Functions\ChanLocs_Cz.mat']); Chanlocs = Chanlocs.TEMP;
end

% Define the EEG components of interest
% Default values for the GUI below 
to_display = cell(4,5);

% Define Components and bounds
ScreenSize=get(0,'ScreenSize');
f = figure('Position', [ScreenSize(3)/2-500/2 ScreenSize(4)/2-500/2 600 600]);
p = uitable('Parent', f,'Data',to_display,'ColumnEdit',true(1,size(to_display,2)),'ColumnName',...
    {'Name','Lower bound (ms)','Upper bound (ms)','Electrode Visual', 'Cluster of electrodes'},'Position',[50 -100 500 400],...
    'CellEditCallBack','CompList = get(gco,''Data'');');
uicontrol('Style', 'text', 'Position', [60 350 450 130], 'String',...
        {['COMPONENTS DEFINITION' newline ''],['The definition of the Period(s) of Interest needs to be structured in the following way:'...
        newline '1st column = Name of the Component (e.g. N2)',...
        newline '2nd column = Corresponding lower bound (e.g. 200ms post-stim.)',...
        newline '3rd column = Corresponding higher bound (e.g. 350ms post-stim.)',...
        newline '4th column = Index for electrode(s) of interest for each component (visualization purpose)',...
        newline '5th column = Cluster of electrodes to compute average voltage amplitude (leave empty if all electrodes)']});
% Wait for t to close until running the rest of the script
waitfor(p)

% If no modifications to the example in the figure
if ~exist('DesignList','var')
    CompList = to_display;
end

% Plotting parameters
% Ticks
XTStr = Epoch(1):50:Epoch(2);
TickSpacing = 50*(SamplingRate/1000); % Needs to be adjusted based on sampling rate
Ticks = 1:TickSpacing:Epoch(2)+abs(Epoch(1))+TickSpacing;

%% COMPUTE GFP
% Taken from "Matlabroutinesforcartoolers" functions 
% https://sites.google.com/site/cartoolcommunity/files
Fields = fieldnames(EEGData); 

% Create a .txt file for outputs
date_name = datestr(now,'dd-mm-yy_HHMM');
fid = fopen([OutputPath '\' OutputFolder '\Results_' date_name '.txt'],'w');
fprintf(fid,'%s\r\n','---- EEG PARAMETERS SUMMARY ----');
fprintf(fid,'%s\r\n',['Epoch (in ms) = ' num2str(Epoch)]);
fprintf(fid,'%s\r\n',['Sampling rate = ' num2str(SamplingRate)]);
fprintf(fid,'%s\r\n',['Number of time frames (TF) = ' num2str(numtimeframes)]);
fprintf(fid,'%s\r\n',['Number of channels = ' num2str(numchannels)]);
fprintf(fid,'%s\r\n',['Peak identification method used : ' Method]);

% Creating the .txt output for averaged GFP
% fid2 = fopen([OutputPath '\' OutputFolder '\AvgGFP_' date_name '.txt'],'w');
        
for n=1:sum(~cellfun(@(x) isempty(x), CompList(:,1))) % For each Comp
    
    % Initialize variables
    Pos = []; Numbers = repmat({[]},[size(AllNames,2) 1]);
    
    % Remove empty lines and store data for the current component
%     EmptIdx = ~cellfun(@(x) isempty(x),CompList(n,:));
%     CompN = CompList(n,EmptIdx); 
    CompN = CompList(n,:);
    
    % Component
    TEMPComp = cellfun(@(x) str2double(x),CompN(2:3)); % Load current Component bounds
    if any(Epoch<0)
        TEMPComp = TEMPComp + abs(Epoch(1)); % Move based on prestim
    end
    TEMPCompinTF = round(TEMPComp/1000*SamplingRate); % Convert to TimeFrames and round to nearest integer
    
    % Electrodes cluster
    if ~isempty(CompN{end})
        ElectClust = str2num(CompN{end});
    else
        TEMPFields = fieldnames(EEGData);
        ElectClust = 1:size(EEGData.(TEMPFields{1}),2); % all electrodes
    end
    
    for k=1:length(Fields) % For each levels
        for m=1:size(EEGData.(Fields{k}),3) % For each subject
            EEGTEMP = squeeze(EEGData.(Fields{k})(:,:,m));
            GFPData.(Fields{k})(:,m)=computegfp(EEGTEMP);

            %% PEAK IDENTIFICATION ON INDIVIDUAL GFP
            if strcmpi(Method,'semi-automatic')
                % Look for maximum/minimum inside Component bounds
                TEMPGFP = GFPData.(Fields{k})(:,m);
                MaxValGFP = max(TEMPGFP(TEMPCompinTF(1):TEMPCompinTF(2)));
                MaxPos = find(TEMPGFP==MaxValGFP);
                
                Fig = figure('units','normalized','outerposition',[0 0 1 1]); % Initialize figure in full screen
                File = strsplit(strrep(Names.(Fields{k}){m},'_','-'),'.'); 
                Title = [CompN{1} ': ' Fields{k} ' ' File{1}]; % Main title (component: level_subject)
                  annotation('textbox',[.4 .5 .5 .5],'String',Title,'FitBoxToText',...
                  'on','FontSize',20,'EdgeColor','w');
                
                % 1. Electrodes level data
                % Color vector
                lineWid = ones(size(EEGTEMP,2),1); % line width
                if length(CompN)>3
%                     ElectToCol = str2double(CompList(~cellfun(@(x) isempty(x),CompList(:,end)),end));
                    str = CompN(~cellfun(@(x) isempty(x),CompN(:,end-1)),end-1);
                    ElectToCol = regexp(str,'(-)?\d+(\.\d+)?(e(-|+)\d+)?','match');
                    ElectToCol = str2double(ElectToCol{:});
                    ElseCol = lines(length(ElectToCol));
                    ColVect = repmat(128/255,size(EEGTEMP,2),3);
                    ColVect(ElectToCol,:) = ElseCol;
                    lineWid(ElectToCol) = repmat(6,length(ElectToCol),1);
                else
                    ColVect = zeros(size(EEGTEMP,2),3);
                end

                subplot(2,3,[1 2]);
                for t = 1:size(EEGTEMP,2)
                    Handles(t) = plot(EEGTEMP(:,t),...
                        'Color',ColVect(t,:),'LineWidth',lineWid(t)); 
                    hold on
                end
                hold off
                title('All electrodes'); % axis tight;
                xlabel('Time (ms)'), ylabel('uV');
                set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr)); 
                axis tight; Yl=ylim(gca); % retrieve auto y-limits   

                % Vertical line
                PosinMS = round(MaxPos/SamplingRate*1000);
                LabelValue = [num2str(PosinMS + Epoch(1)) 'ms ' '[' num2str(MaxPos) 'TF' ']']; % num2str(MaxVal)
                xline(MaxPos,'-',LabelValue,'Color','b','LabelVerticalAlignment','bottom');

                PosRectinTF = [TEMPCompinTF(1) Yl(1) TEMPCompinTF(2)-TEMPCompinTF(1) Yl(2)*0.1]; 
                rectangle('Position',PosRectinTF,'FaceColor',[0 0 1 0.2],'EdgeColor',[0 0 1 0.2]);
                             
                % Extract the handles that require legend entries
                HandlesToLeg = Handles(ElectToCol);
                
                % Legend
                LegendLabel1 = CompList(~cellfun(@(x) isempty(x),CompList(:,end)),end);
                Lgd = legend(HandlesToLeg,LegendLabel1); 
                title(Lgd,'Electrode','Color','k');
            
                % 2. GFP
                subplot(2,3,[4 5]); plot(TEMPGFP, 'Color','k'); 
                title('GFP'); % axis tight;
                title(sprintf('GFP (%f\\muV)',TEMPGFP(MaxPos)));
                axis tight; Yl=ylim(gca); % retrieve auto y-limits
                xlabel('Time (ms)'), ylabel('GFP');
                set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr)); 

                % Vertical line
                xline(MaxPos,'-',LabelValue,'Color','b','LabelVerticalAlignment','bottom');

                % Patches
                PosRectinTF = [TEMPCompinTF(1) Yl(1) TEMPCompinTF(2)-TEMPCompinTF(1) Yl(2)*0.1]; 
                rectangle('Position',PosRectinTF,'FaceColor',[0 0 1 0.2],'EdgeColor',[0 0 1 0.2]);

                % Legend
                LegendLabel2 = {'GFP';CompN{1}};
                Lgd = legend(LegendLabel2); title(Lgd,'Component','Color','k');
                
                % Store parameters
                Param.Ticks=Ticks;
                Param.XTStr=XTStr;
                Param.TEMPCompinTF=TEMPCompinTF;
                Param.Epoch=Epoch;
                Param.SamplingRate=SamplingRate;
                Param.CompN=CompN{1};
                Param.ColVect = ColVect;
                Param.lineWid = lineWid;
                Param.LegendLabel = LegendLabel1;
                Param.ElectToCol = ElectToCol;
                
                % 3. TOPOPLOTS
                slider_plot(Fig,EEGTEMP,MaxPos,Chanlocs,Param,TEMPGFP)
                
                % RESPONSE TABLE
                % Finding position of current file in table
                Pos = find(ismember(AllNames',Names.(Fields{k}){m})==1);
                
                % Matrix to integrate in the following uitable
                Numbers{Pos} = MaxPos;
                to_display = [AllNames',Numbers];

                % Table enabling to adjust the preselected value
                f = figure('Position', [125 125 400 400]);
                p=uitable('Parent', f,'Data',to_display,'ColumnEdit',[false true],'ColumnName',...
                    {'Files', 'Peak GFP (in TF)'},'ColumnWidth', {180,100},'Position',[10 20 350 300],...
                'CellEditCallBack','MaxGFPList = get(gco,''Data'');');
                uicontrol('Style', 'text', 'Position', [80 325 200 50], 'String',...
                        {'Individual GFP peaks','Close the box once you are settled with current file.'});
                
                % If no modifications to the table
                if ~exist('MaxGFPList','var')
                    MaxGFPList = to_display;
                end
                
                % Wait for p to close until running the rest of the script
                waitfor(p);  close gcf;
                
                % Store results
                MaxIndivGFP.(CompN{1}).(Fields{k})(m) = MaxGFPList{Pos,2}; 
                to_display{m,2} = MaxGFPList{Pos,2}; Numbers{Pos} = MaxGFPList{Pos,2};
                clear Fig MaxGFPList;
            end
        end
    end
    
    % Compute mean GFP for each level
    MeanGFP = structfun(@(x) mean(x,2),GFPData,'UniformOutput',0);
    if strcmpi(Method,'semi-automatic')
        TEMP = cell2mat(struct2cell(MaxIndivGFP.(CompN{1})));
        MaxGFPList = TEMP(:);
        MeanPeakGFP = structfun(@(x) mean(x,2),MaxIndivGFP.(CompN{1}),'UniformOutput',0);

        % Compute SD over individual GFP peaks for each level
        SDPeakGFP = structfun(@(x) std(x),MaxIndivGFP.(CompN{1}),'UniformOutput',0);
    end

    % MANUAL INDIVIDUAL GFP PEAKS LOCALISATION (e.g. using CARTOOL)
    if strcmpi(Method,'Manual')

        % Matrix to integrate in the following uitable
        to_display = [AllNames',Numbers];
        
        
        % If excel is installed (easier to copy-paste)
        % see : https://ch.mathworks.com/matlabcentral/answers/631794-how-to-check-if-excel-is-installed
        try
            Excel = matlab.io.internal.getExcelInstance; % Would issue an ActiveX error if excel not installed
            
            % Save a temporary excel file (first column is filenames)
            xlswrite('TempFile.xlsx',to_display,CompN{1});

            % Open excel
            e=actxserver('excel.application');
            eW=e.Workbooks;
            eF=eW.Open([pwd '\TempFile.xlsx']); 
            eS=eF.ActiveSheet;
            e.visible = 1; % If you want Excel visible.

            % Message box stoping code execution
            MessageBox('The code will continue once you press OK','Wait for user input',30,250,70)

            try % to avoid errors if the excel file was saved/closed before the message box
                % edit sheet
                eF.Save;
                eF.Close; % close the file
                e.Quit; % close Excel entirely
            catch
            end
    
            % Load the data and then delete the excel file
            MaxGFPList = readcell([pwd '\TempFile.xlsx'],'sheet',CompN{1});
            MaxGFPList = MaxGFPList(:,1:2); % Only keeping first 2 columns
            delete([pwd '\TempFile.xlsx']);
            
            % Reminder
            clc; % Clear command window
            fprintf('<strong> Current component of interest: %s [%d - %d TF] </strong> \n',...
                CompN{1},TEMPCompinTF(1),TEMPCompinTF(2))
            
            % Extract numbers
            MaxGFPList = cell2mat(MaxGFPList(:,2));
        
        % If excel is not installed
        catch
            
            % Select folders on which to apply analyses
            f = figure('Position', [125 125 400 400]);
            p=uitable('Parent', f,'Data',to_display,'ColumnEdit',[false true],'ColumnName',...
                {'Files', 'Peak GFP'},'ColumnWidth', {180,60},...
            'CellEditCallBack','MaxGFPList = get(gco,''Data'');');
            uicontrol('Style', 'text', 'Position', [20 325 250 50], 'String',...
                    {sprintf('Individual GFP peaks (in TF !) for component : %s',CompN{1}),...
                    'Close the box once you are settled with all files.'});
                
            % Reminder
            clc; % Clear command window
            fprintf('<strong> Current component of interest: %s [%d - %d TF] </strong> \n',...
                CompN{1},TEMPCompinTF(1),TEMPCompinTF(2))

            % Wait for p to close until running the rest of the script
            waitfor(p); 
            
            % Change data type
            MaxGFPList = cellfun(@(x) str2double(x),MaxGFPList(:,2));
            
        end
        
        % Identify numbers out of component's bounds
        OutOfBounds = AllNames(MaxGFPList<TEMPCompinTF(1) | MaxGFPList>TEMPCompinTF(2))';
        
        % Warning for files that are out of component's bounds  
        if ~isempty(OutOfBounds)
            warning('The following files are out of component %s bound:',CompN{1});
            fprintf('%s \n',OutOfBounds{:});
        end        
        
        % Build structure by levels
        for k=1:length(Fields) 
            MaxGFPStructTF.(Fields{k}) = MaxGFPList(contains(Unblinding(:,2),Fields{k}));
        end
        
        % Compute mean GFP for each level
        MeanPeakGFP = structfun(@(x) mean(x),MaxGFPStructTF,'UniformOutput',0);

        % Compute SD over individual GFP for each level
        SDPeakGFP = structfun(@(x) std(x),MaxGFPStructTF,'UniformOutput',0);
        
    % Loading previous analyses
    elseif strcmpi(Method,'Load')
        
        [File, Path] = uigetfile('*.mat',['Load previous data for component: ' ...
            CompN{1}]);
        VarLoad = {'MeanPeakGFP','MeanGFP','SDPeakGFP','MaxGFPList'};
        load([Path File],VarLoad{:})
    end

    %% PEAK IDENTIFICATION ON MEAN GFP
    Fig = figure('units','normalized','outerposition',[0 0 1 1]); % Initialize figure in full screen
    annotation('textbox',[.45 .5 .5 .5],'String',...
        sprintf('Results for component %s',CompN{1}),'FitBoxToText',...
        'on','FontSize',20,'EdgeColor','w');
    
    % Initialize empty fields
    MaxGFP.MeanRange = [];
    MaxGFP.SDRange = [];
    MaxVolt.MeanRange = [];
    MaxVolt.SDRange = [];
    
    % Min & Max GFP value across levels
    MinValGFP = min(cell2mat(struct2cell(MeanGFP)));
    MaxValGFP = max(cell2mat(struct2cell(MeanGFP)));
    
    % Min & Max GFP value across levels
    MinValVolt = min(cell2mat(struct2cell(...
        structfun(@(x) mean(x(:,ElectClust,:),[2,3]), EEGData, 'UniformOutput', 0))));
    MaxValVolt = max(cell2mat(struct2cell(...
        structfun(@(x) mean(x(:,ElectClust,:),[2,3]), EEGData, 'UniformOutput', 0))));
    
    % For .txt output
    fprintf(fid,'\r\n%s\r\n',['----' CompN{1} ' COMPONENT [' ...
        num2str(TEMPCompinTF(1)) '-'  num2str(TEMPCompinTF(2)) ' TF]'  '----']);
%     fprintf(fid2,'\r\n%s\r\n',['----' CompN{1} ' COMPONENT [' ...
%         num2str(TEMPCompinTF(1)) '-'  num2str(TEMPCompinTF(2)) ' TF]'  '----']);
    
    for m=1:length(Fields) % For each level
        
        % Initialize data
        MaxGFP.(Fields{m}).Pos = round(MeanPeakGFP.(Fields{m}));
        MaxGFP.(Fields{m}).Value = MeanGFP.(Fields{m})(MaxGFP.(Fields{m}).Pos);
        EEGTEMP = squeeze(mean(EEGData.(Fields{m})(:,:,:),3));

        %% PLOTTING
        
        % 1) VOLTAGE AMPLITUDE PLOT AND PARAMETERS
        subplot(3,length(Fields),m)
        PlotDat = squeeze(mean(EEGData.(Fields{m})(:,ElectClust,:),[2,3]));
        plot(PlotDat,'Color','k');axis tight;
        xlabel('Time (ms)'), ylabel('Mean Voltage Amplitude [\muV]')
        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr),'ylim',[MinValVolt MaxValVolt]);
        
        % Vertical line (Peak)
        PosinTF = MaxGFP.(Fields{m}).Pos;
        PosinMS = round(PosinTF/SamplingRate*1000);
        LabelValue = [num2str(PosinMS + Epoch(1)) 'ms ' '[' num2str(PlotDat(PosinTF)) ']'];
        xline(MaxGFP.(Fields{m}).Pos,'-',LabelValue,'Color',[1 0 0],'LabelVerticalAlignment','bottom');
        
        % Patches (+- 1SD around the peak)
        TEMP = structfun(@(x) squeeze(mean(x(:,ElectClust,:),2)),EEGData,'UniformOutput',0);
        MaxTEMP = cell2mat(struct2cell(structfun(@(x) max(x,[],[1,2]),TEMP,'UniformOutput',0)));
        MinTEMP = cell2mat(struct2cell(structfun(@(x) min(x,[],[1,2]),TEMP,'UniformOutput',0)));
        PosRectinTF = [PosinTF-SDPeakGFP.(Fields{m}) min(MinTEMP) ...
            2*SDPeakGFP.(Fields{m}) max(MaxTEMP)+abs(min(MinTEMP))];
        rectangle('Position',PosRectinTF,'FaceColor',[1 0 0 0.2],'EdgeColor',[1 0 0 0.2]);
        
        % 2) GFP PLOT AND PARAMETERS
        subplot(3,length(Fields),m+length(Fields))
        plot(MeanGFP.(Fields{m}),'Color','k');axis tight;
        xlabel('Time (ms)'), ylabel('Mean GFP')
        title(sprintf('%s',strrep(Fields{m},'_','-'))) 
        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr),'ylim',[MinValGFP MaxValGFP]);
        
        % Vertical line (Peak)
        PosinTF = MaxGFP.(Fields{m}).Pos;
        PosinMS = round(PosinTF/SamplingRate*1000);
        LabelValue = [num2str(PosinMS + Epoch(1)) 'ms ' '[' num2str(MaxGFP.(Fields{m}).Value) ']'];
        xline(MaxGFP.(Fields{m}).Pos,'-',LabelValue,'Color',[1 0 0],'LabelVerticalAlignment','bottom');
        
        % Patches (+- 1SD around the peak)
        PosRectinTF = [PosinTF-SDPeakGFP.(Fields{m}) 0 2*SDPeakGFP.(Fields{m}) max(MeanPeakGFP.(Fields{m}))];
        rectangle('Position',PosRectinTF,'FaceColor',[1 0 0 0.2],'EdgeColor',[1 0 0 0.2]);
        
        % 3) TOPOPLOTS
        subplot(3,length(Fields),m+2*length(Fields)); 
        Low = round(MaxGFP.(Fields{m}).Pos - SDPeakGFP.(Fields{m}));
        High = round(MaxGFP.(Fields{m}).Pos + SDPeakGFP.(Fields{m}));
        topoplotIndie(mean(EEGTEMP(Low:High,:),1),Chanlocs);
        title(sprintf('Averaged topography between %d and %d TF',Low,High));
        caxis([MinValVolt MaxValVolt]) % Same colorbar limits for all graphs
        colorbar
        
        % Save results in .txt file for each level        
        fprintf(fid,'\r\n%s\r\n',sprintf('LEVEL %d: %s',m,Fields{m}));
        fprintf(fid,'%s\r\n',['Individual GFP standard deviation (in TF) = ' num2str(SDPeakGFP.(Fields{m}))]);
        fprintf(fid,'%s\r\n',['Group-/Condition-averaged max GFP (in TF) = ' num2str(MaxGFP.(Fields{m}).Pos)]);
        
        %% Retrieving individual GFP averaged over the POI identified above
        % Write in the .txt output
%         fprintf(fid2,'\r\n%s\r\n',['individual GFP averaged over the POI: ' ...
%             num2str(MaxGFP.(Fields{m}).Pos) ' ' num2str(SDPeakGFP.(Fields{m})) ...
%             ' (mean +- 1SD)']);
%         fprintf(fid2,'\r\n%s\r\n',sprintf('LEVEL %d: %s',m,Fields{m}));
        
        % Compute individual average GFP over the POI (+- 1 SD)
        MaxGFP.MeanRange = [MaxGFP.MeanRange mean(GFPData.(Fields{m})(Low:High,:))];
        MaxGFP.SDRange = [MaxGFP.SDRange std(GFPData.(Fields{m})(Low:High,:))];    
        
        % Compute individual average voltage amplitude over the POI (+- 1
        % SD) and electrodes cluster of interest
        TEMPDat = squeeze(mean(EEGData.(Fields{m})(Low:High,ElectClust,:),2)); 
        MaxVolt.MeanRange = [MaxVolt.MeanRange mean(TEMPDat)];
        MaxVolt.SDRange = [MaxVolt.SDRange std(TEMPDat)];
        
        % Fill in the .txt output
%         ToSave = [AllNames;num2cell(MaxGFP.MeanRange);num2cell(MaxGFP.SDRange)]';
%         for j=1:length(ToSave);fprintf(fid2,'%s - %.4f %.4f\r\n',ToSave{j,1},ToSave{j,2}, ToSave{j,3});end
%         fprintf(fid2,'\r\n');
    end
    
    % Saving .mat file with results
    save([OutputPath '\' OutputFolder '\GFPResults_' CompN{1} '_' Method '_' date_name '.mat'],'MaxGFP',...
        'MeanPeakGFP','MeanGFP','SDPeakGFP','MaxGFPList')
    
    % Save in outputs in table and export to .xlsx file
    MaxGFPms = MaxGFPList/(SamplingRate/1000)+Epoch(1);
    TabOut = array2table([Unblinding num2cell([MaxGFPList MaxGFPms MaxGFP.MeanRange' MaxGFP.SDRange' ...
        MaxVolt.MeanRange' MaxVolt.SDRange'])]);
    TabOut.Properties.VariableNames = {'Files','Levels','MaxGFPTF','MaxGFPms',...
        'IndivGFPMeanOverPOI','IndivGFPSDOverPOI','IndivVoltageMeanOverPOIClust','IndivVoltageSDOverPOIClust'};
    writetable(TabOut,[OutputPath '\' OutputFolder '\GFPResults' '_' Method '_' date_name '.xlsx'],'Sheet',CompN{1})
    
    % Save figure 
    FigName = [OutputPath '\' OutputFolder '\GFPPeak_' CompN{1} '_' Method '_' date_name '.png'];
    fprintf(fid,'\r\n%s\r\n',['Related figure name: ' FigName]);
    saveas(Fig,FigName); close gcf; clear Fig;    
    
%     % Saving files and their respective component's peak GFP (in TF)
%     fprintf(fid,'\r\n%s\r\n','Filesnames and their respective component peak GFP (in TF)');
%     ToSave = [AllNames' num2cell(MaxGFPList)];
%     for j=1:length(ToSave);fprintf(fid,'%s - %d\r\n',ToSave{j,1},ToSave{j,2});end
%     fprintf(fid,'\r\n');
     clear MaxGFPList;
end

% Close txt file
fclose(fid);
% fclose(fid2);
