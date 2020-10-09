%% GFP PEAK IDENTIFICATION
% This script automatically identifies the group-averaged GFP peak inside a
% period of interest (POI)

% Corentin Wicht
% Laboratory for Neurorehabilitation Science, UNIFR, CH
% 02.06.2020

% TO DO:
% - The variable "TempElect" should change for each component!
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
    '1) semi-automatic' newline '2) manual']},'Settings',1,{'1','1'});
PromptIndivGFP = str2double(PromptIndivGFP);

% Extension
if PromptIndivGFP(1)==1; Extension='.ep'; elseif PromptIndivGFP(1)==2; Extension='.eph'; end

% Method
if PromptIndivGFP(2)==1; Method='semi-automatic'; elseif PromptIndivGFP(2)==2; Method='manual'; end

% Path of your upper folder containing your .EP data
root_folder = uigetdir(p2,'Choose the path of your most upper folder containing your files.');
cd(root_folder)
FileList = dir(['**/*' Extension]);
OutputPath = strsplit(root_folder,'\');
OutputPath = strcat(OutputPath(1:end-1),'\');
OutputPath = strcat(OutputPath{:});

% Sort FileList according to natural order
[~, Idx] = natsort({FileList.name});
FileList = FileList(Idx);

% Define design
% Default values for the GUI below 
to_display = [{'B';'Group';'OH_OH';'PBO_OH'} cell(4,3)];

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
        newline '3rd/4th lines = Name of the levels'...
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

% Reshaping the design structure
Design=[];
Factors = [0 0];
for k=1:size(DesignList,2)
    if strcmpi(DesignList{1,k},'W')
        Factors(1) = Factors(1) + 1;
        Design.Within.(sprintf('Factor%d',Factors(1))){k}=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    elseif strcmpi(DesignList{1,k},'B')
        Factors(2) = Factors(2) + 1;
        Design.Between.(sprintf('Factor%d',Factors(2))){k}=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    end
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
for k=1:length(Fields) % For BS and/or WS category
    Factors = fieldnames(Design.(Fields{k}));
    for j=1:length(Factors) % For each factor
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
to_display = cell(4,3);

% Define Components and bounds
ScreenSize=get(0,'ScreenSize');
f = figure('Position', [ScreenSize(3)/2-500/2 ScreenSize(4)/2-500/2 600 600]);
p = uitable('Parent', f,'Data',to_display,'ColumnEdit',true(1,size(to_display,2)),'ColumnName',...
    {'Name','Lower bound (ms)','Upper bound (ms)'},'Position',[50 -100 500 400],...
    'CellEditCallBack','CompList = get(gco,''Data'');');
uicontrol('Style', 'text', 'Position', [60 350 450 130], 'String',...
        {['COMPONENTS DEFINITION' newline ''],['The definition of the Period(s) of Interest needs to be structured in the following way:'...
        newline '1st column = Name of the Component (e.g. N2)',...
        newline '2nd column = Corresponding lower bound (e.g. 200ms post-stim.)',...
        newline '3rd column = Corresponding higher bound (e.g. 350ms post-stim.)']});
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
        
for n=1:sum(~cellfun(@(x) isempty(x), CompList(:,1))) % For each Comp
    
    % Initialize variables
    Pos = 1; Numbers = repmat({[]},[size(AllNames,2) 1]);
    
    % Remove empty lines and store data for the current component
    EmptIdx = ~cellfun(@(x) isempty(x),CompList(n,:));
    CompN = CompList(n,EmptIdx); 
    
    % Component
    TEMPComp = cellfun(@(x) str2double(x),CompN(2:3)); % Load current Component bounds
    if any(Epoch<0)
        TEMPComp = TEMPComp + abs(Epoch(1)); % Move based on prestim
    end
    TEMPCompinTF = round(TEMPComp/1000*SamplingRate); % Convert to TimeFrames and round to nearest integer
    
    for k=1:length(Fields) % For each levels
        for m=1:size(EEGData.(Levels{k}),3) % For each subject
            EEGTEMP = squeeze(EEGData.(Levels{k})(:,:,m));
            GFPData.(Levels{k})(:,m)=computegfp(EEGTEMP);

            %% PEAK IDENTIFICATION ON INDIVIDUAL GFP
            if strcmpi(Method,'semi-automatic')
                % Look for maximum/minimum inside Component bounds
                TEMPGFP = GFPData.(Levels{k})(:,m);
                MaxVal = max(TEMPGFP(TEMPCompinTF(1):TEMPCompinTF(2)));
                MaxPos = find(TEMPGFP==MaxVal);
                
                Fig = figure('units','normalized','outerposition',[0 0 1 1]); % Initialize figure in full screen
                File = strsplit(strrep(Names.(Levels{k}){m},'_','-'),'.'); 
                Title = [CompN{1} ': ' Fields{k} ' ' File{1}]; % Main title (component: level_subject)
                  annotation('textbox',[.4 .5 .5 .5],'String',Title,'FitBoxToText',...
                  'on','FontSize',20,'EdgeColor','w');
                
                % 1. Electrodes level data
                if any(EmptIdx); TempElect = [CompN(end) CompN(end)]; 
                else; TempElect = CompN(end-2:end); 
                end
                
                subplot(2,3,[1 2]); 
                plot(EEGTEMP, 'Color','k'); 
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
            
                % 2. GFP
                subplot(2,3,[4 5]); plot(TEMPGFP, 'Color','k'); 
                title('GFP'); % axis tight;
                axis tight; Yl=ylim(gca); % retrieve auto y-limits
                xlabel('Time (ms)'), ylabel('GFP');
                set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr)); 

                % Vertical line
                xline(MaxPos,'-',LabelValue,'Color','b','LabelVerticalAlignment','bottom');

                % Patches
                PosRectinTF = [TEMPCompinTF(1) Yl(1) TEMPCompinTF(2)-TEMPCompinTF(1) Yl(2)*0.1]; 
                rectangle('Position',PosRectinTF,'FaceColor',[0 0 1 0.2],'EdgeColor',[0 0 1 0.2]);

                % Legend
                LegendLabel = {'GFP';CompN{1}};
                Lgd = legend(LegendLabel); title(Lgd,'Component','Color','k');
                
                % Store parameters
                Param.Ticks=Ticks;
                Param.XTStr=XTStr;
                Param.TEMPCompinTF=TEMPCompinTF;
                Param.Epoch=Epoch;
                Param.SamplingRate=SamplingRate;
                Param.CompN=CompN{1};
                
                % 3. TOPOPLOTS
                slider_plot(Fig,EEGTEMP,MaxPos,Chanlocs,Param,TEMPGFP)
                
                % RESPONSE TABLE
                % Matrix to integrate in the following uitable
                Numbers{Pos} = MaxPos;
                to_display = [AllNames',Numbers];

                 % Select folders on which to apply analyses
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
                MaxIndivGFP.CompN{1}.(Levels{k})(m) = MaxGFPList{Pos,2}; 
                to_display{m,2} = MaxGFPList{Pos,2}; Numbers{Pos} = MaxGFPList{Pos,2};
                Pos = Pos + 1;
                clear Fig MaxGFPList;
            end
        end
    end
    
    % Compute mean GFP for each level
    MeanGFP = structfun(@(x) mean(x,2),GFPData,'UniformOutput',0);
    if strcmpi(Method,'semi-automatic')
        TEMP = cell2mat(struct2cell(MaxIndivGFP.CompN{1}));
        MaxGFPList = TEMP(:);
        MeanPeakGFP = structfun(@(x) mean(x,2),MaxIndivGFP.CompN{1},'UniformOutput',0);

        % Compute SD over individual GFP peaks for each level
        SDPeakGFP = structfun(@(x) std(x),MaxIndivGFP.CompN{1},'UniformOutput',0);
    end

    % MANUAL INDIVIDUAL GFP PEAKS LOCALISATION (e.g. using CARTOOL)
    if strcmpi(Method,'Manual')

        % Matrix to integrate in the following uitable
        to_display = [AllNames',Numbers];

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
        
        % Convert from TF to ms (round to nearest integer)
        MaxGFPList = cellfun(@(x) str2double(x),MaxGFPList(:,2));
        OutOfBounds = AllNames(MaxGFPList<TEMPCompinTF(1) | MaxGFPList>TEMPCompinTF(2))';
        MaxGFPListMS = round(MaxGFPList/1024*1000 + Epoch(1));
        
        % Warning for files that are out of component's bounds  
        warning('The following files are out of component %s bound:',CompN{1});
        fprintf('%s \n',OutOfBounds{:});
        
        % Build structure by levels
        for k=1:length(Fields) 
            MaxGFPStructTF.(Fields{k}) = MaxGFPList(contains(AllNames,Levels{k}));
        end
        
        % Compute mean GFP for each level
        MeanPeakGFP = structfun(@(x) mean(x),MaxGFPStructTF,'UniformOutput',0);

        % Compute SD over individual GFP for each level
        SDPeakGFP = structfun(@(x) std(x),MaxGFPStructTF,'UniformOutput',0);
    end

    %% PEAK IDENTIFICATION ON MEAN GFP
    Fig = figure('units','normalized','outerposition',[0 0 1 1]); % Initialize figure in full screen
    annotation('textbox',[.45 .5 .5 .5],'String',...
        sprintf('Results for component %s',CompN{1}),'FitBoxToText',...
        'on','FontSize',20,'EdgeColor','w');
    
    % For .txt output
    fprintf(fid,'\r\n%s\r\n',['----' CompN{1} ' COMPONENT [' ...
        num2str(TEMPCompinTF(1)) '-'  num2str(TEMPCompinTF(2)) ' TF]'  '----']);
    
    for m=1:length(Fields) % For each level
        MaxGFP.(Fields{m}).Pos = round(MeanPeakGFP.(Fields{m}));
        MaxGFP.(Fields{m}).Value = MeanGFP.(Fields{m})(MaxGFP.(Fields{m}).Pos);
    
        %% PLOTTING
        
        % Main plot and parameters
        subplot(2,2,m)
        plot(MeanGFP.(Levels{m}),'Color','k');axis tight;
        ylim([0 max(MeanGFP.(Levels{m}))])
        xlabel('Time (ms)'), ylabel('Mean GFP')
        title(sprintf('%s',strrep(Fields{m},'_','-'))) 
        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
        
        % Vertical line (Peak)
        PosinMS = round(MaxGFP.(Fields{m}).Pos/SamplingRate*1000);
        LabelValue = [num2str(PosinMS + Epoch(1)) 'ms ' '[' num2str(MaxGFP.(Fields{m}).Value) ']'];
        xline(MaxGFP.(Fields{m}).Pos,'-',LabelValue,'Color',[1 0 0],'LabelVerticalAlignment','bottom');
        
        % Patches (+- 1SD around the peak)
        PosRectinTF = [PosinMS-SDPeakGFP.(Fields{m}) 0 2*SDPeakGFP.(Fields{m}) max(MeanPeakGFP.(Fields{m}))];
        rectangle('Position',PosRectinTF,'FaceColor',[1 0 0 0.2],'EdgeColor',[1 0 0 0.2]);
       
        % TOPOPLOTS
        subplot(2,2,m+2); 
        Low = round(MaxGFP.(Fields{m}).Pos - SDPeakGFP.(Fields{m}));
        High = round(MaxGFP.(Fields{m}).Pos + SDPeakGFP.(Fields{m}));
        topoplotIndie(mean(EEGTEMP(Low:High,:),1),Chanlocs);
        title(sprintf('Averaged topography between %d and %d TF',Low,High));
        
        % Save results in .txt file for each level        
        fprintf(fid,'\r\n%s\r\n',sprintf('LEVEL %d: %s',m,Fields{m}));
        fprintf(fid,'%s\r\n',['Individual GFP standard deviation (in TF) = ' num2str(SDPeakGFP.(Fields{m}))]);
        fprintf(fid,'%s\r\n',['Group-/Condition-averaged max GFP (in TF) = ' num2str(MaxGFP.(Fields{m}).Pos)]);
    end
    
    % Save figure 
    FigName = [OutputPath '\' OutputFolder '\GFPPeak_' CompN{1} '_' Method '_' date_name '.png'];
    fprintf(fid,'\r\n%s\r\n',['Related figure name: ' FigName]);
    saveas(Fig,FigName); close gcf; clear Fig;    
    
    % Saving files and their respective component's peak GFP (in TF)
    fprintf(fid,'\r\n%s\r\n','Filesnames and their respective component peak GFP (in TF)');
    ToSave = [AllNames' num2cell(MaxGFPList)];
    for j=1:length(ToSave);fprintf(fid,'%s - %d\r\n',ToSave{j,1},ToSave{j,2});end
    fprintf(fid,'\r\n');
    clear MaxGFPList;
end

% Close txt file
fclose(fid);