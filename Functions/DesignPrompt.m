%% PROMPT
clear variables

% DESIGN
if strcmpi(answer{4},'n')
    % Default values for the GUI below 
    to_display = [{'B';'Group';'OH_OH';'PBO_OH'} cell(4,3)];

    % Select folders on which to apply analyses
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

    % Saving the averaging parameters into a structure
    uisave('DesignList','TimeFreq_Design.mat')   
else
    % Open the design parameters
    uiopen('TimeFreq_Design.mat')
end

% Frequency bands definition
if strcmpi(answer{5},'n')
    % Default values for the GUI below 
    to_display = [[{'Theta';'Beta'};cell(10,1)] [{'4 7';'13 30'};cell(10,1)]];

    % Select folders on which to apply analyses
    ScreenSize=get(0,'ScreenSize');
    f = figure('Position', [ScreenSize(3)/2-500/2 ScreenSize(4)/2-500/2 500 500]);
    p = uitable('Parent', f,'Data',to_display,'ColumnEdit',true,'ColumnName',...
        {'Bands', 'Boundaries'},'CellEditCallBack',...
        'BandsList = get(gco,''Data'');');
    uicontrol('Style', 'text', 'Position', [20 350 400 150], 'String',...
            {['FREQUENCY BANDS DEFINITION' newline ''],...
            ['The definition of the frequency bands needs to be structured in the following way:'...
            newline ...
            newline '1st column = Name of the frequency bands'...
            newline '2nd column = Low and high frequency boundaries'...
            newline]});
    % Wait for t to close until running the rest of the script
    waitfor(p)

    % If no modifications to the example in the figure
    if ~exist('BandsList','var')
        BandsList = to_display;
    end

    % Saving the averaging parameters into a structure
    uisave('BandsList','TimeFreq_Bands.mat')   

else
    % Open the design parameters
    uiopen('TimeFreq_Bands.mat')
end

% Optional: Importing TF decomposed data
ImportTF = questdlg(['OPTIONAL' newline newline 'Would you like to import time-frequency (TF) decomposition matrix?'...
    newline 'This will bypass the first half of the script (computation of TF decomposition), which is very slow.'], ...
	'Optional TF decomposition import', ...
	'Yes','No','No');

%% Saving all data in mat file
clear f p ScreenSize

save([pwd '\TimeFreq_Param.mat'])