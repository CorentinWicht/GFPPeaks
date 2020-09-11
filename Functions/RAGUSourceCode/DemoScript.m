%% Tabula rasa in Matlab and an informed user on the console :-)
clearvars
clearvars -global
close all
clc

% The following script demonstrates the currently implemented script
% capacities of Ragu. The idea is to permit the user to run the
% time-consuming parts of the analysis, namely the randomization stats in
% an unsupervised mode.
%
% The figure handle 'RaguHandle' is used to keep a hand on all the data and
% is always the argument following the particular function to be invoked.
% Just keep this as it is... 

%% Open Ragu
%-------------------------------------------------------------------------
RaguHandle = Ragu;

% Load a particular Ragu file
%-------------------------------------------------------------------------
DataFileName = 'E:\Dropbox (PUK)\Docs\MATLAB\Ragu\Demo\Sample Maria Normalized.mat';
Ragu('LoadFile',RaguHandle,DataFileName);

%% Clear all previously applied filters and baseline corrections
%-------------------------------------------------------------------------
Ragu('UndoDataProcessing',RaguHandle);

% Apply a baseline correction
%-------------------------------------------------------------------------
StartTime =   0; % in milliseconds
EndTime   = 200;

Ragu('RemoveBaseline',RaguHandle,StartTime,EndTime);

%% Filter the data
%-------------------------------------------------------------------------
LowCut      =  1; % in Hertz
HighCut     = 30; % if a particular filter is not required, make the 
Notch       = 50; % corresponding variable empty (e.g. Notch = [];)

FilterOrder = 4;

Ragu('FilterData',RaguHandle,LowCut,HighCut,FilterOrder,Notch);

%% Now, we set all the options for the randomization part
%-------------------------------------------------------------------------
nRuns      = 10;  % The number of randomization runs
pThreshold = 0.05;  % The p-threshold (for display purposes only)
Normalize  = false; % Whether to L2 normalize the data for the TANOVA
NoFactXing = true;  % Set whether factor crossings are permitted when
                    % the data is randomized
Ragu('SetRandomizationOptions',RaguHandle,nRuns,pThreshold,Normalize,NoFactXing);

%% Now we set the analysis time period
%-------------------------------------------------------------------------
StartTime = []; % This is the onset of the analysis window. An empty matrix 
                % sets the start time to the beginning of the data

EndTime   = []; % This is the ffset of the analysis window. An empty matrix 
                % sets the end time to the end of the data

DoAverage = false;  % Sets whether the computations shall be done on the
                    % average of the time period, or time-point-wise

Ragu('SetAnalysisWindow',RaguHandle,StartTime,EndTime,DoAverage);
                    
%% This does the TCT
%-------------------------------------------------------------------------
OutputPDF = 'e:\Dropbox (PUK)\Desktop\TCTResult.pdf';
Ragu('ComputeTCT',RaguHandle,OutputPDF)
                    
%% This is the computation of the TANOVA                    
%-------------------------------------------------------------------------
OutputPDF = 'e:\Dropbox (PUK)\Desktop\TANOVAResult.pdf';
Ragu('ComputeTanova',RaguHandle,OutputPDF);

%% This is the computation of the GFP stats
%-------------------------------------------------------------------------
OutputPDF = 'e:\Dropbox (PUK)\Desktop\GFPResult.pdf';
Ragu('ComputeGFPStats',RaguHandle,OutputPDF);

%% Finally, we save the Ragu file, including all the results
FileToSave = 'e:\Dropbox (PUK)\Desktop\RaguDataWithResults.mat';
Ragu('SaveFile',RaguHandle,FileToSave);

%% And we close the whole thing
Ragu('Byebye',RaguHandle);