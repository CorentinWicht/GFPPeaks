function varargout = Ragu_Import_Help(varargin)
% RANDOMIZER_HELP M-file for Ragu_Import_Help.fig
%      RANDOMIZER_HELP, by itself, creates a new RANDOMIZER_HELP or raises the existing
%      singleton*.
%
%      H = RANDOMIZER_HELP returns the handle to a new RANDOMIZER_HELP or the handle to
%      the existing singleton*.
%
%      RANDOMIZER_HELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RANDOMIZER_HELP.M with the given input arguments.
%
%      RANDOMIZER_HELP('Property','Value',...) creates a new RANDOMIZER_HELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_Import_Help_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_Import_Help_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_Import_Help

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Last Modified by GUIDE v2.5 25-Nov-2014 12:56:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_Import_Help_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_Import_Help_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Ragu_Import_Help is made visible.
function Ragu_Import_Help_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_Import_Help (see VARARGIN)

% Choose default command line output for Ragu_Import_Help
%get(hObject)
handles.output = hObject;

type = varargin{1};

switch type
    case 1
        set(handles.HelpText,'String',HelpTextGeneral());
   case 2
        set(handles.HelpText,'String',HelpTextImport());
   case 3
        set(handles.HelpText,'String',HelpTextDesign());
   case 4
        set(handles.HelpText,'String',HelpTextLoadSave());
   case 5
        set(handles.HelpText,'String',HelpTextOptions());
   case 6
        set(handles.HelpText,'String',HelpTextAnalysis());
   case 7
        set(handles.HelpText,'String',HelpTextView());
        
        
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ragu_Import_Help wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Ragu_Import_Help_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function HelpText_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function HelpText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt = HelpTextImport()

txt = {...
''...
'                                             Importing your data'...
'                                             -------------------'...
'Importing the EEG/ERP data',...
'--------------------------',...
'Before you can go into the analysis of your data, you have to import all the EEG/ERP data and specify the design of the experiment. Before you can import anything, make sure you have all your EEG/ERP data organized as follows:'...
'- All data is in the same subdirectory. Optimally, this subdirectory contains only the data that you want to import,'...
'  nothing else.'...
'- For each subject and condition, there is an ASCII text file with the data. Each line of this file contains'...
'  the measurements of all channels at one moment in time, each column in the file contains all moments of time'...
'  of one channel. The program used the Matlab function �load� to read the file. In case of trouble, look up Matlab�s'...
'  documentation on the function �load�. Filenames should not contain spaces.'...
'- All files of a single subject should have a common part in the name that uniquely identifies the subject.'...
'- All files of a single condition should have a common part in the name that uniquely identifies the condition.'...
'- All files have the same amount of channels and time-points, and the channels are always in the same order.'...
''...
'Once you have organized the data as described above, you can proceed to import the data. In the dialog that opens when you click the �Import EEG/ERP data� button, you have to set the following parameters.'...
'- The subdirectory where all the data is located. Use the �Change directory� button to navigate there. You can have',...
'  a look at the content of the directory by clicking the �View directory� button.'...
'- The search mask. Enter a string with Windows-wildcarts in it that identifies the files of one condition for all'...
'  subjects. To illustrate this, lets assume you have data from 5 subjects, and 3 conditions. The data of subject 1',...
'  condition A is in a file called �Sub1_Cond_A.asc�, condition B of the same subject is in file �Sub1_Cond_B.asc�',...
'  condition B of subject 2 is in a file called �Sub2_Cond_A.asc�, and so on. Wildcart that identify the files of all',...
'  subjects for one condition would for example be �*Cond_A.asc�, or �*Cond_B.asc�. The tool is based on Matlab�s �dir�',...
'  command, so check their documentation in case you run into problems or have further questions.'...
'- Tags for the different conditions. Enter the part of the filenames that uniquely identifies a single condition in the'...
'  field �new tag� and press �Add�. In the above example, this would e.g. be �Cond_A�. Repeat until all conditions that'...
'  you want to import are listed. In the above example, we�d have to add �Cond_A�, �Cond_B� and �Cond_C�.'...
'- Click Go. Based on the information that you have provided, the program is now able to construct all the filenames it'...
'  can expect and load these files.'...
''...
'When you have successfully imported the ERP data, you can have a glance at your date using the �View EEG/ERP data� button. You should get a list of the files that you have imported, ordered as a two-dimensional subject by condition matrix. Double-clicking on a filename displays the ERP data of that subject and condition in a butterfly plot.'...
'',...
'Defining the montage',...
'--------------------',...
'You can optionally define the montage of your data, which will later allow you to better visualize the results. When you work with Analyzer and you have my Matlab plugin, you can simply save the Channel structure to a Matlab file and load this file as your montage. Alternatively, you can load a text file with the xyz coordinates of the electrodes. (x<0 is left x > 0 is right, y>0 is anterior, y < 0 is posterior, z > 0 is superior). The file format is identical with the SXYZ files of the sLORETA package. It has the number of channels on the first line followed by a line of information for each channel. This line contains the x, y, and z coordinates of each channel, and the electrode label. Channels must be in the same order as in the data.'...
''...
'Here is an example:'...
'70'...
'-26.81 84.06 -10.56 Fp1'...
'29.41 83.74 -10.04 Fp2'...
'-48.05 51.87 39.87 F3'...
'...'...
''...
'Data options'...
'------------'...
'Specify the interval between adjacent points on the x-axis (typically the sampling interval) and the onset of your data on the x-axis to get proper axes in your displays. You can also set the units on the x-axis.'...
''...
'Inspecting your data'...
'--------------------'...
'When you click on "View" EEG/ERP menu, you can inspect whether your data was correctly loaded. You will see a list of all files that you have loaded. Mark one or several datasets you which to inspect, and click the "Show" button. Ragu will show you a butterfly plot of the mean of the data you have selected. When you click on the butterfly plot, you obtain a map of the data.'...
''...
};

function txt = HelpTextDesign()

txt = {...
'',...
'                                       Define the design of the experiment'...
'                                       -----------------------------------'...
''...
'Specifying the within subject design'...
'------------------------------------'...
'If you convinced yourself that the data is correctly imported, you can proceed to specify your design. In the dialog that opens when you hit the �Within subject design� button, first choose whether you have one or two within subject factors. Then you proceed to assign the data that you have imported to different levels of your factors. To specify the levels of factor 1 click the F1 button. To specify the levels of factor 2, click the F2 button. To set the level of a condition for a given factor, select the condition in the list-box and use the �+� and �-� buttons to modify the level. In the subsequent analysis of the effect of a given factor, conditions that have the same level will be assumed to have the same characteristics of this factor. Let�s illustrate this. Assume you have recorded ERP data with rare and frequent stimuli in a series of subjects. Each subject performed the experiment once under drug A, once under drug B and once under drug C. We have thus 6 conditions, Rare_DrugA, Rare_DrugB, Rare_DrugC, Frequent_DrugA, Frequent_DrugB and Frequent_DrugC. There are thus a factor drug, and a factor stimulus type. Let�s make drug factor 1 and stimulus type factor 2. To specify the design, we thus have to do the following. Select factor one. Set all conditions belonging to drug A (e.g. Rare_DrugA and Frequent_DrugA) to the same level (i.e. to 1). Set all conditions belonging to drug B to another level (i.e. to 2). Set all all conditions belonging to drug C to a third level (i.e. to 3). Factor one is now defined. Click �F2� to specify the levels for the second factor. In our case set all conditions belonging the rare stimulus type (Rare_DrugA, Rare_DrugB and Rare_DrugC) to the same level of factor 2 (i.e. to 1).Finally, set all conditions belonging to the frequent stimulus type to another value (i.e. 2). The number you use to specify the levels are not important, but your design has to be orthogonal and complete. If there are no repeated measures, leave all conditions on the same level. You can also exclude certain conditions from the analysis using the exclude botton. '...
'There is a checkbox to indicate if the first within subject factor is rank/interval scaled. In this case, a within subject TANCOVA will be computed, Otherwise, a TANOVA will be computed, assuming a categorical within subject factor.'...
''...
'Specifying the between subject design'...
'-------------------------------------'...
'When you click on the �Subject Grouping / Behavior� button, you can define the feature that distinguishes the subjects. This can either be data on some rank or interval scale (e.g reaction time, or age), or some categorical information, like gender or diagnosis. The graph shows the subjects feature. To change the value of a subject, enter the value that to want to set in the value box, make sure you have �Update� checked and click on the subject or subjects you want. The subject�s feature is updated. To inspect the data, uncheck the �Update� box. As before, if the feature is categorical, all subjects having the same feature are assumed to belong to the same group, the values itself don�t matter. Click �Done� when you�re finished to store your work. If there are no groups, leave all subjects on the same level. You can also exclude subjects using the exclude checkbox.'... 
'There is a checkbox to indicate if the between subject data is rank/interval scaled. In this case, a between subject TANCOVA will be computed, Otherwise, a TANOVA will be computed, assuming a categorical between subject grouping factor.'...
''...
};

function txt = HelpTextLoadSave()

txt = {...
'',...
'                                           Storing data and results'...
'                                           ------------------------'...
''...
'Saving and loading of the data'...
'------------------------------'...
'Use the save and load buttons to save and load all the data. All data, design information and results are stored and can later be retrieved. To save different analyses based on the same data, load a file with the data and design, and save it under a different name after the analyses. All data is stored in Matlab format, so you can open, inspect, manipulate (and destroy) the data using Matlab.'...
'In addition, you have several other options to save data.'...
};

function txt = HelpTextOptions()

txt = {...
''...
'                                          Setting global parameters'...
'                                          -------------------------'...
'Randomization options',...
'---------------------'...
'Here, you can set the following options:'...
'- Number of randomizations: The more, the better, but the more the longer. If you�re not sure of what you�re doing, start'...
'  with small numbers (like 200), for publishable results, make it 5000.'...
'- Normalize: This makes sets all data to equal variance across electrode before the analysis. Useful if you want to'...
'  separate effect of topographic differences from effects of amplitude differences. If you normalize your data, it'...
'  it is a good practice to test for the normalization factor, i.e. the GFP.' ...
'- Threshold: Where to consider an effect as significant. This options influences the results of the post-hoc'...
'  TANOVA/TANCOVA hit statistics'...
'- No factor crossing: If in a two-factorial design, one effect has much larger effects than the other, a fully free'...
'  randomization will make it difficult for the smaller factor to become significant. You can restrain the possible'...
'  permutations such that each factor competes only against permutations of itself, and not of the other factor.'...
''...
'Clear console window'...
'--------------------'...
'The console window is useful to track the causes if something went wrong. If you encounter problems, have a look there for possible hints. You can clear this window to get rid of older messages.'...
''...
'Pretty maps'...
'--------------------'...
'There is a faster and a prettier way to display maps. Your choice...'...
};

function txt = HelpTextAnalysis()

txt = {...
''...
'                                             Analyzing your data'...
'                                             -------------------'...
'Analysis of topographic consistency'...
'-----------------------------------'...
'This analysis tests the null-hypothesis that the Global Field Power (GFP) of the data one condition / Group at one moment could have been observed if there was no consistent topography across subjects (Koenig et al., 2010). The resulting graph shows, for each group and condition the GFP of the Grand-Mean (black line), the significance of the test (gray areas), and the significance threshold chosen. If you have applied a duration threshold criterion (see Koenig et al., 2010 for computational details), periods that also meet the criterion are shown in green. You can click on the graph to get a momentary map. Using a right-click on this map gives you further options like t-maps against null.'...
''...
'TANOVA/TANCOVA'...
'--------------'...
'Establishes whether topographic differences between the levels of the conditions could have been observed by chance. Includes the interaction between the factors and group.'...
'The display shows the p-value, white areas indicate p-values below the threshold that you have set. All the ideas behind this are laid out in detail in our book Electrical Neuroimaging.'...
'For a detailed description of the functionality of the TANOVA result display, we refer to the 2011 article that came with the software.'...
''...
'GFP/RMS'...
'-------'...
'Tests for differences in the GFP, using the same mechanics as the TANOVA. This is an important additional step if you have normalized your data for the TANOVA.'...
''...
'Hit statistics'...
'--------------'...
'Count: Once you have run a TANOVA/TANCOVA, you can test whether the number of significant timer periods that you have found across the time interval may also have occurred by chance. The resulting histograms show the count of the number of significant time points under the null-hypothesis, the red dot indicates the number of significant time-points that you have actually encountered. The p-values are in the titles of the graphs. '...
''...
'Duration'...
'--------'...
'Similarly as above, you can ask how long a continuous period of a significant effect would last under the null-hypothesis (this may for example depend on the filters you have applied before the statistics). The result shows the probability (y-axis) of the duration of the of a significant effect (x-axis) under the null-hypothesis, the significance level (set in the randomization options) is indicated by the red line. The green line indicates above what duration the null-hypothesis that an effect of some duration has been observed by chance is rejected. When you compute these critical duration and redisplay the results of the TANOVA/TANCOVA, effects meeting this criterion ae shown in green.'...
'All the ideas behind this are laid out in detail in our book Electrical Neuroimaging.'...
''...
'Global p-AUC test'...
'-----------------'...
'This tests asks if the distribution of p-values across all tested time points is significantly smaller than one would expect under the null-hypothesis.'...
''...
'Microstate analysis'...
'-------------------'...
''...
'Compute microstate maps:'...
'- - - - - - - - - - - - '...
'Here, you can compute microstate maps and test for differences in microstate features as described in Koenig et al, 2014.'...
'The computation of microstate maps offers several options: '...
'- You can chose between the "old" k-means and faster and more robust AAHC algorithm (see Murray et al. Brain Topogr.'...
'  (2008) 20:249�264 for details). If you chose the k-means you have to set the number of re-initializations, since the'...
'  algorithm may get caught in a local minimum.'...
'- You can then set a fixed number of microstates to fit, or you can use the cross-validation procedure in Koenig et al.,'...
'  2014 to find the optimal number of classes. If you do so, you have to set the size of the training group (on how many'...
'  subjects each model shall be constructed), and how many cross validation runs you which to compute. You will receive a'...
'  display that shows the explained variance as function of number of classes in the training and test-set. The relevant'...
'  information is in the test set, search for the number of classes where this curve starts to show a plateau. This is a'...
'  good choice for the optimal number of microstates.'...
'- You can smooth the labelling of the microstate assignment across time. See Pascual-Marqui et al., IEEE Trans Biomed Eng.'...
'  1995; 42:658�65 for details.'...
''...
'Loading and saving microstate maps:'...
'- - - - - - - - - - - - - - - - - -'...
'Uses plain ASCII text files.'...
''...
'Microstate statistics'...
'- - - - - - - - - - -'...
'Computes microstate statistics as developped in Koenig et al., 2014.'...
''...
'Understanding the display'...
'- - - - - - - - - - - - -'...
'The microstate display shows on the right top area the obtained microstate maps. Below, you see the assignment of the data to the microstates, the color indicates class, the height indicates explained GFP and the black line is the total GFP of the data. You can click o the button below each map to select a particular microstate class for further analysis.'...
'On the left side, there are boxes for the different microstate parameters. By clicking on a cell of one of these boxes, the microstate assignment corresponding to that cell is shown on the right, with the corresponding values in the title of each graph. If you have already computed the statistics, these cells will contain the optained p-values. If a test cannot be computed because a condition does not show a microstate class, the cell shows NaN.'...
''...
'The features being computed are:'...
'- Onset: Time the class is observed first.'...
'- Offset: Time the class is observed last.'...
'- Duration: Total time the class is observed.'...
'- Area under the curve: What it says, for the given class.'...
'- Center of gravity: In time, for the given class.'...
'- Mean GFP: Also obvious.'...
};


function txt = HelpTextGeneral()


txt = {...
''...
'                                             General information'...
'                                             -------------------'...
''...
'Introduction                                                                          Ragu Version: November 2014',... 
'------------'...
'Ragu lets you compute TANOVAs and TANCOVAs of EEG and ERP scalp field data. It can accommodate two within- and one between-subject factors, each with potentially several levels, and will also test for interactions between these factors. If the between subject factor is rank or interval scaled, a TANCOVA can be computed. Ragu  is Matlab based and stores all results in Matlab format. You can therefore access all the results also directly with Matlab.'...
'Ragu  is free, but we offer only limited support and expect you to quote the suggested papers. We are open for suggestions to improve the software where we find it useful for a wider community. We take absolutely no responsibility for the results obtained with Ragu , and the software may not be fool-prove.',...
'',...
'References'...
'----------'...
'In general, quote'...
'- Koenig T, Melie-Garcia L. (2009) Statistical Analysis of Multichannel Scalp Field Data. In Michel, Koenig, Brandeis,'...
'  Gianotti, Wackermann: Elecrical NeuroImaging. Cambridge University Press'...
''...
'For the TANOVA, quote'...
'- Wirth M, Horn H, Koenig T, Razafimandimby A, Stein M, Mueller T, et al. (2008) The early context effect reflects'...
'  activity in the temporo-prefrontal semantic system �evidence from electrical neuroimaging of abstract and'...
'  concrete word reading. Neuroimage. 42:423-436.'...
''...
'For the TANCOVA, quote'...
'- Koenig T, Melie-Garcia L, Stein M, Strik W, Lehmann C. (2008) Establishing correlations of scalp field maps with other'...
'  experimental variables using covariance analysis and resampling methods. Clin Neurophysiol. 119:1262-70.'...
''...
'For the Topographic consistency test, quote'...
'-  Koenig T, Melie-Garc�a L (2010) A method to determine the presence of event-related fields using randomization'...
'   tests. Brain Topography, 23:233-242.'...
''...
'For the microstate statistics, quote'...
'- Koenig T, Stein M, Grieder M, Kottlow M (2014) A tutorial on data-driven methods for statistically assessing ERP' ...
'  topographies. Brain Topography, 27(1):72-83.'...
''...
};

function txt = HelpTextView()

txt = {...
''...
'                                                Viewing results'...
'                                                ---------------'...
''...
'The view menu lets you recall the result windows from previous analysis. See the help on the analysis menue to understand the displays.'...
''...
'Furthermore, there is an interface to compute t-maps and sLORETA images. Choose a title for the contrast and the groups and/or conditions you want to contrast into the positive and negative parts of the contrast. If you want to compute double differences, use a baseline. In the parameters, you can set the time window and the scaling of the resulting maps. When clicking the "Show t-map" or "Show p-map", a new window opens with all the relevant information. You can create as many of these contrast windows as you which, such that you can compare different results.'... 
''...
'The sLORETA part computes contrasts in the inverse space based on the sLORETA implementation of the KEY Institute, and relies on that software. If you specify the path to the sLORETA system files and the file with the pseudoinverse, you can do contrasts in the inverse.'...
''...
'In addition, you can switch on a data cursor to read out specific results.'...
};






% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
