function varargout = RaguCollapseXAxis(varargin)
% RAGUCOLLAPSEXAXIS MATLAB code for RaguCollapseXAxis.fig
%      RAGUCOLLAPSEXAXIS, by itself, creates a new RAGUCOLLAPSEXAXIS or raises the existing
%      singleton*.
%
%      H = RAGUCOLLAPSEXAXIS returns the handle to a new RAGUCOLLAPSEXAXIS or the handle to
%      the existing singleton*.
%
%      RAGUCOLLAPSEXAXIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGUCOLLAPSEXAXIS.M with the given input arguments.
%
%      RAGUCOLLAPSEXAXIS('Property','Value',...) creates a new RAGUCOLLAPSEXAXIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RaguCollapseXAxis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RaguCollapseXAxis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RaguCollapseXAxis

% Last Modified by GUIDE v2.5 05-Mar-2019 14:06:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RaguCollapseXAxis_OpeningFcn, ...
                   'gui_OutputFcn',  @RaguCollapseXAxis_OutputFcn, ...
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


% --- Executes just before RaguCollapseXAxis is made visible.
function RaguCollapseXAxis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RaguCollapseXAxis (see VARARGIN)

% Choose default command line output for RaguCollapseXAxis
handles.output = hObject;

txtX = 'ms';
data = num2cell(zeros(3,2));
for i = 1:size(data,1)
    data{i,3} = txtX;
end

set(handles.uitable1,'Data',data,'ColumnName',{'Start','End','Unit'});

guidata(hObject, handles);

% UIWAIT makes RaguCollapseXAxis wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RaguCollapseXAxis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in ButtonAdd.
function ButtonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.uitable1,'Data');
nRows = size(data,1);
data(nRows+1,1) = {data{nRows,2}+1};
data(nRows+1,2) = {data{nRows+1,1}+1};
data(nRows+1,3) = data(nRows,3);
set(handles.uitable1,'Data',data);

% --- Executes on button press in ButtonLess.
function ButtonLess_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.uitable1,'Data');
if size(data,1) > 1
    data(end,:) = [];
    set(handles.uitable1,'Data',data);
end


% --- Executes on button press in ButtonOK.
function ButtonOK_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud.Latencies = get(handles.uitable1,'Data');
ud.Latencies(:,3) = [];
ud.NewLabel  = get(handles.editNewLabel,'String');
set(handles.output,'Userdata',ud);
uiresume(handles.figure1);


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function editNewLabel_Callback(hObject, eventdata, handles)
% hObject    handle to editNewLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNewLabel as text
%        str2double(get(hObject,'String')) returns contents of editNewLabel as a double


% --- Executes during object creation, after setting all properties.
function editNewLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNewLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
