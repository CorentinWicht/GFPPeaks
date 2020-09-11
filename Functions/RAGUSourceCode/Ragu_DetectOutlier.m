function varargout = Ragu_DetectOutlier(varargin)
% RAGU_DETECTOUTLIER MATLAB code for Ragu_DetectOutlier.fig
%      RAGU_DETECTOUTLIER, by itself, creates a new RAGU_DETECTOUTLIER or raises the existing
%      singleton*.
%
%      H = RAGU_DETECTOUTLIER returns the handle to a new RAGU_DETECTOUTLIER or the handle to
%      the existing singleton*.
%
%      RAGU_DETECTOUTLIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGU_DETECTOUTLIER.M with the given input arguments.
%
%      RAGU_DETECTOUTLIER('Property','Value',...) creates a new RAGU_DETECTOUTLIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_DetectOutlier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_DetectOutlier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_DetectOutlier

% Last Modified by GUIDE v2.5 15-Sep-2016 12:27:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_DetectOutlier_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_DetectOutlier_OutputFcn, ...
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


% --- Executes just before Ragu_DetectOutlier is made visible.
function Ragu_DetectOutlier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_DetectOutlier (see VARARGIN)

% Choose default command line output for Ragu_DetectOutlier
handles.output = hObject;

ud.V = varargin{4}.V;
if ~isempty(varargin{4}.Design)
    NotUsedConditions = isnan(varargin{4}.Design(:,1));
    ud.V(:,NotUsedConditions,:,:) = [];
end
ud.IndFeature = varargin{4}.IndFeature;
ud.OrgIndex = 1:numel(ud.IndFeature);
ud.UseMe = ~isnan(ud.IndFeature);
ud.Labels = varargin{4}.Names(:,1);
set(handles.Ragu_DetectOutlier,'UserData',ud);

UpdatePlot(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ragu_DetectOutlier wait for user response (see UIRESUME)
uiwait(handles.Ragu_DetectOutlier);

function UpdatePlot(handles)

ud = get(handles.Ragu_DetectOutlier,'UserData');

DataToWorkWith = ud.V(ud.UseMe,:,:,:);
IndFeatureToWorkWith = ud.IndFeature(ud.UseMe);
Dims = size(DataToWorkWith);

groups = unique(IndFeatureToWorkWith);

Centering = eye(Dims(1)) - 1/ Dims(1);

DataToWorkWith = Centering * reshape(DataToWorkWith,[Dims(1),prod(Dims(2:end))]);

cov = DataToWorkWith * DataToWorkWith';

[v,d] = eigs(cov,2);

axis(handles.MDSAxis);

lbl = {'*g','+y','oB','.r','xc','sm','dk','vw'};


if numel(groups) <= 8
    cla(handles.MDSAxis);
    hold on
    for i = 1:numel(groups)
        idx = IndFeatureToWorkWith == groups(i);
        plot(handles.MDSAxis,v(idx,1),v(idx,2),lbl{i});
    end
else
    plot(handles.MDSAxis,v(:,1),v(:,2),'*g');
end
    dummy = v(:,1:2);
mx = max(abs(dummy(:))) * 1.05;
axis(handles.MDSAxis,[-mx mx -mx mx]);
ud.points = v;
ud.mx = mx;
set(handles.MDSAxis,'Color',[0.5 0.5 0.5]);


set(handles.Ragu_DetectOutlier,'UserData',ud);


% --- Outputs from this function are returned to the command line.
function varargout = Ragu_DetectOutlier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles)
    varargout{1} = [];
else
    ud = get(handles.Ragu_DetectOutlier,'UserData');
    set(handles.output,'UserData',ud.UseMe);
    varargout{1} = handles.output;
end


% --- Executes on button press in ButtonSelect.
function ButtonSelect_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.Ragu_DetectOutlier,'UserData');
plot(ud.points(:,1),ud.points(:,2),'*g')
axis(handles.MDSAxis,[-ud.mx ud.mx -ud.mx ud.mx]);
set(handles.MDSAxis,'Color',[0.5 0.5 0.5]);

[x,y] = ginput(1);

coords(:,1) = ud.points(:,1) - x;
coords(:,2) = ud.points(:,2) - y;

[dummy,idx] = min(sum(coords.^2,2));
MappedIndex = ud.OrgIndex(ud.UseMe);
ud.ToBeRemoved = MappedIndex(idx);
set(handles.CaseName,'String',ud.Labels{MappedIndex(idx)});
hold on
plot(ud.points(idx,1),ud.points(idx,2),'or','MarkerFaceColor','r');
hold off

axis(handles.MDSAxis,[-ud.mx ud.mx -ud.mx ud.mx]);
set(handles.MDSAxis,'Color',[0.5 0.5 0.5]);


set(handles.ButtonRemove,'Enable','on');
%set(handles.ButtonSelect,'Enable','off');
set(handles.Ragu_DetectOutlier,'UserData',ud);



function CaseName_Callback(hObject, eventdata, handles)
% hObject    handle to CaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CaseName as text
%        str2double(get(hObject,'String')) returns contents of CaseName as a double


% --- Executes during object creation, after setting all properties.
function CaseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonRemove.
function ButtonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.Ragu_DetectOutlier,'UserData');
ud.UseMe(ud.ToBeRemoved) = false;
set(handles.ButtonRemove,'Enable','off');
%set(handles.ButtonSelect,'Enable','on');
set(handles.CaseName,'String','');
set(handles.Ragu_DetectOutlier,'UserData',ud);

UpdatePlot(handles);


% --- Executes on button press in ButtonDone.
function ButtonDone_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.Ragu_DetectOutlier);


% --- Executes on button press in ButtonAutoSelect.
function ButtonAutoSelect_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonAutoSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.Ragu_DetectOutlier,'UserData');
plot(ud.points(:,1),ud.points(:,2),'*g')
axis(handles.MDSAxis,[-ud.mx ud.mx -ud.mx ud.mx]);
set(handles.MDSAxis,'Color',[0.5 0.5 0.5]);
p_val = str2double(get(handles.Edit_p,'String'));
MappedIndex = ud.OrgIndex(ud.UseMe);

dist2 = mahal(ud.points,ud.points);

[n,p] = size(ud.points);

critDist = ACR(p,n,p_val);

[mx,idx] = max(dist2);

if mx > critDist
    ud.ToBeRemoved = MappedIndex(idx);
    set(handles.CaseName,'String',ud.Labels{MappedIndex(idx)});;
    hold on
    plot(ud.points(idx,1),ud.points(idx,2),'or','MarkerFaceColor','r');
    hold off

    axis(handles.MDSAxis,[-ud.mx ud.mx -ud.mx ud.mx]);
    set(handles.MDSAxis,'Color',[0.5 0.5 0.5]);

    set(handles.ButtonRemove,'Enable','on');
    set(handles.Ragu_DetectOutlier,'UserData',ud);
else
    msgbox('No (further) outlier detected');
end
    



function Edit_p_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_p as text
%        str2double(get(hObject,'String')) returns contents of Edit_p as a double


% --- Executes during object creation, after setting all properties.
function Edit_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
