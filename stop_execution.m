function varargout = stop_execution(varargin)
% STOP_EXECUTION M-file for stop_execution.fig
%      STOP_EXECUTION, by itself, creates a new STOP_EXECUTION or raises the existing
%      singleton*.
%
%      H = STOP_EXECUTION returns the handle to a new STOP_EXECUTION or the handle to
%      the existing singleton*.
%
%      STOP_EXECUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STOP_EXECUTION.M with the given input arguments.
%
%      STOP_EXECUTION('Property','Value',...) creates a new STOP_EXECUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stop_execution_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stop_execution_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stop_execution

% Last Modified by GUIDE v2.5 07-May-2006 16:39:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stop_execution_OpeningFcn, ...
                   'gui_OutputFcn',  @stop_execution_OutputFcn, ...
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


% --- Executes just before stop_execution is made visible.
function stop_execution_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stop_execution (see VARARGIN)

% Choose default command line output for stop_execution
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stop_execution wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stop_execution_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


