% Copyright (c) 2008, Department of Information Engineering, University of Padova.
% All rights reserved.
% 
% This file is part of YourTubes.
% 
% YourTubes is free software: you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% YourTubes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with YourTubes.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
%  @date 28/04/2008
%  @author Filippo Zanella <filippo.zanella@dei.unipd.it>
%          Daniele Beninato <beninato@dei.unipd.it>
% ===================================================================================

function varargout = YourTubes(varargin)
% YOURTUBES M-file for YourTubes.fig
%      YOURTUBES, by itself, creates a new YOURTUBES or raises the existing
%      singleton*.
%
%      H = YOURTUBES returns the handle to a new YOURTUBES or the handle to
%      the existing singleton*.
%
%      YOURTUBES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in YOURTUBES.M with the given input arguments.
%
%      YOURTUBES('Property','Value',...) creates a new YOURTUBES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before YourTubes_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to YourTubes_OpeningFcn via varargin.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @YourTubes_OpeningFcn, ...
    'gui_OutputFcn',  @YourTubes_OutputFcn, ...
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

addpath(genpath(pwd));    % Add the sub-directory to Matlab search $PATH
global DEBUG; DEBUG = 0;  % Activate all debug lines


% --------------------------------------------------------------------
function YourTubes_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% --- Executes just before YourTubes is made visible.
global DEBUG;
if (DEBUG), fprintf('YourTubes_OpeningFcn\n'); end

splash('logo.bmp',3000);

% uiwait(handles.YourTubes); % UIWAIT makes YourTubes wait for user response (see UIRESUME)

handles.output = hObject;    % Choose default command line output for YourTubes

set(handles.slider1, 'Enable', 'off');
set(handles.slider2, 'Enable', 'off');
set(handles.slider3, 'Enable', 'off');
set(handles.edit11, 'Enable', 'off');
set(handles.edit22, 'Enable', 'off');
set(handles.edit33, 'Enable', 'off');
set(handles.pushbutton, 'Enable', 'off');
set(handles.pushbutton, 'Visible', 'off');

handles.hYourTubes = handles.YourTubes; 
handles.hScrollPanelFig = [];
handles.hScrollPanelIm = [];
handles.hOverviewPan = [];
handles.hMeasureFig = [];

handles.filePrint = [];
handles.zoomInButton = [];
handles.zoomOutButton = [];
handles.roi = [];
%handles.measurePipes = [];

sliderValue = get(handles.slider1,'Value');
set(handles.edit11, 'String', sliderValue);
handles.th1 = sliderValue;
sliderValue = get(handles.slider2,'Value');
set(handles.edit22, 'String', sliderValue);
handles.th2 = sliderValue;
sliderValue = get(handles.slider3,'Value');
set(handles.edit33, 'String', sliderValue);
handles.th3 = sliderValue;

guidata(hObject,handles);        % Update handles structure

toolbar = findall(hObject,'type','uitoolbar');
delete(toolbar);
createToolbar(hObject,handles);  % Create Toolbar

%save YourTubesHandles handles;

setappdata(hObject,'YourTubesProperties',get(hObject));
%%% List of appdata %%%
% measurefig - struct (IMG,width,height)
% tubesfig - struct (IMG,width,height)
% polyfig - struct (ROI,x,y,xi,yi,roix,roiy)
% YourTubesProperties
% IMG
% OverviewListeners
% apiSP
% numberPipes
% indet - struct (Xi,Yi,Xf,Yf)


% --------------------------------------------------------------------
function choiceImage(hObject, eventdata, handles, fileName) %#ok<INUSL>
global DEBUG;
if (DEBUG), fprintf('choiceImage\n'); end

hYourTubes = handles.hYourTubes;

IMG = imread(fileName); % Original image
setappdata(hYourTubes,'IMG',IMG);

hOriginalFig = createFigure('Original Image');
set(hOriginalFig,'Visible','off');
hOriginalIm = imshow(IMG,'InitialMagnification','fit');
hScrollPanel = imscrollpanel(hOriginalFig,hOriginalIm);
hScrollPanelFig = ancestor(hScrollPanel,'figure'); % = hOriginalFig
hScrollPanelIm = hOriginalIm;

% apiScrollPanel.replaceImage(IMG); % Image sobstitution
apiScrollPanel = iptgetapi(hScrollPanel);
setappdata(hYourTubes,'apiSP',apiScrollPanel);
apiScrollPanel.setMagnification(apiScrollPanel.findFitMag()); 

drawnow; % DON'T MOVE THIS LINE: is a workaround to geck 268506

hOverviewPan = imoverviewpanel(hYourTubes,hScrollPanelIm);
set(hOverviewPan,'Units','Normalized','Position',[0 .5 1 .5]);

set(hYourTubes,'Colormap',get(hScrollPanelFig,'Colormap'));
set(hYourTubes,'Renderer',get(hScrollPanelFig,'Renderer'));

linkFig = linkprop([hScrollPanelFig hYourTubes],'Colormap');
setappdata(hYourTubes, 'OverviewListeners', linkFig);

iptwindowalign(hYourTubes, 'right',hScrollPanelFig, 'left');
iptwindowalign(hYourTubes, 'top',hScrollPanelFig, 'top');

set(hOriginalFig,'Visible','on');

updateZoomButtons(apiScrollPanel.getMagnification())
magCallbackID = apiScrollPanel.addNewMagnificationCallback(@updateZoomButtons);

set(hYourTubes, 'DeleteFcn', @(src, varargin) apiScrollPanel.removeNewMagnificationCallback(magCallbackID));
reactToImageChangesInFig(hScrollPanelIm,@deleteFcn);

handles.hScrollPanelIm = hScrollPanelIm;
handles.hScrollPanelFig = hScrollPanelFig;
handles.hOverviewPan = hOverviewPan;
guidata(hYourTubes,handles);

set(handles.roi,'Enable','on');


% --------------------------------------------------------------------
function roi(varargin)
global DEBUG;
if (DEBUG), fprintf('roi\n'); end

hObject = findobj('Tag','YourTubes');
handles = guidata(hObject);
if ishandle(handles.roi)
    hroi = handles.roi;
    set(hroi,'Enable','off')
end
try % This exception must be changed
    [x,y,ROI,xi,yi,roix,roiy] = poly(handles.hScrollPanelFig);
    setappdata(hObject,'ROI',ROI);
    pFig.x = x;
    pFig.y = y;
    pFig.xi = xi;
    pFig.yi = yi;
    pFig.roix = roix;
    pFig.roiy = roiy;
    setappdata(hObject,'polyfig',pFig);
catch
    close all;
    eid = sprintf('roi:wrongPolyAcquisition');
    msg = 'Region of interest must be a parallelogram.';
    error(eid,'%s',msg);
end

if ishandle(handles.roi)
    hroi = handles.roi;
    set(hroi,'Enable','on')
end

IMG = getappdata(hObject,'IMG');
ROI = getappdata(hObject,'ROI');

wb = waitbar(0,'Finding region of interest. Please wait...');
IMG_ROI = uint8(255*ones(size(IMG,1),size(IMG,2),3));
for i=1:size(IMG,1)
    for k=1:size(IMG,2)
        if(ROI(i,k)==1)
            IMG_ROI(i,k,1:3) = IMG(i,k,1:3);
        end
    end
    waitbar(i/size(IMG,1));
end
close(wb);
wa = getWorkArea();
w = wa.width;
h = wa.height;
IMG_RES = IMG;
IMG_ROI_RES = IMG_ROI;
while (w<=size(IMG_RES,1) || h<=size(IMG_RES,2))
    IMG_RES = imresize(IMG_RES, 0.25);
    IMG_ROI_RES = imresize(IMG_ROI_RES, 0.25);
end

hRoiFig = createFigure('Region of interest');
set(hRoiFig,'Visible','off');
imshow(IMG_ROI_RES,'InitialMagnification',100);
hold on;
h = imshow(IMG_RES,'InitialMagnification',100);
set(h, 'AlphaData', 0.4);

iptwindowalign(hObject, 'right',hRoiFig, 'left');
iptwindowalign(hObject, 'top',hRoiFig, 'top');

transformImage;

set(hRoiFig,'Visible','on');

button = questdlg('Does this ROI satisfy your choice?','Confirm ROI','Yes','No','Yes');
switch button
    case 'Yes',
        set(handles.slider1, 'Enable', 'on');
        set(handles.slider2, 'Enable', 'on');
        set(handles.slider3, 'Enable', 'on');
        set(handles.edit11, 'Enable', 'on');
        set(handles.edit22, 'Enable', 'on');
        set(handles.edit33, 'Enable', 'on');
        set(handles.pushbutton, 'Enable', 'on');
        detectPipes;
    case 'No',
        delete(hRoiFig);
end

%set(handles.measurePipes,'Enable','on');


% --------------------------------------------------------------------
function transformImage(varargin)
global DEBUG;
if (DEBUG), fprintf('transformImage\n'); end

hObject = findobj('Tag','YourTubes');
IMG = getappdata(hObject,'IMG'); % Image IMG
% ROI = getappdata(hObject,'ROI'); % Image ROI
pFig = getappdata(hObject,'polyfig');
tFig = getappdata(hObject,'tubesfig');

% The points acquired using ROI are associated to the refrement system
roix = pFig.roix; % abscissa of parallelogram [pxl]
roiy = pFig.roiy; % ordinate of parallelogram [pxl]
% xi = pFig.xi;   % abscissa of parallelogram [axes]
% yi = pFig.yi;   % ordinate of parallelogram [axes]
[Xi,Yi,Xf,Yf] = associatePoints(roix,roiy);

in.Xi = Xi; in.Yi = Yi; in.Xf = Xf; in.Yf = Yf;
setappdata(hObject,'indet',in);

% Removing unuseful details (to reduce time computation)
ymax = max(Yi); xmax = max(Xi);
ymin = min(Yi); xmin = min(Xi);
rect = [xmin ymin (xmax-xmin) (ymax-ymin)];
IMG_CROP = imcrop(IMG, rect);

Xi = Xi - xmin; % Translating abscissa
Yi = Yi - ymin; % Translating ordinate

basePoints = [Xf Yf];
inputPoints = [Xi Yi];
T = cp2tform(inputPoints,basePoints,'projective');
[IMG_CROP_HOM, Xoffset, Yoffset] = imtransform(IMG_CROP,T);

rect = [Xf(1)-Xoffset(1) Yf(1)-Yoffset(1) (Xf(4) - Xf(1)) (Yf(2) - Yf(1))];
IMG_CROP_HOM_CROP = imcrop(IMG_CROP_HOM, rect);
tFig.IMG = IMG_CROP_HOM_CROP;

if DEBUG,
    save IMG_CROP_HOM_CROP IMG_CROP_HOM_CROP;
    createFigure('Image croped'); imshow(IMG_CROP);
    createFigure('Image croped 2-D spatial transformation'); imshow(IMG_CROP_HOM);
    createFigure('Image croped 2-D spatial transformation croped'); imshow(IMG_CROP_HOM_CROP);
end

setappdata(hObject,'tubesfig',tFig);


% --------------------------------------------------------------------
function [Xi,Yi,Xf,Yf] = associatePoints(roix,roiy)
global DEBUG;
if (DEBUG), fprintf('associatePoints\n'); end

Xf = [0; 0; 1000; 1000];
Yf = [0; 300; 300; 0];

xy = [roix roiy];
xySort = sortrows(xy,1);
xySortDiv1 = xySort(1:2,:);
xySortDiv2 = xySort(3:4,:);
xySortDiv1Sort = sortrows(xySortDiv1,2);
xySortDiv2Sort = sortrows(xySortDiv2,2);
x1 = xySortDiv1Sort(1,1);
y1 = xySortDiv1Sort(1,2);

x2 = xySortDiv1Sort(2,1);
y2 = xySortDiv1Sort(2,2);

x3 = xySortDiv2Sort(2,1);
y3 = xySortDiv2Sort(2,2);

x4 = xySortDiv2Sort(1,1);
y4 = xySortDiv2Sort(1,2);

Xi = [x1; x2; x3; x4];
Yi = [y1; y2; y3; y4];


% --------------------------------------------------------------------
function detectPipes(varargin)
global DEBUG;
if (DEBUG), fprintf('detectPipes\n'); end

hObject = findobj('Tag','YourTubes');
handles = guidata(hObject);
tFig = getappdata(hObject,'tubesfig');
IMG_CROP_HOM_CROP = tFig.IMG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% BEGIN ALGHORITM FOR TUBES DETECTION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ibw = rgb2gray(IMG_CROP_HOM_CROP);
Ibw = Ibw(1:100,:);
Idouble = im2double(Ibw);

% % Executing SVD compression
% [U,Sigma,V] = svd(Idouble);
% C = 10;
% Idouble = U(:,1:C)*Sigma(1:C,1:C)*V(:,1:C)';

l = size(Ibw);

Isum = sum(Ibw);
%Ivar = var(Idouble).*100;
Istd = std(Idouble).*100;
%Media = mean(Istd);
Igradsum = gradient(Isum);
%Igradvar = gradient(Ivar).*10;
%MediaGrad = mean(Igradsum);
Igradstd = gradient(Istd).*10;
%Ilapsum = gradient(Igradsum);
%Ilapvar = gradient(Igradvar);

SogliaStd = handles.th1;
SogliaGradSumInf = handles.th2;
SogliaGradSum = handles.th3;

tokenVar = ones(1,l(2)).*100;
for k2 = 1:l(2)
    if Istd(k2) < SogliaStd
        tokenVar(k2) = 0;
    end
end

if DEBUG, stairs(tokenVar+20,'y'); end

TubeNum = 0;
TubeDiam = [];
TubePos = [];
InTube = false;
wb = waitbar(0,'Tubes detection. Please wait...');
stop = l(2)-10;
start = 11;
for k2 = start:stop   
    if tokenVar(k2) == 0 && ~InTube
        InTube = true;
        TubeNum = TubeNum + 1;
        TubeDiam(TubeNum)= 0;
        TubePos(TubeNum) = k2;
    end
    if tokenVar(k2) == 0 && InTube
        if ((-SogliaGradSumInf) > Igradsum(k2) || Igradsum(k2) > SogliaGradSum) && max(abs(Igradsum((k2-2):(k2+5))))==abs(Igradsum(k2)) && TubeNum > 0 && TubeDiam(TubeNum) > 5
            TubeNum = TubeNum + 1;
            TubeDiam(TubeNum)= 0;
            TubePos(TubeNum) = k2;
        else
            TubeDiam(TubeNum) = TubeDiam(TubeNum) + 1;
        end
    end
    if tokenVar(k2) == 100
        InTube = false;
        if TubeNum > 1 && TubeDiam(TubeNum) < 10
            TubeNum = TubeNum - 1;
        end
    end
    waitbar(k2/(abs(stop-start)));
end
close(wb);
     
Pipes = zeros(TubeNum,3);
for k3 = 1:TubeNum
    diam = TubeDiam(k3)/l(2)*100;
    posy = diam/2;
    posx = TubePos(k3)/l(2)*100 + diam/2;
    Pipes(k3,:) = [diam posx posy];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% END ALGHORITM %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setappdata(hObject,'numberPipes',TubeNum);
setappdata(hObject,'tubesfig',tFig);

[IMG_MEA imPos] = createMeasureFig(Pipes);  %#ok<NASGU>


% --------------------------------------------------------------------
function [IMG_MEA mfPos] = createMeasureFig(Pipes)
global DEBUG;
if (DEBUG), fprintf('createMeasureFig\n'); end

UseAxes = 1; 

hObject = findobj('Tag','YourTubes');
handles = guidata(hObject);

%if(isempty(findobj('Tag','MeasureFig')))
    hMeasureFig = createFigure('Measurement Image');
    set(hMeasureFig,'Tag','MeasureFig');
    set(hMeasureFig,'Resize','off');
    set(hMeasureFig,'Visible','off');
    set(hMeasureFig,'Position',getInitialPosition);
%else
%    hMeasureFig = findobj('Tag','MeasureFig');
%end

handles.hMeasureFig = hMeasureFig;
guidata(hObject,handles);

[hp,huit] = createMeasureList(hMeasureFig,Pipes); %#ok<NASGU>
hpPos = get(hp,'Position'); %[pxl]

mFig = getappdata(hObject,'measurefig');
tFig = getappdata(hObject,'tubesfig');
numberPipes = getappdata(hObject,'numberPipes');

hmfPos = get(hMeasureFig,'Position'); %[pxl] getpixelposition(hMeasureFig);
%FigWidth = hmfPos(3); FigHeight = hmfPos(4);

space = 0.03*hmfPos(3); %[pxl] 

mFig.width = round(hmfPos(3) - 2*space); %[pxl]

maxD = max(Pipes(:,1)); % [units]
maxD = round(maxD*mFig.width/100); % [pxl]
gap = 15; % lasco [pxl]
mFig.height = round(maxD + gap); % [pxl]

IMG_MEA = 255*ones(mFig.height,mFig.width,3);
for i=[1,2,size(IMG_MEA,1)-1,size(IMG_MEA,1)]
    for k=1:size(IMG_MEA,2)
        IMG_MEA(i,k,1:3)=0;
    end
end
for k=[1,2,size(IMG_MEA,2)-1,size(IMG_MEA,2)]
    for i=1:size(IMG_MEA,1)
        IMG_MEA(i,k,1:3)=0;
    end
end
mFig.IMG = IMG_MEA;

x_mf = space;
y_mf = (hpPos(2) + hpPos(4)) + space;
w_mf = size(IMG_MEA,2);
h_mf = size(IMG_MEA,1);
mfPos = [x_mf, y_mf, w_mf, h_mf];

for n=1:numberPipes
    dEllipse = round(Pipes(n,1)*w_mf/100);
    xEllipse = 1 + round(x_mf + Pipes(n,2)*w_mf/100 - dEllipse/2);
    yEllipse = 1 + round(y_mf + Pipes(n,3)*w_mf/100 - dEllipse/2);

    annotation(hMeasureFig,'ellipse','LineWidth',1.5,...
        'FaceColor',[0.7529 0.7529 0.7529],...
        'Unit','pixel',...
        'Position',[xEllipse yEllipse dEllipse dEllipse],...
        'Color',[rand(1) rand(1) rand(1)]);
    % Create textbox
    annotation(hMeasureFig,'textbox','String',{num2str(n)},'FitHeightToText','on',...
    'LineStyle','none',...
    'Unit','pixel',...
    'HorizontalAlignment','center',... 
    'Color',[128 8 8]./255,...
    'Position',[round(xEllipse+dEllipse/2-1), round(yEllipse+dEllipse+20), 0,0]);
end

tFig.width = mFig.width; %[pxl]
tFig.height = mFig.height; %[pxl]
IMG_TUB = tFig.IMG;
IMG_TUB_RES = imresize(IMG_TUB,[tFig.height tFig.width]);
for i=[1,2,size(IMG_TUB_RES,1)-1,size(IMG_TUB_RES,1)]
    for k=1:size(IMG_TUB_RES,2)
        IMG_TUB_RES(i,k,1:2)=145;
        IMG_TUB_RES(i,k,3)=255;
    end
end
for k=[1,2,size(IMG_TUB_RES,2)-1,size(IMG_TUB_RES,2)]
    for i=1:size(IMG_TUB_RES,1)
        IMG_TUB_RES(i,k,1:2)=145;
        IMG_TUB_RES(i,k,3)=255;
    end
end

x_tf = space;
y_tf = y_mf + h_mf + space;
w_tf = size(IMG_TUB_RES,2);
h_tf = size(IMG_TUB_RES,1);
tfPos = [x_tf, y_tf, w_tf, h_tf];

% xPipes = [];
% yPipes = [];
% dPipes = [];
% for n=1:numberPipes
%     xPipes = x+round(Pipes(n,2)*w/100); %[pxl]
%     yPipes = y+round(Pipes(n,3)*w/100); %[pxl]
%     dPipes = round(Pipes(n,1)*w/100);   %[pxl]
% end

set(hMeasureFig,'Position',[hmfPos(1), hmfPos(2), hmfPos(3), h_mf+h_tf+4*space+hpPos(4)]);
XLim = [0 100];
YLim = [0 round(mFig.height*100/mFig.width)];
imAxes = axes('Units','pixels','Position',mfPos,'YLim',YLim,'XLim',XLim);
box('on');
if(~UseAxes)
    imshow(IMG_MEA,'Parent',imAxes);
end
tfAxes = axes('Units','pixels','Position',tfPos);
imshow(IMG_TUB_RES,'Parent',tfAxes);

setappdata(hObject,'measurefig',mFig);
setappdata(hObject,'tubesfig',tFig);

iptwindowalign(hObject, 'vcenter', hMeasureFig, 'vcenter');
set(hMeasureFig,'Visible','on');


% --------------------------------------------------------------------
function [hp,huit] = createMeasureList(hMeasureFig,Pipes)
global DEBUG;
if (DEBUG), fprintf('createMeasureList\n'); end

rows = size(Pipes,1);
cols = size(Pipes,2);

data = cell(rows, cols);
for r = 1:rows
    for c = 1:cols
        data{r,c} = sprintf('%g',Pipes(r,c));
    end
end

[hp,huit] = table(hMeasureFig,data,[]);
set(hp,'BorderType','beveledin');


%-------------------------------------------------------------------
function deleteFcn(varargin)
global DEBUG;
if (DEBUG), fprintf('deleteFcn\n'); end

hObject = findobj('Tag','YourTubes');
handles = guidata(hObject);
if ishandle(handles.hOverviewPan), delete(handles.hOverviewPan); end
if ishandle(handles.filePrint), set(handles.filePrint,'Enable','off'); end
if ishandle(handles.zoomInButton), set(handles.zoomInButton,'Enable','off'); end
if ishandle(handles.zoomOutButton), set(handles.zoomOutButton,'Enable','off'); end
if ishandle(handles.roi), set(handles.roi,'Enable','off'); end
%if ishandle(handles.measurePipes), set(handles.measurePipes,'Enable','off'); end


% --------------------------------------------------------------------
function updateZoomButtons(mag)
global DEBUG;
if (DEBUG), fprintf('updateZoomButtons\n'); end

hObject = findobj('Tag','YourTubes');
handles = guidata(hObject);
apiSP = getappdata(hObject,'apiSP');

if ishandle(hObject)
    if mag <= apiSP.getMinMag();
        set(handles.zoomOutButton,'Enable','off');
    else
        set(handles.zoomOutButton,'Enable','on');
    end

    if mag>=1024 % arbritrary big choice
        set(handles.zoomInButton,'Enable','off');
    else
        set(handles.zoomInButton,'Enable','on');
    end
end


% --------------------------------------------------------------------
function createToolbar(hObject,handles)
global DEBUG;
if (DEBUG), fprintf('createToolbar\n'); end

toolbar =  uitoolbar(hObject);
[iconRoot,iconRootMATLAB] = ipticondir; %#ok<NASGU>

zoomInIcon = makeToolbarIconFromPNG(fullfile(iconRoot,'overview_zoom_in.png'));
handles.zoomInButton = createToolbarPushItem(toolbar,zoomInIcon,{@zoomIn},'Zoom in');

zoomOutIcon = makeToolbarIconFromPNG(fullfile(iconRoot,'overview_zoom_out.png'));
handles.zoomOutButton = createToolbarPushItem(toolbar,zoomOutIcon,{@zoomOut},'Zoom out');

roiIcon = makeToolbarIconFromPNG(fullfile(iconRoot,'point.png'));
handles.roi = createToolbarPushItem(toolbar,roiIcon,{@roi},'ROI acquisition');

%measurePipesIcon = makeToolbarIconFromGIF(fullfile(iconRoot,'distance_tool.gif'));
%handles.measurePipes = createToolbarPushItem(toolbar,measurePipesIcon,{@detectPipes},'Pipes measurment');

guidata(hObject,handles);

set(handles.zoomInButton,'Enable','off');
set(handles.zoomOutButton,'Enable','off');
set(handles.roi,'Enable','off');
%set(handles.measurePipes,'Enable','off');


%-------------------------------------------------------------------
function item = createToolbarPushItem(toolbar,icon,callback,tooltip)
global DEBUG;
if (DEBUG), fprintf('createToolbarPushItem\n'); end

item = uipushtool(toolbar,...
    'Cdata',icon,...
    'TooltipString',tooltip,...
    'Tag',lower(tooltip),...
    'ClickedCallback',callback);


% --------------------------------------------------------------------
function zoomIn(varargin)
global DEBUG;
if (DEBUG), fprintf('zoomIn\n'); end

hObject = findobj('Tag','YourTubes');
apiSP = getappdata(hObject,'apiSP');
newMag = findZoomMag('in',apiSP.getMagnification());
apiSP.setMagnification(newMag);


% --------------------------------------------------------------------
function zoomOut(varargin)
global DEBUG;
if (DEBUG), fprintf('zoomOut\n'); end

hObject = findobj('Tag','YourTubes');
apiSP = getappdata(hObject,'apiSP');
newMag = findZoomMag('out',apiSP.getMagnification());
apiSP.setMagnification(newMag);


% --------------------------------------------------------------------
function hFig = createFigure(figName)
global DEBUG;
if (DEBUG), fprintf('createFigure\n'); end

if(ischar(figName))
hFig = figure('Toolbar','none',...
    'Menubar','none',...
    'IntegerHandle','on',... %'IntegerHandle','off',...
    'NumberTitle','on',...
    'Tag','YourTubesFig',...
    'Name',figName,...
    'Visible','on',...
    'HandleVisibility','on');%'HandleVisibility','callback',...
else
    eid = sprintf('createFigure:%s:invalidFigName',figName);
    msg = 'Figure name must be a string.';
    error(eid,'%s',msg);
end


% --------------------------------------------------------------------
function initPos = getInitialPosition
global DEBUG;
if (DEBUG), fprintf('getInitialPosition\n'); end

wa =getWorkArea();
SCALE1 = 0.7;
SCALE2 = 0.7;
w = SCALE1*wa.width;
h = SCALE2*wa.height;
x = wa.left + (wa.width - w)/2;
y = wa.bottom + (wa.height - h)/2;
initPos = round([x y w h]);


% --------------------------------------------------------------------
function file_open_Callback(hObject, eventdata, handles)  %#ok<DEFNU,INUSL>
global DEBUG;
if (DEBUG), fprintf('file_open_Callback\n'); end

[fileName, userCanceled] = imgetfile;
if userCanceled
    return;
else
    file_close_all_Callback(hObject, eventdata, handles);
    set(handles.filePrint,'Enable','on');
    choiceImage(hObject, eventdata, handles, fileName);
end


% --------------------------------------------------------------------
function file_print_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
global DEBUG;
if (DEBUG), fprintf('file_print_Callback\n'); end

printImageToFigure(gcf);


% --------------------------------------------------------------------
function file_close_Callback(hObject, eventdata, handles)  %#ok<DEFNU,INUSL>
global DEBUG;
if (DEBUG), fprintf('file_close_Callback\n'); end

close(handles.YourTubes);
figs = findall(0,'Type','figure','Tag','YourTubesFig');
close(figs);


% --------------------------------------------------------------------
function file_close_all_Callback(hObject, eventdata, handles) %#ok<INUSL>
global DEBUG;
if (DEBUG), fprintf('file_close_all_Callback\n'); end

figs = findall(0,'Type','figure','Tag','YourTubesFig');
close(figs);
set(handles.slider1, 'Enable', 'off');
set(handles.slider2, 'Enable', 'off');
set(handles.slider3, 'Enable', 'off');
set(handles.edit11, 'Enable', 'off');
set(handles.edit22, 'Enable', 'off');
set(handles.edit33, 'Enable', 'off');
set(handles.pushbutton, 'Enable', 'off');
set(handles.pushbutton, 'Visible', 'off');
if ishandle(handles.filePrint), set(handles.filePrint,'Enable','off'); end
if ishandle(handles.zoomInButton), set(handles.zoomInButton,'Enable','off'); end
if ishandle(handles.zoomOutButton), set(handles.zoomOutButton,'Enable','off'); end
if ishandle(handles.roi), set(handles.roi,'Enable','off'); end
%if ishandle(handles.measurePipes), set(handles.measurePipes,'Enable','off'); end


% --------------------------------------------------------------------
function help_about_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
global DEBUG;
if (DEBUG), fprintf('help_about_Callback\n'); end
about;


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
global DEBUG;
if (DEBUG), fprintf('file_menu_Callback\n'); end


% --------------------------------------------------------------------
function help_menu_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
global DEBUG;
if (DEBUG), fprintf('help_menu_Callback\n'); end


% --------------------------------------------------------------------
function varargout = YourTubes_OutputFcn(hObject, eventdata, handles) %#ok<INUSL>
% --- Outputs from this function are returned to the command line.
global DEBUG;
if (DEBUG), fprintf('YourTubes_OutputFcn\n'); end

varargout{1} = handles.output; % Get default command line output from handles structure


% --------------------------------------------------------------------
function slider1_Callback(hObject, eventdata, handles)  %#ok<INUSL,DEFNU>
global DEBUG;
if (DEBUG), fprintf('slider1_Callback\n'); end

sliderValue = get(hObject,'Value');
sliderValueRound = round(sliderValue);
if(sliderValueRound > sliderValue)
    sliderValue = sliderValueRound;
else
    sliderValue = sliderValueRound+0.5;
end
set(handles.slider1, 'String', sliderValue);
set(handles.edit11, 'String', sliderValue);
handles.th1 = round(sliderValue/10)*10;
guidata(hObject,handles); 

detectPipes;

% --------------------------------------------------------------------
function slider1_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
global DEBUG;
if (DEBUG), fprintf('slider1_CreateFcn\n'); end

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function slider2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
global DEBUG;
if (DEBUG), fprintf('slider2_Callback\n'); end

sliderValue = get(hObject,'Value');
set(handles.slider2, 'String', round(sliderValue/10)*10);
set(handles.edit22, 'String', round(sliderValue/10)*10);
handles.th2 = round(sliderValue/10)*10;
guidata(hObject,handles); 

detectPipes;

% --------------------------------------------------------------------
function slider2_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
global DEBUG;
if (DEBUG), fprintf('slider2_CreateFcn\n'); end

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function slider3_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
global DEBUG;
if (DEBUG), fprintf('slider3_Callback\n'); end

sliderValue = get(hObject,'Value');
set(handles.slider3, 'String', round(sliderValue/10)*10);
set(handles.edit33, 'String', round(sliderValue/10)*10);
handles.th3 = round(sliderValue/10)*10;
guidata(hObject,handles); 
% get(hObject,'Min') and get(hObject,'Max') to determine range of slider

detectPipes;

% --------------------------------------------------------------------
function slider3_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
global DEBUG;
if (DEBUG), fprintf('slider3_CreateFcn\n'); end

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function edit11_Callback(hObject, eventdata, handles)
global DEBUG;
if (DEBUG), fprintf('edit11_Callback\n'); end

val = str2double(get(hObject,'String'));
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.slider1,'Min') && ...
   val <= get(handles.slider1,'Max')
   set(handles.slider1,'Value',val);
else
   handles.th1 = round(val);
   guidata(hObject,handles); 
end

detectPipes;

% --------------------------------------------------------------------
function edit11_CreateFcn(hObject, eventdata, handles)
global DEBUG;
if (DEBUG), fprintf('edit11_CreateFcn\n'); end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit22_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
global DEBUG;
if (DEBUG), fprintf('edit22_Callback\n'); end

val = str2double(get(hObject,'String'));
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.slider2,'Min') && ...
   val <= get(handles.slider2,'Max')
   set(handles.slider2,'Value',val);
else
   handles.th2 = round(val/10)*10;
   guidata(hObject,handles); 
end

detectPipes;

% --------------------------------------------------------------------
function edit22_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
global DEBUG;
if (DEBUG), fprintf('edit22_CreateFcn\n'); end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit33_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
global DEBUG;
if (DEBUG), fprintf('edit33_Callback\n'); end

val = str2double(get(hObject,'String'));
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.slider3,'Min') && ...
   val <= get(handles.slider3,'Max')
   set(handles.slider3,'Value',val);
else
   handles.th3 = round(val/10)*10;
   guidata(hObject,handles); 
end

detectPipes;

% --------------------------------------------------------------------
function edit33_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
global DEBUG;
if (DEBUG), fprintf('edit33_CreateFcn\n'); end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function pushbutton_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
global DEBUG;
if (DEBUG), fprintf('pushbutton_Callback\n'); end

detectPipes;
