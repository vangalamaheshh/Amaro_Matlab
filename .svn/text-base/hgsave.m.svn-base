function hgsave(varargin)
% HGSAVE  Saves an HG object hierarchy to a MAT file.
% 
% HGSAVE(H, 'Filename') saves the objects identified by handle array H
% to a file named 'Filename'.  If 'Filename' contains no extension,
% then the extension '.fig' is added.  If H is a vector, none of the
% handles in H may be ancestors or descendents of any other handles in
% H.
%
% HGSAVE('Filename') saves the current figure to a file named
% 'Filename'.
%
% HGSAVE(..., 'all') overrides the default behavior of excluding
% non-serializable objects from the file.  Items such as the default
% toolbars and default menus are marked as non-serializable, and are
% normally not saved in FIG-files, because they are loaded from
% separate files at figure creation time.  This allows the size of
% FIG-files to be reduced, and allows for revisioning of the default
% menus and toolbars without affecting existing FIG-files. Passing
% 'all' to HGSAVE insures non-serializable objects are also saved.
% Note: the default behavior of HGLOAD is to ignore non-serializable
% objects in the file at load time, and that can be overridden using
% the 'all' flag with HGLOAD.
%
% HGSAVE(..., '-v6') saves a FIG-file that can be loaded by versions
% prior to MATLAB 7. When creating a figure to be saved and used in a
% version prior to MATLAB 7 use the 'v6' option to the plotting
% commands. See the help for PLOT, BAR and other plotting commands for
% more information.
%
% See also HGLOAD, SAVE.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%   Copyright 1984-2007 The MathWorks, Inc.
%   D. Foti  11/10/97

error(nargchk(1, inf, nargin,'struct'));

% pull off handle + 'filename,' or just 'filename'
if ischar(varargin{1})
  h = gcf;
  filename = varargin{1};
  first_pass_through = 2;
else
  h = varargin{1};
  if any(~ishandle(h))
    error('MATLAB:hgsave:InvalidHandle','Invalid handle');
  end
  filename = varargin{2};
  first_pass_through = 3;
end

% any trailing args that start with '-' are passed to save
save_args = nargin+1;
while save_args > first_pass_through && ...
      ischar(varargin{save_args-1}) && ...
      varargin{save_args-1}(1) == '-'
  save_args = save_args - 1;
end

[path, file, ext] = fileparts(filename);

% fileparts returns everything from the last . to the end of the
% string as the extension so the following test will catch
% an extension with 0, 1, or infinity dots.
% for example, all these filenames will have .fig added to the end:
%  foo.
%  foo..
%  foo.bar.
%  foo.bar...

if isempty(ext) || strcmp(ext, '.')
  filename = fullfile(path, [file , ext, '.fig']);
end

% Revision encoded as 2 digits for major revision,
% 2 digits for minor revision, and 2 digits for
% patch revision.  This is the minimum revision 
% required to fully support the file format.
% e.g. 070000 means 7.0.0

% if saving a figure and plotedit, zoom, camera toolbar,rotate3d or
% brushing are on, save their states and
% turn them off before saving
% and if scribe clear mode callback appdata
% exists, remove it.
plotediton = zeros(length(h),1);
rotate3dstate = cell(length(h),1);
zoomstate = cell(length(h),1);
datacursorstate = false(length(h),1);
panstate = false(length(h),1);
scmcb = cell(length(h),1);
camtoolbarstate = zeros(length(h),1);
camtoolbarmode = cell(length(h),1);
% brushing = false(length(h),1);

olddata = cell(length(h),1);
for i = 1:length(h)
  if strcmp(get(h(i), 'type'), 'figure')
    camtoolbarstate(i) = cameratoolbar(h(i),'GetVisible');
    plotediton(i) = plotedit(h(i), 'isactive');
    rotate3dstate{i} = getappdata(h(i),'Rotate3dOnState');
    zoomstate{i} = getappdata(h(i),'ZoomOnState');
    datacursorstate(i) = strcmp(datacursormode(h(i),'ison'),'on');
    panstate(i) = pan(h(i),'ison');
%     brushing(i) = brush(h(i),'ison');
    s = getappdata(h(i),'ScribeClearModeCallback');
    if camtoolbarstate(i)
        camtoolbarmode = cameratoolbar(h(i),'GetMode');        
        cameratoolbar(h(i),'save');
    end
    if ~isempty(s) && iscell(s)
        scmcb{i} = s;
        rmappdata(h(i),'ScribeClearModeCallback');
    end
    if plotediton(i)
        plotedit(h(i),'off');
    end
    if ~isempty(rotate3dstate{i})
        rotate3d(h(i),'off');
    end
    if ~isempty(zoomstate{i})
        zoom(h(i),'off');
    end    
    if datacursorstate(i)
        datacursormode(h(i),'off');
    end
    % Serialize data tip information:
    hDCM = datacursormode(h(i));
    hDCM.serializeDatatips;
    if panstate(i)
        pan(h(i),'off');
    end
%     if brushing(i)
%         brush(h(i),'off')
%     end
  end
  fch = findall(h(i));
  for k=1:length(fch)
    hh = handle(fch(k));
    if ismethod(hh,'preserialize')
      olddata{i}(2*k-1:2*k) = {hh,preserialize(hh)};
    end
  end
end

%If axes are linked, we ned to capture this in the saved file
allAxes = unique(findall(h,'Type','axes'));
l = length(allAxes);
linkage = [];
for i = 1:l
    %For all the axes which are linked to other axes, obtain the handle to the
    %linkprop objects
    if isappdata(allAxes(i),'graphics_linkaxes')
        temp_link = getappdata(allAxes(i),'graphics_linkaxes');
        if ishandle(temp_link), linkage = [linkage temp_link];end
    end
end
linkage = unique(linkage);
targets = [];
for i = 1:length(linkage)
    param = '';
    t = get(linkage(i),'Targets');
    props = get(linkage(i),'PropertyNames');
    if any(strcmp(props,'XLim'))
        param = strcat(param,'x');
    end
    if any(strcmp(props,'YLim'))
        param = strcat(param,'y');
    end
    for j = 1:length(t)
        %Only store this information if the target is being saved
        if any(allAxes == t(j))
            setappdata(t(j),'graphics_linkaxes_targets',i);
            setappdata(t(j),'graphics_linkaxes_props',param);
        else
            t(j) = handle(-500);
        end
    end
    t(~ishandle(t)) = [];    
    targets = [targets t];
end

% Serialize any attached behavior objects
% This is a work around until HG supports MCOS serialization
hWithBehaviors = localSerializeBehaviorObjects(h);

% Serialize the Annotation property.
% This is a work around until HG supports composite object / UDD
% serializeation.
localSerializeAnnotations(h);

%If we fail here, remove the additional application data
try
    hgS_070000 = handle2struct(h, varargin{first_pass_through:save_args-1}); %#ok
catch
    for i = 1:length(targets)
        rmappdata(targets(i),'graphics_linkaxes_targets');
        rmappdata(targets(i),'graphics_linkaxes_props');
    end
    rethrow(lasterror);
end

% Remove temporary behavior serialization data
localClearBehaviorSerialization(hWithBehaviors);
% Remove temporary annotaion serialization data.
localClearAnnotationSerialization(h);

% restore plotedit, zoom, camera toolbar and rotate3d states if saving a figure
for i = 1:length(h)
  oldobj = olddata{i};
  for k=1:2:length(oldobj)
    if ismethod(oldobj{k},'postserialize')
      postserialize(oldobj{k:k+1});
    end
  end
  if strcmp(get(h(i),'type'),'figure')
    % if the camera toolbar was on, restore it
    if camtoolbarstate(i)
        cameratoolbar(h(i),'toggle');
        cameratoolbar(h(i),'SetMode',camtoolbarmode);
    end
    % if plotedit was on, restore it
    if plotediton(i)
      plotedit(h(i),'on');
    end
    % if rotate3d was on, restore it
    if ~isempty(rotate3dstate{i})
      rotate3d(h(i),rotate3dstate{i});
    end
    % if zoom was on, restore it
    if ~isempty(zoomstate{i})
      zoom(h(i),zoomstate{i});
    end
    if datacursorstate(i)
      datacursormode(h(i),'on');
    end
    % Remove any appdata that was created by the serialization
    hDCM = datacursormode(h(i));
    hDCM.clearDatatipSerialization;
    if panstate(i)
      pan(h(i),'on');
    end
%     if brushing(i)
%       brush(h(i),'on');
%     end
    % if there was a scribeclearmodecallback, reset it
    if ~isempty(scmcb{i})
      setappdata(h(i),'ScribeClearModeCallback',scmcb{i});
    end   
  end
end

%Clean up
for i = 1:length(targets)
    rmappdata(targets(i),'graphics_linkaxes_targets');
    rmappdata(targets(i),'graphics_linkaxes_props');
end


save(filename, 'hgS_070000',varargin{save_args:end});

%-------------------------------------------------%
function ret = localSerializeBehaviorObjects(h)

% Find all the objects with non-empty behavior objects
% struct = empty struct, the default value of 'behavior' property
ret = findall(h,'-and','-not',{'Behavior',struct},'-function',@localDoSerialize);
 
%-------------------------------------------------%
function localSerializeAnnotations(h)

% Find all the objects that have an annotation property:
ret = findall(h,'-function',@(h)(isprop(handle(h),'Annotation')));
% If the property has an a serialized annotation, skip it. Otherwise,
% serialize its annotations:
for i = 1:numel(ret)
    if isappdata(ret(i),'SerializedAnnotationV7')
        continue;
    end
    hA = get(ret(i),'Annotation');
    if ~isa(hA,'hg.Annotation')
        continue;
    end
    % The convention for the "Annotation" property is that each property is
    % a handle. The handle contains the state we are interested in:
    hP = get(hA);
    serProp = structfun(@localhandle2struct,hP,'UniformOutput',false);
    setappdata(ret(i),'SerializedAnnotationV7',serProp);
end


%-------------------------------------------------%
function localClearBehaviorSerialization(h)
% Remove temporary appdata serialization 

for n = 1:length(h)
   if ishandle(h(n)) && isappdata(h(n),'SerializedBehaviorV7')
      rmappdata(h(n),'SerializedBehaviorV7');
   end
end

%-------------------------------------------------%
function localClearAnnotationSerialization(h)
% Remove temporary appdata serialization 

for n = 1:length(h)
   if ishandle(h(n)) && isappdata(h(n),'SerializedAnnotationV7')
      rmappdata(h(n),'SerializedAnnotationV7');
   end
end

%-------------------------------------------------%
function [ret]= localDoSerialize(h)
% For the supplies handle, find all the behavior objects that 
% support serialization. Serialize each object in appdata

ret = false;
b = get(h,'Behavior');
if ~isempty(b)
    
    % Find behavior objects with 'Serialize' = true
    b = struct2cell(b);
    b = [b{:}];
    b = find(b,'Serialize',true);
   
    if ~isempty(b)
        appdata = struct;
   
        % Loop through each behavior object and create a structure 
        % that represents the state
        count = 1;
        for n = 1:length(b)
            
            % Serialize behavior object as a structure in appdata
            try %#ok
                s = localhandle2struct(b(n));
                if ~isempty(s) 
                   appdata(count).class = class(b(n));
                   appdata(count).properties = s;
                   count = count + 1;
                end
            end
        end
        setappdata(double(h),'SerializedBehaviorV7',appdata);
        ret = true;
    end
end

%-------------------------------------------------%
function s = localhandle2struct(hThis)
% Converts a generic UDD handle to a structure for serialization

s = [];
hCls = classhandle(hThis);
hProp = get(hCls,'Properties');

% Loop through properties
for n = 1:length(hProp)
    p = hProp(n);
    if ishandle(p)
        propname = get(p,'Name');
        propval = get(hThis,propname);

        % Serialize any properties that are public set, non-default
        if isequal(p.AccessFlags.Serialize,'on') && ...
           isequal(p.AccessFlags.PublicSet,'on') && ...
           ~isequal(get(p,'FactoryValue'),p)
              s.(propname) = propval;
        end
    end
end


