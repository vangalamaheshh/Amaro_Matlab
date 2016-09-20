function [Cnew] = preproc_scaledata(Craw,dataislog,operation,autosomes_only)
%PREPROC_SCALEDATA Scale data by column mean, median, or something fancy
%
%   [CNEW] = PREPROC_SCALEDATA(CRAW,DATAISLOG,operation,MASK) scales each column
%   of sample data in struct CRAW.dat to the data column's XX(Craw.dat) value and
%   returns the scaled data in CNEW. XX is some function, can be 'median', 'mean', or a function handle.
%   DATAISLOG is used to specify whether
%   the input data is log transformed (1) or linear (0). Default for
%   DATAISLOG is 1.  CRAW can be a struct or a cell array of structures.  
%---
% $Id$
% $Date: 2011-11-22 11:26:53 -0500 (Tue, 22 Nov 2011) $
% $LastChangedBy: schum $
% $Rev$

if ~exist('dataislog','var')
    dataislog = 1;
end

if ~iscell(Craw)
    Cnew = {Craw};
    wascell = 0;
else
    Cnew = Craw;
    wascell = 1;
end
%! schum note: if using disk fields, Cnew is sharing them with Craw

if ischar(operation)
    if strcmpi('mean',operation)
        fh = @(x,idx) nanmean(x(idx,:),1);
    elseif strcmpi('median',operation)
        fh = @(x,idx) nanmedian(x(idx,:),1);
    else
        error('unknown operation');
    end
elseif isempty(operation)
    fh = @(x,idx) nanmean(x(idx,:),1);
else
    fh = operation;
end

if ~exist('autosomes_only','var') || isempty(autosomes_only)
    autosomes_only = false;
end


for k = 1:length(Cnew)

    verbose('Scaling data using method: ', char(operation));
    
    if autosomes_only
        mask = find(Cnew{k}.chrn < 23); %!!!HUMGEN specific to human
    else
        mask = 1:getsize(Cnew{k},'dat',1);
    end

    if dataislog
        Cnew{k} = itrfcn2(Cnew{k},'dat','dat',1,@(x,m) x - repmat(fh(x,m),size(x,1),1),mask);
    else
        Cnew{k} = itrfcn2(Cnew{k},'dat','dat',1,@(x,m) x ./ repmat(fh(x,m),size(x,1),1),mask);
    end
end

if ~wascell
    Cnew = Cnew{1};
end


