function [X,XX] = load_dRanger_input(fname,params)

if ~exist('params','var'), params =[]; end
% params = impose_default_value(params,'dedup',true);

fprintf('Loading %s... ',fname);tic
load(fname);
if ~exist('X','var'), error('%s does not contain "X"',fname); end

if isstruct(X)        % struct (old style)

  XX = X;

  X = [XX.namenumber XX.chr1 XX.strand1 XX.start1 XX.end1 XX.chr2 XX.strand2 ...
       XX.start2 XX.end2 XX.qual1 XX.qual2 XX.rgrp];

elseif isnumeric(X)   % matrix (new style)

  XX=[];
  XX.namenumber = X(:,1);
  XX.chr1 = X(:,2);
  XX.strand1 = X(:,3);
  XX.start1 = X(:,4);
  XX.end1 = X(:,5);
  XX.chr2 = X(:,6);
  XX.strand2 = X(:,7);
  XX.start2 = X(:,8);
  XX.end2 = X(:,9);
  XX.qual1 = X(:,10);
  XX.qual2 = X(:,11);
  XX.rgrp = X(:,12);

else
  error('Unknown format for "X"');
end

toc
