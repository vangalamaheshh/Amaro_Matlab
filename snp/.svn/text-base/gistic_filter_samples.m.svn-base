function [M1,idx]=gistic_filter_samples(M,SI,select_type)
%GISTIC_FILTER_SAMPLES - Filter data structure for desired samples
%
%   [M1,IDX] = GISTIC_FILTER_SAMPLES(M,SI,SELECT_TYPE) returns the data
%   structure M1 filtered on SELECT_TYPE and the indices IDX 
%   of input data structure M used to form M1.  
%   The output data structure M1 is the intersection of "good" 
%   data (indicated by 'empty','yes', 'y', or no value in SI.good) and data 
%   matching the requested SELECT_TYPE.  SELECT_TYPE is either indices (of SI), 
%   or it is a cell array of types sample types matching values in SI.type.
%
%   Updates:
%
%       21 Sept 08: Modified select good samples; error-catch if no 'good
%       field'.  Added help and documentation.  -jdobson@broad.mit.edu

%% Select good samples

%%%%%%%%REPLACE WITH
% tmp={SI.good};
% tmp=lower(cellfun_any('replace_empty(x,''EMPTY'')',tmp));
% good=grep('^(empty|yes|y)$',tmp,1);

%%%%%%%%%THIS
if isfield(SI,'good')
    empidx = find(cellfun(@isempty,{SI.good}));
    [SI(empidx).good] = deal('EMPTY');
else
    [SI.good] = deal('EMPTY');
end

goodmatches= regexpi({SI.good},'(EMPTY)|(y)|(yes)');
good = find(cell2mat(goodmatches));  %gives indices of "good" data

%%%%%%%%%%%%%

%% Filter on select type

if isnumeric(select_type)
  idx=select_type;
else
  %----- select samples to use (from sample info)
  g=[];
  dirname='';
  if ischar(select_type)  %turn select type into a cell array
    clear tmp;
    tmp{1}=select_type;
    select_type=tmp;
  end
  for i=1:length(select_type)
    % give indices of samples matching select type
    g=union(g,strmatch(lower(select_type{i}),lower(strvcat(SI(:).type)),'exact')); 
    dirname=[dirname '_' select_type{i}];
  end
  %- choose also the ploidy=2 samples (for normalization)
  p2=findstrings(strvcat(SI(:).ploidy),'2')'; 
  if isempty(g)
    M1=[];
    idx=[];
  else
    % select all those that are g or p2
    idx=union(g,p2);
  end
end

idx=unique(intersect(idx,good));  %make sure all selected samples are good


%% Update SIs and M1

if isempty(idx)
  M1=[];
else
  SIs=SI(idx);
  M1=add_sample_info(M,SIs,'array');  %add_sample_info updates M1
  M1=rmfield_if_exists(M1,{'sis','supacc','supdat','supdesc','orig'});  
end




