function M1=add_sample_info(M,SI,map_using,platform)
%ADD_SAMPLE_INFO Return data structure with sample information added.  
%
%   M1 = ADD_SAMPLE_INFO(M,SI,MAP_USING,PLATFORM) returns a data structure
%   M1 that includes sample information from the sample information file.
%   The only samples (columns) included in data structure M1 are those that
%   are represented in the input sample info structure SI.  ADD_SAMPLE_INFO
%   also updates the supplemental information fields .supacc, .supdesc, and
%   .supdat of the data structure M. MAP_USING is an optional input that
%   sets information matching by 'name' or 'array' (default is 'name').
%

%idx = which rows to use in SI
idx=1:length(SI);
% %%%% OLD CODE -- should platform be added as a .gsdesc field?
% if (0) % FIXME: old code
%   if isfield(SI,'good')
%     idx=strmatch('yes',{SI.good});
%   else
%     idx=1:length(SI);
%   end
% 
%   if exist('platform','var')
%     platform_idx=strmatch(platform,{SI.platform});
%     idx=intersect(idx,platform_idx);
%   end
% end

% match names
if ~exist('map_using','var')
  map_using='name';
end

switch map_using
 case 'name'
  [nms,ui,uj]=unique(strvcat(SI(idx).name),'rows');
  [Mtmp,m1,m2]=match_string_sets(M.sdesc,nms);
 case 'array'
  [nms,ui,uj]=unique(strvcat(SI(idx).array),'rows');
  [Mtmp,m1,m2]=match_string_sets(M.sdesc,nms);
 otherwise
  error('no such field');
end

SI2=SI(idx(ui(m2)));

M1=reorder_D_cols(M,m1);
clear M;



% removed this to keep names as array 
%if ~strcmp(map_using,'name')
%  M1.sdesc={SI2.name};
%end
verbose('Keeping names as in M',10);

M1.sis=SI2;

%  MOVE THIS TO HAPPEN RIGHT BEFORE BATCH CORRECTION (batch numbers should
%  be assigned after plates merge)
% 
% % add Batch (now the batch is in the sample info file)
% [ub,ui,uj]=unique(strvcat(SI2.batch),'rows');
% M1=add_D_sup(M1,'BATCH','Batch',uj','cols');


% add Normal label to supdat
% (will divide by these samples)
v=zeros(size(M1.dat,2),1);
if isfield(SI2,'tumor_normal')
    % TCGA sample info has tumor_normal field
    v(strcmpi('Normal',{SI2.tumor_normal}))=1;
elseif isfield(SI2,'ploidy')
    % GCM normals have a ploidy of 2
    v(strcmpi('2',{SI2.ploidy}))=1;
else
    warning('No normal field in sample info.');
end
if ~any(v)
    warning('Data contains no normals.');
end
M1=add_D_sup(M1,'N','Normal',v','cols');

% add Control label to supdat if type sample info column exists
v=zeros(size(M1.dat,2),1);
if isfield(SI2,'type')
    v(strcmpi('control',{SI2.type}))=1;
else
    verbose('Did not find a control field... assuming none',10);
end
M1=add_D_sup(M1,'CTRL','Control',v','cols');

% add Gender label to supdat if gender sample info column exists
if isfield(SI2,'gender')
    v=ones(size(M1.dat,2),1);
    v(strcmp('F',{SI2.gender}))=2;
    M1=add_D_sup(M1,'GENDER: 1-M/2-F','GENDER: 1-Male/2-Female',v','cols');
else
    verbose('Did not find a gender field... assuming empty',10);
end

% %<<<<<<<<<<THIS SHOULD BE MOVED TO BEFORE NORMALIZATION
% % add Core label
% v=ones(size(M1.dat,2),1);
% if isfield(SI2,'core')
%   [cores,ci,cj]=unique(strvcat(SI2(:).core),'rows');
% else
%   cores='GENERAL';
%   cj=v;
% end
% M1=add_D_sup(M1,make_sup_list('CORE',cores),make_sup_list('CORE',cores),cj','cols');
% %<<<<<<<<<<<  

% add Past_QC label
v=zeros(size(M1.dat,2),1);
if isfield(SI2,'past_qc')
    v(strcmp('yes',{SI2.past_qc})) = 1;
else
    verbose('Did not find a past_qc field... assuming empty',10);
end
M1=add_D_sup(M1,'FORCE','FORCE',v','cols');

% add Rep label
v=zeros(size(M1.dat,2),1);
if isfield(SI2,'rep')
    v(strcmp('yes',{SI2.rep})) = 1;
else
    verbose('Did not find a rep field... assuming empty',10);
end

M1=add_D_sup(M1,'REP','Replicate',v','cols');

% add SI2 index
M1=add_D_sup(M1,'SI2','SI2 index',(1:length(SI2)),'cols');
