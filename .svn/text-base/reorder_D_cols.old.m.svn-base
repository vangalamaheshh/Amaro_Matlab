function D=reorder_D_cols(D,varargin)

D=add_history(D,mfilename,varargin{:});

if length(varargin)==1
  idx=varargin{1};
else
  idx=find_D_cols(D,varargin{1},varargin{2});
end

try
  D.dat=D.dat(:,idx);
catch %memory problems
  D.dat(:,setdiff(1:size(D.dat,2),idx))=[];
end

if isfield(D,'sdesc')
  if iscell(D.sdesc)
    D.sdesc=D.sdesc(idx);
  else
    D.sdesc=deblank(D.sdesc(idx,:));
  end
end

if isfield(D,'affy_call')
  D.affy_call=D.affy_call(:,idx);
end
if isfield(D,'affy_calls')
  D.affy_calls=D.affy_calls(:,idx);
end
if isfield(D,'scans')
  D.scans=D.scans(idx);
end
if isfield(D,'supdat')
  D.supdat=D.supdat(:,idx);
end
if isfield(D,'sup')
  D.sup=D.sup(idx);
end
if isfield(D,'sscale') & ~isempty(D.sscale)
  D.sscale=D.sscale(idx,:);
end
if isfield(D,'sdat') D.sdat=D.sdat(:,idx); end
if isfield(D,'lsdat') D.lsdat=D.lsdat(:,idx); end
if isfield(D,'fdat') D.fdat=D.fdat(:,idx); end
if isfield(D,'scale2vec')
  D.scale2vec=D.scale2vec(idx);
end

vector_fields={'residx','origidx','gcm_name','cbs_rl','sis','used_normals','medians',...
               'peaks','joins','final_dist','sscaling'};
if isfield(D,'sample_fields')
  vector_fields=[ D.sample_fields];
end
for j=1:length(vector_fields)
  if isfield(D,vector_fields{j})
    x=getfield(D,vector_fields{j});
    x=x(idx);
    D=setfield(D,vector_fields{j},x);
  end
end

matrix_fields={'prectrls','smooth','cbs','cbs_fixed','sm1','sm2','sm2j','sm3','raw','si','flag','hmm','pre_control_sum','sitab','ref','level','scores'};
if isfield(D,'matrix_fields')
  matrix_fields=[ D.matrix_fields];
end
for j=1:length(matrix_fields)
  if isfield(D,matrix_fields{j})
    x=getfield(D,matrix_fields{j});
    x=x(:,idx);
    D=setfield(D,matrix_fields{j},x);
  end
end

