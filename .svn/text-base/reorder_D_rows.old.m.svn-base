function D=reorder_D_rows(D,varargin)


D=add_history(D,mfilename,varargin{:});

if length(varargin)==0
  error('missing argument');
elseif length(varargin)==1
  idx=varargin{1};
else
  idx=find_D_rows(D,varargin{1},varargin{2});
end

D.dat=D.dat(idx,:);
if isfield(D,'gdesc')
  if iscell(D.gdesc)
    D.gdesc=D.gdesc(idx);
  else
    D.gdesc=deblank(D.gdesc(idx,:));
  end
end

if isfield(D,'gsymb')
  if iscell(D.gsymb)
    D.gsymb=D.gsymb(idx);
  else
    D.gsymb=deblank(D.gsymb(idx,:));
  end
end

if isfield(D,'gacc')
  if iscell(D.gacc)
    D.gacc=D.gacc(idx);
  else
    D.gacc=deblank(D.gacc(idx,:));
  end
end

if isfield(D,'id')
  if iscell(D.id)
    D.id=D.id(idx);
  else
    D.id=deblank(D.id(idx,:));
  end
end

if isfield(D,'symb')
  if iscell(D.symb)
    D.symb=D.symb(idx);
  else
    D.symb=deblank(D.symb(idx,:));
  end
end

if isfield(D,'affy_call')
  D.affy_call=D.affy_call(idx,:);
end

if isfield(D,'adat')
  D.adat=D.adat(idx,:,:);
end

if isfield(D,'gsupdat')
  if size(D.gsupdat,1)==size(D.dat,1) % rows match
    warning('recommend using a column per row in dat');
    D.gsupdat=D.gsupdat(idx,:);
  else
    D.gsupdat=D.gsupdat(:,idx);    
  end
end

vector_fields={'cyto','gorigidx','marker','chr','chrn','cM','pos','score','grg','start','end','ll','cyton','cyto_stain', ...
               'armn','chrloc','refgene','refgene_idx','grange','chrarmn','n_collapse'};
if isfield(D,'gene_fields')
  vector_fields=[ D.gene_fields];
end
for j=1:length(vector_fields)
  if isfield(D,vector_fields{j})
    x=getfield(D,vector_fields{j});
    x=x(idx);
    D=setfield(D,vector_fields{j},x);
  end
end

matrix_fields={'smooth','sm1','sm2','sm3','sm2j','raw','cbs','cbs_fixed','flag','hmm','affy_calls','ref','ads','qv','fxa','fxa_extra','pvs'};
if isfield(D,'matrix_fields')
  matrix_fields=[ D.matrix_fields];
end
for j=1:length(matrix_fields)
  if isfield(D,matrix_fields{j})
    x=getfield(D,matrix_fields{j});
    x=x(idx,:);
    D=setfield(D,matrix_fields{j},x);
  end
end
