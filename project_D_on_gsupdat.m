function P=project_D_on_gsupdat(D,proj_type)

P.sdesc=D.sdesc;

if isfield(D,'supacc')
  P.supacc=D.supacc;
end

if isfield(D,'supdesc')
  P.supdesc=D.supdesc;
end

if isfield(D,'supdat')
  P.supdat=D.supdat;
end

P.gacc=D.gsupacc;
P.gdesc=D.gsupdesc;

if ischar(proj_type)
  tmp=struct('method',proj_type);
  proj_type=tmp;
end

switch proj_type.method
 case 'sum'
  P.dat=(D.dat'*(D.gsupdat'>0))';
 case 'sum_of_normalized'
  D=preprocess_D(D,'row_center_and_normalize');
  P.dat=(D.dat'*(D.gsupdat'>0))'; % D.gsupdat*D.dat;
 case 'mean_of_col_normalized'
  D=preprocess_D(D,'row_center');
  D=preprocess_D(D,'col_center_and_normalize');
  P.dat=D.gsupdat*D.dat./repmat(sum(D.gsupdat,2),1,size(D.dat,2));
 otherwise
  error('no such method');
end

