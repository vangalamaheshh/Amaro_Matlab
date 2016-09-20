function l=lsf(pathname,prevrun,use_compiler)


[s,w]=unix('which bsub');
spos=findstr(w,'not found');
if ~isempty(spos)
  l.has_lsf=0;
else
  l.has_lsf=1;
  rand('state',sum(100*clock));
  if nargin>=2 && ~isempty(prevrun)
    l.lsf_uid=['LSF_' prevrun];
  else
    l.lsf_uid=['LSF_' rand_uid(8)];
  end
end

if l.has_lsf
  if nargin>=1
    l.lsf_path=pathname;
  else
    l.lsf_path='./';
  end 
else
  verbose('No bsub function');
end

l.jobs=[];
l.maxjobs=Inf;
l.res={};

if exist('use_compiler','var')
  l.use_compiler=use_compiler;
else
  l.use_compiler=0; % compile m-files before running
end

l=class(l,'lsf');

