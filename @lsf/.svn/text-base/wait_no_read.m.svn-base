function n_running=wait_no_read(l,hvec,return_at_num,maxsec)
if (nargin==1) | isempty(hvec)
  hvec=find(l.jobs>0);
  verbose(['waiting and not reading for ' num2str(length(hvec)) ' jobs']);
end
if nargin<4
  maxsec=Inf;
end

tic;
nextdisp=10;
while ~isempty(hvec) & (toc<maxsec)
  % Poll for results
  new_hvec=hvec;
  for i=1:length(hvec)
      basename=[ l.lsf_path get_basename(l,hvec(i))];
      if exist([basename '_finished.txt'],'file')
        new_hvec=setdiff(new_hvec,hvec(i));
      end
  end
  if length(new_hvec)<=return_at_num
    return
  end
  hvec=new_hvec;
%  disp([ num2str(new_hvec) ' are running' ]);
  if toc>nextdisp
    lsf_update(l);
    nextdisp=nextdisp+10;
  end
  pause(0.1);
end

