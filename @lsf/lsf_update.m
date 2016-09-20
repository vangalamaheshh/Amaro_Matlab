function lsf_update(l)

[s,r]=unix(['date >> ' l.lsf_path l.lsf_uid '.log']);
%[s,r]=unix(['bjobs >> ' l.lsf_path l.lsf_uid '.log']);
st=[];
for i=1:length(l.jobs)
  if l.jobs(i)>0
    st=[st ' ' num2str(l.jobs(i))];
  end
end

[s,r]=unix(['bjobs ' st ' >> ' l.lsf_path l.lsf_uid '.log']);


if (0)
  figure(1); clf;
  text(0,1,strvcat(datestr(now),r),'VerticalAlignment','top');
  axis off;
  drawnow;
end
