function GI=read_genome_info_file(fname,mode)
if ~exist('mode','var')
  d=read_dlm_file(fname);
  GI.dat=d(2:end);
  GI.header=d{1};
elseif strcmp(mode,'fast')
  f=fopen(fname);
  d=textscan(f,'%s%s%s%s%s%s','bufSize',10000000);
  fclose(f);
  d=cat(2,d{:});
  GI.header=d(1,:);
  GI.dat=d(2:end,:);
else
  error('no such mode');
end
