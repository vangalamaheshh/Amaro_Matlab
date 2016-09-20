function chip=read_chip_file(fname,chip_name)

tab=read_table(fname,[],[char(9) ','],1,'bufsize',10000000);

c={'symb',{'Symbol'};'probeset',{'Probe'}; 'desc',{'Title','Description'}};

tabs=tab.headers{1};
tabsi=find(~cat(1,cellfun('isempty',tabs)));

chip=[];
for i=1:size(c,1)
  cols=[];
  for j=1:length(c{i,2})
    cols=[cols find(~cat(1,cellfun('isempty',regexpi(tabs,lower(c{i,2}{j})))))];
  end
  m{i}=cols;
  if length(m{i})>1
    disp('Matched more than one col');
  else
    chip=setfield(chip,c{i,1},tab.dat{m{i}});
  end
end

if exist('chip_name','var')
  chip.chip=cellstr(repmat(chip_name,length(chip.probeset),1));
end

notmatch=setdiff(tabsi,cat(1,m{:}));
if ~isempty(notmatch)
  disp(['Did not match cols: ' num2str(notmatch)]);
end

