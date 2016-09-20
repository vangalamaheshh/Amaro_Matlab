function add_symbol_to_gene_set(fname_in,annot_in,fname_out)

if ~exist(fname_in,'file')
  disp([ fname_in ' is missing']);
  return;
end

if ~exist(annot_in,'file')
  disp([ annot_in ' is missing']);
  return;
end

if strcmp(fname_in((end-2):end),'txt')
  pw_xls=read_dlm_file(fname_in,char(9));
  pw_xls=cat(1,pw_xls{:});
else
  [n,pw_xls,raw]=xlsread(fname_in);
  pw_xls=make_text_xlsread(raw);
end

tt=pw_xls(2:end,:);
xx=cell(size(tt,1),1);
for i=1:size(tt,1)
  xx{i}=sprintf(repmat('%s ; ',1,size(tt,2)),tt{i,:});
end
yy=pw_xls(2:end,2);

annot=read_dlm_file(annot_in,char(9));
maxlen=max(cellfun('length',annot));
for i=1:length(annot)
  if length(annot{i})<maxlen
    annot{i}=[ annot{i} cell(1,maxlen-length(annot{i}))];
  end
end
annot=cat(1,annot{:});

% if there is a symbol column .... do nothing

if any(~cellfun('isempty',regexpi(pw_xls(1,:),'symbol')))
  f=fopen(fname_out,'w');
  for i=1:size(pw_xls,1)
    fprintf(f,['%s' repmat('\t%s',1,size(pw_xls,2)-1) '\n'],pw_xls{i,:});
  end    
  fclose(f);
else
  f=fopen(fname_out,'w');
  fprintf(f,['%s' repmat('\t%s',1,size(pw_xls,2)) '\n'],pw_xls{1,:},'Gene Symbol');
  for i=2:size(pw_xls,1)
    if ~isempty(yy{i-1})
      pos1=strmatch(yy{i-1},annot(:,1));
    else 
      pos1=[];
    end
    if ~isempty(xx{i-1})
      pos2=strmatch(xx{i-1},annot(:,1));
    else
      pos2=[];
    end
    pos=union(pos1,pos2);
    disp(pos');
    if ~isempty(pos)
      onepos=find(~cellfun('isempty',regexp(annot(pos,3),'1')));
      disp(pos(onepos)');
      if isempty(onepos)
        fprintf(f,['%s' repmat('\t%s',1,size(pw_xls,2)) '\n'],pw_xls{i,:},'---');
      else
        fprintf(f,['%s' repmat('\t%s',1,size(pw_xls,2)) '\n'],pw_xls{i,:},annot{pos(onepos(1)),4});
      end
    end
  end    
  fclose(f);
end




