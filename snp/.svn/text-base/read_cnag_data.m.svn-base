function C=read_cnag_data(d,C,save_name)

if ischar(d)
  fname=d;

  try
    fprintf(1,['Reading ' fname ' ...']);
%    [n,t]=xlsread(fname);
    [t,fid]=read_dlm_file(fname,char(9),2);
%    n=textscan(fid,'%s%s%s%s%s%s%s','bufSize',1000000);
    n=textscan(fid,'%f%f%f%f%f%f%f','bufSize',1000000);
    fclose(fid);
  catch
    disp('error....skipping file');
    C=[];
    return
    %    error('error reading excel file');
  end
  if (1)
    posn=find(n{2}>0);
    posn(posn>min(cellfun('length',n)))=[];
    C.pos=n{2}(posn);
    C.chrn=n{1}(posn);
    C.dat=n{4}(posn); % taking raw CNAG data
    %  C.flag=str2num(strvcat(n{3}(posn)));
    %  C.hmm=str2num(strvcat(n{5}(posn)));
    C.affy_call=n{7}(posn);
    C.sdesc=t{1}{1};    
  else
    posn=setdiff(find(str2num(strvcat(n{2}))>0),find(~cat(1,cellfun('isempty',regexp(n{4},'.*IN.*')))));
    C.pos=str2num(strvcat(n{2}(posn)));
    C.chrn=str2num(strvcat(n{1}(posn)));
    C.dat=str2num(strvcat(n{4}(posn))); % taking raw CNAG data
    %  C.flag=str2num(strvcat(n{3}(posn)));
    %  C.hmm=str2num(strvcat(n{5}(posn)));
    C.affy_call=str2num(strvcat(n{7}(posn)));
    C.sdesc=t{1}{1};
  end
  med=median(C.dat,1);
  disp(['subtracting ' num2str(med) ' and multiplying by 2']);
  C.dat=C.dat-med;
  C.dat=C.dat*2;
  fprintf(1,'%s\n',' Done.'); 
else
  if ~iscell(d)
    d={d.name};
  end
  for i=1:length(d)
    if ~exist('C','var') || isempty(C)
      C=read_cnag_data(d{i});
    elseif isfield(C,'sdesc') && ~isempty(findstrings(C.sdesc,d{i}))
      disp([ d{i} ' already exists']);
    else
      c=read_cnag_data(d{i});
      if ~isempty(c)
        C=unite_Ds({C,c},'cols');
      end
    end
    if mod(i,10)==0 && exist('save_name','var')
      save(save_name,'C');
    end
  end
  if exist('save_name','var')
    save(save_name,'C');
  end
end

