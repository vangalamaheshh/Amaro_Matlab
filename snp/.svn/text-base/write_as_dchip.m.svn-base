function write_as_dchip(fname,C,remove_top_line,ref_samples,gp_format,force_log,thresh_val)
% gp_format=-1 for snp format gp_format=-2 for snp format w/ chr and pos
%
%   History:
%          02 Oct 07:  Added error-catching for out-of-range affy_call
%          values.  Jen Dobson (jdobson@broad.mit.edu) 
%---
% $Id$ 
% $Date$ 
% $LastChangedBy$
% $Rev$

if ~exist('gp_format','var') 
  gp_format=0;
end

if ~exist('force_log','var') 
  force_log=0;
end

if ~exist('thresh_val','var') 
  thresh_val=0.1;
end


if gp_format>=0
  if floor(C.pos(min(find(~isnan(C.pos)))))==C.pos(min(find(~isnan(C.pos)))) % not in Mb
    C.pos=C.pos/1e6;
  end
end


if (abs(mean(C.dat(:)))>0.5 && (force_log==0)) || (force_log==1)
  % not log2
  if ~gp_format
    disp(['Taking log2(x)-1 of data after thresholding at ' num2str(thresh_val)]);
    C = preproc_log2trans(C,1,thresh_val);
  end
elseif (force_log==0) || (force_log==2)
  % log2
  if gp_format
    disp('Taking 2^(x+1) of data');
    C.dat=2.^(C.dat+1);
  end
end

if exist('ref_samples','var') && ~isempty(ref_samples)
  C.sdesc=[C.sdesc ref_samples];
  C.dat=[C.dat zeros(getsize(C,'dat',1),length(ref_samples))];
  C.dat=2.^(C.dat+1);
end

if ischar(C.sdesc)
  nms=cellstr(C.sdesc);
else
  nms=C.sdesc;
end


n=getsize(C,'dat',2);

f=fopen(fname,'w');

if ~exist('remove_top_line','var')
  remove_top_line=0
end

if remove_top_line==0
  fprintf(f,'Data type: Copy number\r\n');
end

nms2=cell(2*length(nms),1);
for i=1:length(nms)
  nms2{i*2-1}=nms{i};
  nms2{i*2}=[ nms{i} ' Call'];
end

if gp_format>0
  fmt=['Marker\tChromosome\tPosition (b)' repmat('\t%s',1,n*2) '\r\n'];
  fprintf(f,fmt,nms2{:});
  fmt=['%s\t%s\t%d' repmat('\t%f\t%s',1,n) ...
       '\r\n'];
  if ~isfield(C,'affy_calls')
    C.affy_calls=nan(getsize(C,'dat'));
  end
  C.affy_calls(isnan(C.affy_calls))=4;
  bad_calls = find(C.affy_calls<1 | C.affy_calls>4);
  if ~isempty(bad_calls)
      warning('One or more affy calls out-of-range.  Replacing with ''No Call''.')
      C.affy_calls(bad_calls) = 4
  end
  if gp_format==2 %affy calls (AA,BB,AB,NoCall)
    call_str={'AA','BB','AB','NoCall'};
  else
    call_str={'A','B','AB','NoCall'};
  end
elseif gp_format==0
  fmt=['Marker\tChromosome\tPosition (Mb)\tGenetic (cM)\tScore' repmat('\t%s',1,n) ...
       '\r\n'];
  fprintf(f,fmt,nms{:});
  fmt=['%s\t%s\t%f\t-100\t-100' repmat('\t%f',1,n) ...
       '\r\n'];
elseif gp_format==-1
  fmt=['probe set' repmat('\t%s',1,n*2) '\r\n'];
  fprintf(f,fmt,nms2{:});
  if ~isfield(C,'affy_calls')
    C.affy_calls=nan(getsize(C,'dat'));
  end  
  C.affy_calls(isnan(C.affy_calls))=4;
  bad_calls = find(C.affy_calls<1 | C.affy_calls>4);
  if ~isempty(bad_calls)
      warning('One or more affy calls out-of-range.  Replacing with ''No Call''.')
      C.affy_calls(bad_calls) = 4;
  end
  if gp_format==-1
    call_str={'A','B','AB','NoCall'};
  else
    call_str={'AA','BB','AB','NoCall'};    
  end
  fmt=['%s' repmat('\t%4.3f\t%s',1,n) '\r\n'];
elseif gp_format==-2
  if ~isfield(C,'affy_calls')
    C.affy_calls=nan(getsize(C,'dat'));
  end  
  C.affy_calls(isnan(C.affy_calls))=4;
  bad_calls = find(C.affy_calls<1 | C.affy_calls>4);
  if ~isempty(bad_calls)
      warning('One or more affy calls out-of-range.  Replacing with ''No Call''.')
      C.affy_calls(bad_calls) = 4;
  end
  fmt=['Marker\tChromosome\tPosition' repmat('\t%s',1,n) '\r\n'];
  fprintf(f,fmt,nms{:});
  fmt=['%s\t%s\t%d' repmat('\t%f',1,n) '\r\n'];  
end

if ischar(C.marker)
  disp('Converting .marker to cell array');
  C.marker=cellstr(C.marker);
end

if ~isfield(C,'chr')
  if isfield(C,'chrn')
    C.chr=num2chromosome(C.chrn);
  else
    if gp_format~=-1
      error('No chromosome field in data structure')
    end
  end
end
  
if isfield(C,'chr') && ischar(C.chr)
  disp('Converting .chr to cell array');
  C.chr=cellstr(C.chr);
end

if gp_format==-2


  for i=1:getsize(C,'dat',1)
    dat_str = sprintf('\t%f',C.dat(i,:));
    out_str = sprintf('%s\t%s\t%d%s\r\n',C.marker{i},C.chr{i},C.pos(i),dat_str);
    fprintf(f,'%s',out_str);
    if mod(i,100000)==0
      fprintf(1,'.%d',i);
    end
  end

else
  for i=1:getsize(C,'dat',1)
    if gp_format>0
      c=mat2cell(C.dat(i,:),1,ones(getsize(C,'dat',2),1)); 
      c(2,:)=call_str(C.affy_calls(i,:));
      fprintf(f,fmt,C.marker{i},C.chr{i},round(C.pos(i)*1e6),c{:});
    elseif gp_format==0
      fprintf(f,fmt,C.marker{i},C.chr{i},C.pos(i),C.dat(i,:));
    elseif gp_format<0
      c=mat2cell(C.dat(i,:),1,ones(getsize(C,'dat',2),1)); 
      c(2,:)=call_str(C.affy_calls(i,:));
      fprintf(f,fmt,C.marker{i},c{:});   
    end
    if mod(i,10000)==0
      fprintf(1,'.%d',i);
    end
  end
end
fprintf(1,'.%d',i);
fprintf(1,'\n');

fclose(f);

