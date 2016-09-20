function write_soft_file(fname,D,sample_info,series_name,series_info,write_dat)
% WRITE_SOFT_FILE

if ~exist('write_dat','var')
  write_dat=1;
end

fid=fopen(fname,'w');
%platform
% MISSING

nl=sprintf(newline);
%samples
strs1=get_D_string(D,1:size(D.dat,2),'A ${CC} sample from the ${TT}.');
strs2=get_D_string(D,1:size(D.dat,2),'${SSC}');
strs3=get_D_string(D,1:size(D.dat,2),'${SP}');
for i=1:size(D.dat,2)
  fprintf(fid,'^SAMPLE=%s%s',deblank(D.sdesc(i,:)),nl);
  
  fprintf(fid,'!Sample_title = %s%s',deblank(D.sdesc(i,:)),nl);
  fprintf(fid,'!Sample_type = single channel%s',nl);
  fprintf(fid,'!Sample_organism = %s%s',strs3{i},nl);
  fprintf(fid,'!Sample_description = %s%s',strs1{i},nl);
  fprintf(fid,'!Sample_target_source = %s%s',strs2{i},nl);
  % must have also !Sample_platform_id
  fnms=fieldnames(sample_info);
  for j=1:length(fnms)
    x=getfield(sample_info,fnms{j});
    if iscell(x)
      for k=1:length(x)
        fprintf(fid,'!Sample_%s = %s%s',fnms{j},x{k},nl);
      end
    else
      fprintf(fid,'!Sample_%s = %s%s',fnms{j},x,nl);
    end
  end
  
  % DATA
  if write_dat
    fprintf(fid,'#ID_REF = %s',nl);
    if isfield(D,'mvdat')
      for vi=1:size(D.values,1)
        fprintf(fid,['#%s = %s%s'],D.values{vi,1},D.values{vi,2},nl);        
      end
      fprintf(fid,'ID_REF');
      for vi=1:size(D.values,1)
        fprintf(fid,'\t%s',D.values{vi,1});
      end      
      fprintf(fid,'%s',nl);
      for j=1:size(D.dat,1)
        fprintf(fid,'%d',j);
        for vi=1:size(D.values,1) 
          fprintf(fid,'\t%8.2f',D.mvdat(j,i,vi));
        end
        fprintf(fid,'%s',nl);
      end      
    else
      if isfield(D,'values')
        fprintf(fid,['#%s = %s%s'],D.values{1,1},D.values{1,2},nl);
        fprintf(fid,'ID_REF\t%s%s',D.values{1,1},nl);
      else
        fprintf(fid,'#VALUE = Median Fluorescent Intensity%s',nl);
        fprintf(fid,'ID_REF\tVALUE%s',nl);
      end        
      for j=1:size(D.dat,1)
        fprintf(fid,'%d\t%8.2f%s',j,D.dat(j,i),nl);
      end
    end
  end
end

%series
fprintf(fid,'^SERIES=%s%s',series_name,nl);
% must have title, type, description
fnms=fieldnames(series_info);
for j=1:length(fnms)
  x=getfield(series_info,fnms{j});
  if iscell(x)
    for k=1:length(x)
      fprintf(fid,'!Series_%s = %s%s',fnms{j},x{k},nl);
    end
  else
    fprintf(fid,'!Series_%s = %s%s',fnms{j},x,nl);
  end
end
for i=1:size(D.dat,2)
  fprintf(fid,'!Series_sample_id = %s%s',deblank(D.sdesc(i,:)),nl);
end

fclose(fid);
