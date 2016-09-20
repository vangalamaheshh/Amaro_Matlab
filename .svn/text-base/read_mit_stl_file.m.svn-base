function [stl,field]=read_mit_stl_file(fname)
%
% miRNA Sample Translation File Format
% suffix: .stl
%
% Line1: Head line containing annotation fields
% Line2 and on: Sample information
%
% Column1: PCR_Name
% Column2: Labeling_Name
% Column3: Sample_Name
% Column4: Use
% Column5 to an empty column: 0/1 based sample annotation
% 
% After empty column: Other annotations

f=read_dlm_file(fname);

field=[];
if length(f)>=3
  for j=1:length(f{1})
    field(j).abr=deblank(f{1}{j});
    if ~isempty(f{2}{j})
      field(j).full=deblank(f{2}{j});
    end
    if ~isempty(f{3}{j})
      field(j).isnum=strcmp(deblank(f{3}{j}),'N');
    end    
  end
end


for i=4:length(f)
  for j=1:length(f{i})
    if field(j).isnum
      if isempty(deblank(f{i}{j}))
        stl{i-3,j}=NaN;
      else
        stl{i-3,j}=str2num(f{i}{j});
      end
    else
      if isempty(deblank(f{i}{j}))
        stl{i-3,j}=[];
      else
        stl{i-3,j}=deblank(f{i}{j});
      end
    end
  end
  for j=(length(f{i})+1):length(field)
    if field(j).isnum
      stl{i-3,j}=NaN;
    else
      stl{i-3,j}=[];
    end
  end
end
