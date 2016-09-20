function color_names=read_color_names(fname)
if ~exist('fname','var')
  fname='~gadgetz/matlab/color_names.txt';
end

d=read_dlm_file(fname);
d=cat(1,d{:});

cols=[str2num(strvcat(d(:,2))) str2num(strvcat(d(:,3))) str2num(strvcat(d(:,4))) ]/255;
color_names=[ lower(d(:,1)) mat2cell(cols,ones(size(d,1),1))];
