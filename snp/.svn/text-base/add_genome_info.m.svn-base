function D=add_genome_info(D,GI,isXba)

if isXba
  r=1:length(D.marker);
else
  r=(size(GI.dat,1)-length(D.marker)+1):(size(GI.dat,1));
end

disp(range(strvcat(D.marker)-strvcat(GI.dat(r,1))));

D.chr=GI.dat(r,2);
pos_w_nan=GI.dat(r,3);
idx_w_nan=find(cellfun('isempty',pos_w_nan));
pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
D.pos=str2num(strvcat(pos_w_nan));
