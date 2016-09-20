function [lst,map_chrn,map_pos]=mark_cnv(mapping,xrt)

xrt_chrn=chromosome2num(xrt.dat(:,3));
%xrt_st=str2num(strvcat(cellfun_any('replace_empty(x,NaN)',xrt.dat(:,4))));
xrt_st = xrt.dat(:,4);
xrt_st(cellfun(@isempty,xrt_st)) = {'nan'};
xrt_st = cellfun(@str2num,xrt_st);
%xrt_en=str2num(strvcat(cellfun_any('replace_empty(x,NaN)',xrt.dat(:,5))));
xrt_en = xrt.dat(:,5);
xrt_en(cellfun(@isempty,xrt_en)) = {'nan'};
xrt_en = cellfun(@str2num,xrt_en);

tmp=grep('---',mapping.dat(:,4),1);
mapping.dat(tmp,4)=cellstr(repmat('NaN',length(tmp),1));
map_chrn=chromosome2num(mapping.dat(:,4));
tmp=find(cellfun('isempty',mapping.dat(:,5)));
if ~isempty(tmp)
  mapping.dat(tmp,5)=cellstr(repmat('NaN',length(tmp),1));
end
tmp=grep('---',mapping.dat(:,5),1);
mapping.dat(tmp,5)=cellstr(repmat('NaN',length(tmp),1));
map_pos=str2num(strvcat(mapping.dat(:,5)));

lst=[];
for i=1:length(xrt_st)
  if ~isnan(xrt_chrn(i)+xrt_st(i)+xrt_en(i))
    idx=find(map_chrn==xrt_chrn(i) & map_pos>=xrt_st(i) & map_pos<=xrt_en(i));
    if ~isempty(idx)
      lst=[lst; mapping.dat(idx,1) cellstr(repmat(xrt.dat(i,1),length(idx),1))];
    end
  end
  if mod(i,100)==0
    disp(i)
  end
end
