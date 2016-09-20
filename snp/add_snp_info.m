function M=add_snp_info(M,GI,range,match_method)
%ADD_SNP_INFO add snp information to data structure.
%
%   M = ADD_SNP_INFO(M,GI,RANGE,MATCH_METHOD)
%           MATCH_METHOD  < {'same_order','by_pos','by_name'}
%           RANGE: index vector of cells to use in GI file, where GI.dat{RANGE}{2} matches the
%           M.marker

if isempty(range)
  range=1:size(GI.dat,1);
end

if ~exist('match_method','var')
  match_method='same_order';
end

switch match_method
 case 'same_order'
  M.chr=GI.dat(range,2);
  pos_w_nan=GI.dat(range,3);
  idx_w_nan=find(cellfun('isempty',pos_w_nan));
  pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
  M.pos=str2num(strvcat(pos_w_nan));
  M=add_chrn(M);
 case 'by_pos'
  if ~isfield(M,'chrn')
    M=add_chrn(M);
  end
  pos1=M.chrn*1e11+M.pos;
  
  chr=GI.dat(range,2);
  pos_w_nan=GI.dat(range,3);
  idx_w_nan=find(cellfun('isempty',pos_w_nan));
  pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
  pos=str2num(strvcat(pos_w_nan));
  pos2=chromosome2num(chr)*1e11+pos;
  
  [spos1,si1]=sort(pos1);
  [spos2,si2]=sort(pos2);
  n1=length(pos1);
  n2=length(pos2);
  m1=1;
  m2=1;
  match=nan(n1,1);
  while (m1<=n1) && (m2<=n2)
    if spos1(m1)==spos2(m2)
      match(m1)=m2;
      m1=m1+1;
      m2=m2+1;
    else
      m2=m2+1;
    end
  end
  
  M.marker(si1)=GI.dat(si2(match),1);
  M.marker=as_column(M.marker);
  if ~isfield(M,'chr')
    M.chr=num2chromosome(M.chrn);
  end
 case 'by_name'
  [Mt,m1,m2]=match_string_sets_hash(M.marker,GI.dat(:,1));
  M=reorder_D_rows(M,m1);
  M.chr=GI.dat(m2,2);
  pos_w_nan=GI.dat(m2,3);
  idx_w_nan=find(cellfun('isempty',pos_w_nan));
  pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
  M.pos=str2num(strvcat(pos_w_nan));
  M=add_chrn(M);
end
