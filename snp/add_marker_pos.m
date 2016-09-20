function C=add_marker_pos(C,GI)

if length(C.marker)==59015  % Xba
  gi_idx=1:59015;
elseif length(C.marker)==57299 % Hind
  gi_idx=(59015+1):size(GI.dat,1);
else
  error('Not Xba, nor Hind');
end
% range(strvcat(C.marker)-strvcat(GI.dat(gi_idx,1)))

% [dum,m1,m2]=match_string_sets(C.marker,GI.dat(:,1));

C.chr=GI.dat(gi_idx,2);
pos_w_nan=GI.dat(gi_idx,3);
idx_w_nan=find(cellfun('isempty',pos_w_nan));
pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
C.pos=str2num(strvcat(pos_w_nan));
C=reorder_D_rows(C,find(~isnan(C.pos)));
C=add_chrn(C);
