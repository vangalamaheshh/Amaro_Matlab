function D=merge_xba_hind(xba_file,hind_file,out_file)

disp(['reading ' xba_file]);
X=read_modelled_data_file(xba_file);

disp(['reading ' hind_file]);
H=read_modelled_data_file(hind_file);

gi_file='/xchip/data/gadgetz/snp/Rameen/data/release_120k_genome_info_ming_sorted.txt';
disp(['reading genome info from ' gi_file]);
GI=read_genome_info_file(gi_file);
GI.dat=cat(1,GI.dat{:});

Xr=1:length(X.marker);
Hr=(length(Xr)+1):size(GI.dat,1);

% add chromosome string
X.chr=GI.dat(Xr,2);
pos_w_nan=GI.dat(Xr,3);
idx_w_nan=find(cellfun('isempty',pos_w_nan));
pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
X.pos=str2num(strvcat(pos_w_nan));
X=add_chrn(X);

H.chr=GI.dat(Hr,2);
pos_w_nan=GI.dat(Hr,3);
idx_w_nan=find(cellfun('isempty',pos_w_nan));
pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
H.pos=str2num(strvcat(pos_w_nan));
H=add_chrn(H);

[M,m1,m2]=match_string_sets(regexprep(X.sdesc,'100KX',''),regexprep(H.sdesc,'100KH',''));
disp([ 'matched ' num2str(length(m1)) ' samples']);
X1=reorder_D_cols(X,m1);
H1=reorder_D_cols(H,m2);

keyboard
% join and interleave X and H
D=unite_Ds({X1,H1},'rows');
D=order_by_pos(D);
if exist('out_file','var') && ~isempty(out_file)
    disp(['writing ' out_file]);
    write_as_dchip(out_file, D);
end
