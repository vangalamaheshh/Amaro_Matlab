function [q,st,s,Y2,cn_vs_snp_noise]=copy_signal_to_noise(C,params,cyto)

s=calc_copy_quality_score(C,'medianabs',0);

if ~exist('params','var') || isempty(params)
  params.n_thresh=500;
end

C=add_cyto(C,cyto);
C.chrarmn=C.armn+C.chrn*2;
collapse_type='nanmedian';

Y = copyD(C);

Y=collapse_D(Y,'chrarmn',collapse_type); % 'mode'

Y=reorder_D_rows(Y,find(Y.n_collapse>params.n_thresh));

Y2=copyD(Y);

% Y2.dat=2.^(Y.dat+1); % move to raw space

%st = nanstd(Y2.dat,0,1);
st=itrfcn1(Y2,'dat',1,@nanstd,0,1);

q=st./s;

% Calculate noise at the interface between CN and SNP probes (mokelly,
% 080513)
probe_names = strvcat(C.marker);
probe_types = probe_names(:,1:3); % First 3 letters distinguish CN from SNP
num_probes = length(probe_names);
% Find indexes of interfaces
cn_vs_snp_interface_idx = (probe_types(1:num_probes-1)~=probe_types(2:num_probes))';
% Find chromosome interfaces, remove from cn_vs_snp list
chr_interface_idx = C.chrn(1:num_probes-1)~=C.chrn(2:num_probes);
cn_vs_snp_interface_idx = cn_vs_snp_interface_idx & ~chr_interface_idx;
% Calculate diffs at the interfaces, and take the medianabs
interface_diffs = diff(C.dat(cn_vs_snp_interface_idx, :));
cn_vs_snp_noise = nanmedian(abs(interface_diffs),1)/(sqrt(2)*norminv(0.75));


deleteDfiles(Y);