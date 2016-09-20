function [Lt Ln R T N Z] = load_cnqc_job(jobdir,P)
% [Lt Ln R T N Z] = load_cnqc_job(jobdir)
%
% given jobdir for a CopyNumberQCReport job, e.g.
%     /xchip/cga/gdac-prod/genepattern/jobResults/511633
%
% loads the copy-number QC data for examination

if ~exist('P','var'), P=[]; end

x = load_lines([jobdir '/gp_execution_log.txt']);

tumor_lanelist = decell(gpj(regexprep(grep('tumor.lanelist',x),'.*\s(\S*)$','$1')));
normal_lanelist = decell(gpj(regexprep(grep('normal.lanelist',x),'.*\s(\S*)$','$1')));
tumor_rcl = decell(gpj(regexprep(grep('tumor.rcl',x),'.*\s(\S*)$','$1')));
normal_rcl = decell(gpj(regexprep(grep('normal.rcl',x),'.*\s(\S*)$','$1')));
lane_blacklist = decell(gpj(regexprep(grep('lane.blacklist',x),'.*\s(\S*)$','$1')));
region_list = decell(gpj(regexprep(grep('region.list',x),'.*\s(\S*)$','$1')));
normals_db = decell(gpj(regexprep(grep('normals.db',x),'.*\s(\S*)$','$1')));
tumor_seg = decell(gpj(regexprep(grep('tumor.seg',x),'.*\s(\S*)$','$1')));
normal_seg = decell(gpj(regexprep(grep('normal.seg',x),'.*\s(\S*)$','$1')));

[Lt Ln R T N Z]  = load_and_process_CNQC_data(tumor_rcl,normal_rcl,...
  tumor_lanelist,normal_lanelist,lane_blacklist,region_list,normals_db,tumor_seg,normal_seg,P);

if nargout<6
  tmp = Lt;
  Lt = [];
  Lt.Lt = tmp; clear tmp;
  Lt.Ln = Ln; clear Ln;
  Lt.R = R; clear R;
  Lt.T = T; clear T;
  Lt.N = N; clear N;
  Lt.Z = Z; clear Z;
end
