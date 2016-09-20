function [C,bad_normals_cell]=preproc_histqctumors(Cin,base_dir,output_dir,show_hist,hist_delta,hist_normal_sig)
%PREPROC_HISTQCNORMALS Histogram quality control on normals.
%
%   [COUT,BAD_NORMALS_CELL] = PREPROC_HISTQCTUMORS
%    (CIN, BASE_DIR, OUTPUT_DIR, SHOW_HIST, HIST_DELTA, HIST_NORMAL_SIG)
%   
%   Required inputs: CIN, BASE_DIR, OUTPUT_DIR
%   Defaults: SHOW_HIST = 0;  HIST_DELTA = 0.01; HIST_NORMAL_SIG = 0.05;
%
%  About the hist_qc in supdat:  -1 means no qc done on normals, 0 means ok
%  normal, n>0 means normal set as bad normal in qc step n;
%
%
        
%       Revisions
%           5 Dec 07:  Added work-around for datastruct can't produce list
%           from C.sis.ploidy
%---
% $Id$
% $Date: 2008-06-04 10:20:01 -0400 (Wed, 04 Jun 2008) $
% $LastChangedBy: jdobson $
% $Rev$

verbose('Performing histogram_qc on tumors.',10)

varlist1 = {Cin,base_dir,output_dir,show_hist,hist_delta,hist_normal_sig};
required= [1,1,1,0,0,0];
defaults = {'err','err','err','0','0.01','0.05'};

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx}])
    end
end




if ~iscell(Cin)
    Cout{1} = Cin;
    wascell = 0;
else 
    Cout = Cin
    wascell = 1;
end

clear Cin




for k = 1:length(Cin)
    
    C = Cin{k};

    bad_tumors_cell = {};

    one_peak=histogram_qc(C,[base_dir 'hist_qc_tumor.out.ps'],1,show_hist,1,0,hist_delta,hist_normal_sig);

    if isa(C,'datastruct')  %Work around since can't produce LIST
        ploidy = C.sis.ploidy;
    else
        ploidy = {C.sis.ploidy};
    end
    
    normals = strmatch({'2'},ploidy);

    tumors = setdiff(1:size(C.dat,2),normals);

    bad_tumors = intersect(one_peak,tumors);

    v = C.supdat(find_supid(C,'hist_qc'),:,:); %extra dimension

    v(bad_tumors,:,:) = 3;

    C.supdat(find_supid(C,'hist_qc'),:) = v;

    thow_out_n_names = C.sdesc(bad_tumors);
% 
%     save([output_dir 'flag_bad_tumors.mat'],'bad_tumors','throw_out_n_names');


    C = add_history(C,mfilename,bad_normals,throw_out_n_names);

    Cout{k} = C;

    bad_normals_cell{k} = bad_normals;

end

if ~wascell
    tmp = Cout{1};
    Cout = tmp;
end


