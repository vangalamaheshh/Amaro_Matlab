function gp_gistic_from_seg(varargin)
%
% -b base_dir -seg segfile  -mk markers_file 
%  -refgene refgenefile [-alf array_list_file] [-cnv cnv_file] 
%  [-ta t_amp] [-td t_del] [-alldat has_allele_data] 
%  [-excols snpfile_has_extra_columns] [-snphead snp_header_rows] 
%  [-qvt qv_thresh] [-alf array_list_file] [-js join_segment_size(def=5)] 
%  [-ext extension] [-cap cap_vals(def=[-Inf Inf])]
%
%  COMPILE INSTRUCTIONS:
%       0) use matlab7.5.  There are library installation problems that made this
%       unsuccessful with matlab 7.6  To use matlab 7.5 type 'reuse
%       matlab75' and then 'matlab' at the unix prompt.
%
%       1) run 'mbuild -setup'.  Verify that the mbuildopts.sh that gets
%       installed in ~/.matlab/R2007b references the gcc in /util/gcc-4.1.1
%        (in instances where CC is set to 'gcc', make sure that the gcc is
%        '/util/gcc-4.1.1/bin/gcc'.)
%
%       2) go to the directory where you want to compile and type:
%           mcc -v -m -w enable gp_gistic_from_seg
%
%       3) After compilation, you are ready to execute.  Set the
%       LD_LIBRARY_PATH environment variable by typing:
%           'setenv LD_LIBRARY_PATH
%           /broad/tools/apps/matlab75/sys/os/glnxa64:/broad/tools/apps/ma
%           tlab75/bin/glnxa64:/util/gcc-4.1.1/lib:$LD_LIBRARY_PATH'
%
%       4) execute the compiled binary
%  
%---
% $Id$
% $Date: 2008-08-12 12:35:31 -0400 (Tue, 12 Aug 2008) $
% $LastChangedBy: jdobson $
% $Rev$


if isempty(varargin)
    disp('Usage: gp_gistic -b base_dir -seg segmentation_file ')
    disp('-mk markers_file, -refgene ref_gene_file, [-alf array_list_file(def:empty)]')
    disp('[-cnv cnv_file] [-ta amplifications_threshold(def=.1)] [-td deletions_threshold(def=.1)]')
    disp('[-js join_segment_size(def=4)] [-ext extension] [-qvt qv_thresh(def=0.25)]')
    disp('[-rx remove_x(def=1)] [-v verbosity_level(def=0)]') 
    return
end

method_st = 'GISTIC';

a = handle_args({'b','seg','mk','refgene','alf','cnv','ta','td','js','ext',...
    'qvt','rx','v'},varargin);

base_dir = a.b;
if isempty(base_dir)
    error('Must supply a base directory.')
end

if ~exist(base_dir,'dir')
    error('Base directory does not exist');
end


segfile = a.seg;
if isempty(segfile)
    error('Must supply a segmentation file.');
end

markersfile = a.mk;
if isempty(markersfile)
    error('Must supply a markers file.');
end

refgenefile = a.refgene;
if isempty(refgenefile)
    error('Must supply a refgene file.');
end

if ~isempty(a.alf)
    array_list_file = a.alf;
else
    array_list_file = [];
end

if ~isempty(a.cnv)
    cnv_file = a.cnv;
else
    cnv_file = [];
end


if isempty(a.ta)
    t_amp = .1;
else
    t_amp = str2double(a.ta);
end

if isempty(a.td)
    t_del = .1;
else
    t_del = str2double(a.td);
end



if isempty(a.js)
    join_segment_size = 4;
else
    join_segment_size = str2double(a.js);
end



if isempty(a.ext)
    ext=[];
else
    ext=a.ext;
end


if isempty(a.qvt)
    qv_thresh = 0.25;
else
    qv_thresh = str2double(a.qvt);
end


save_seg_data = 0;
res = [];

if isempty(a.rx)
    removex = 1;
else
    removex = str2double(a.rx);
end
% 
% if isempty(a.fa)
     focal_analysis = 0;
% else
%     focal_analysis = str2double(a.fa);
% end

caps = [Inf Inf];


if ~isempty(a.v)
    set_verbose_level(str2double(a.v));
end



% 
% if ~isempty(a.lcap)
%     caps(1) = str2double(a.lcap);
% end
% 
% if ~isempty(a.ucap)
%     caps(2) = str2double(a.ucap);
% end

    

%% RUN_GISTIC_FROM_SEG

verbose('Parameters:',10)
verbose([      'base_dir = ' base_dir],10);
verbose([      'segfile = ' segfile],10);
verbose([      'markersfile = ' markersfile],10);
verbose([      'refgenefile = ' refgenefile],10);
verbose([      'array_list_file = ' array_list_file],10);
verbose([      'cnv_file = ',cnv_file],10);
verbose([      't_amp = ',num2str(t_amp)],10);
verbose([      't_del = ',num2str(t_del)],10);
verbose([      'join_segment_size = ', num2str(join_segment_size)],10);
verbose([      'ext = ', ext],10);
verbose([      'qv_thresh = ',num2str(qv_thresh)],10);
verbose([      'remove x = ', num2str(removex)],10);

param_struct = run_gistic_from_seg(base_dir,segfile,markersfile,refgenefile,array_list_file,cnv_file,...
    t_amp,t_del,join_segment_size,ext,qv_thresh,save_seg_data,res,removex,focal_analysis,caps,[],1);

close all

