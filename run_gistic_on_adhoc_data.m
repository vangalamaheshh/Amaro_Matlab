function run_gistic_on_adhoc_data (outdir,gistic_type,center_tag,seg_file, markers_file, array_list_file, ref_gene_file,cnv_file,t_amp,t_del, join_segment_size)
argnames = {'outdir','gistic_type','center_tag','seg_file','markers_file','array_list_file','ref_gene_file','cnv_file','t_amp','t_del','join_segment_size'};
version.last_changed_revision = '$LastChangedRevision$';
version.last_changed_date = '$LastChangedDate$';

if nargin<1
    %gistic_type = 'standard';
    gistic_type = 'focal';
    outdir = '/xchip/tcga_scratch2/gsaksena/gistic_henrik_combined_2009-05-27';
end


%Create output directory
timestamp = num2cell(clock);
timestamp_str = sprintf('%4d-%02d-%02d__%02d-%02d-%02.0f',timestamp{:});
outdir_run = [outdir '/' gistic_type '_run_' timestamp_str ];
mkdir (outdir_run);

%logging at start
logfilename = 'logfile.txt';
logfile = [outdir_run '/' logfilename];
[fid,msg] = fopen(logfile,'a'); 
if fid==-1
    throw(MException('run_gistic_on_adhoc_data:CantOpenFile',...
        ['Tried to open ' logfile ' but ' msg]));
end
fn = mfilename('fullpath');
d = dir([fn '.m']);

fprintf(fid,'===============================================\n');
fprintf(fid,'===============================================\n');
fprintf(fid,'output directory: %s\n',outdir_run);
fprintf(fid,'filename: %s\n',fn);
fprintf(fid,'last modified on disk: %s\n',d.date)
fprintf(fid,'last modified in repostory: %s\n',version.last_changed_date);
fprintf(fid,'last version in repository: %s\n',version.last_changed_revision);
fprintf(fid,'\n');
fprintf(fid,'arguments: \n');
for i=1:nargin
    fprintf(fid,'%s: %s\n',argnames{i},num2str(eval(argnames{i})));
end
fprintf(fid,'\n');
fprintf(fid,'stack trace\n');
st = dbstack('-completenames');
for i = 1:length(st)
    fprintf(fid,'file:%s,   name:%s,   line:%s\n',st(i).file,st(i).name,num2str(st(i).line));
end
fprintf(fid,'\n');
starttime = now;
fprintf(fid,'starting run at %s',datestr(starttime));
fprintf(fid,'\n');
fprintf(fid,'-------\n');
fprintf(fid,'\n');

fclose(fid);

%TODO add checks for correct input files 


[fid,msg] = fopen(logfile,'a'); 
if fid==-1
    throw(MException('run_gistic_on_adhoc_data:CantOpenFile',...
        ['Tried to open ' logfile ' but ' msg]));
end


try
%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the actual work...
    %run_gistic_on_submitted_data_work (center,center_tag,batches,outdir_run,gistic_type,cnv_file,fid)
    if nargin<3
        run_gistic_on_adhoc_data_work(fid,outdir_run,gistic_type);
    else 
        run_gistic_on_adhoc_data_work(fid,outdir_run,gistic_type,center_tag,seg_file, markers_file, array_list_file, ref_gene_file,cnv_file,t_amp,t_del, join_segment_size);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(fid);
    % log that work was completed
    [fid,msg] = fopen(logfile,'a');
    if fid==-1
        throw(MException('run_gistic_on_adhoc_data:CantOpenFile',...
            ['Tried to open ' logfile ' but ' msg]));
    end

    fprintf(fid,'cat %s; echo completed_log\n',logfile);
    fprintf(fid,'cd %s; echo completed_dir\n',outdir_run);
    
    fclose(fid);

    %update symbolic link to passing directory
    cd (outdir);
    cmd_str = sprintf('rm -f latest_%s_run',outdir_run,gistic_type);
    [status1,result1]=system (cmd_str);

    cmd_str = sprintf('ln -s %s latest_%s_run',outdir_run,gistic_type);
    [status2,result2]=system (cmd_str);
    if status1~=0 || status2~=0
        throw(MException('run_gistic_on_adhoc_data:LinkFailed',['link failed:  ' result1 result2]));
    end
catch ME
    fprintf(1,'%s\n',ME.getReport('extended','hyperlinks','on'));
    %dump exception to logfile.
    [fid,msg] = fopen(logfile,'a');
    if fid==-1
        throw(MException('run_gistic_on_adhoc_data:CantOpenFile',...
            ['Tried to open ' logfile ' but ' msg]));
    end
    fprintf(fid,'%s\n',ME.getReport('extended','hyperlinks','off'));
    fclose(fid);
end
    

%TODO add checks for correct output files

%logging at end
[fid,msg] = fopen(logfile,'a');
if fid==-1
    throw(MException('run_gistic_on_adhoc_data:CantOpenFile',...
        ['Tried to open ' logfile ' but ' msg]));
end

fprintf(fid,'\n');
fprintf(fid,'-------\n');
fprintf(fid,'\n');
endtime = now;
fprintf(fid,'ending run at %s\n',datestr(endtime));
totaltime = num2cell(datevec(endtime-starttime));
totaltime_str = sprintf('%4d-%02d-%02d__%02d-%02d-%02.2f',totaltime{:});
fprintf(fid,'total runtime: %s\n',totaltime_str);
fclose(fid)

end





%%

function run_gistic_on_adhoc_data_work(logfid,outdir_run,gistic_type,center_tag,seg_file, markers_file, array_list_file, ref_gene_file,cnv_file,t_amp,t_del, join_segment_size)
if nargin<4
    center_tag = 'henrik';
    seg_file = '/xchip/tcga_scratch2/ov_jamboree_2009/cn/MSCN,wCBS/segData,txt/TCGA,OV,mscn,wCBS,Chr1-23,1.0.0/All,MSCN,wCBS.seg.txt';
    markers_file = '/xchip/tcga_scratch2/ov_jamboree_2009/cn/MSCN,wCBS/annotationData2,txt/combined.markers2.txt';
    array_list_file = '/xchip/tcga_scratch2/ov_jamboree_2009/cn/MSCN,wCBS/segData,txt/TCGA,OV,mscn,wCBS,Chr1-23,1.0.0/samples_filtered2.txt';
    ref_gene_file = '/xchip/tcga/Annotation_Files/UCSC_hg18/hg18_with_miR_20080407.mat';
    cnv_file = '/xchip/tcga/gbm/analysis/mokelly/080429_convert_CNV_to_BED/CNV.verified_080606.combined.Ovarian.081202.txt';
    %gistic_type = 'standard';
    %gistic_type = 'focal';
    %outdir = '/xchip/tcga_scratch2/gsaksena/gistic_henrik_combined_2009-05-27';
    t_amp = .12;
    t_del = .12;
    join_segment_size = 3;
end

gistic_params_struct.t_amp = t_amp;
gistic_params_struct.t_del = t_del;
gistic_params_struct.join_segment_size = join_segment_size;
gistic_params_struct.qv_thresh=0.25;
gistic_params_struct.save_seg_data= 1;
gistic_params_struct.res=0.0001; % was 0.001 (may change in the future)
gistic_params_struct.remove_X=1;

fprintf(logfid,'seg_file: %s\n',seg_file);
fprintf(logfid,'markers_file: %s\n',markers_file);
fprintf(logfid,'array_list_file: %s\n',array_list_file);
fprintf(logfid,'ref_gene_file: %s\n',ref_gene_file);
fprintf(logfid,'cnv_file: %s\n',cnv_file);
fprintf(logfid,'\n');
 
 
fprintf(logfid,'t_amp: %s\n',num2str(gistic_params_struct.t_amp));
fprintf(logfid,'t_del: %s\n',num2str(gistic_params_struct.t_del));
fprintf(logfid,'join_segment_size: %s\n',num2str(gistic_params_struct.join_segment_size));
fprintf(logfid,'qv_thresh: %s\n',num2str(gistic_params_struct.qv_thresh));
fprintf(logfid,'save_seg_data: %s\n',num2str(gistic_params_struct.save_seg_data));
fprintf(logfid,'res: %s\n',num2str(gistic_params_struct.res));
fprintf(logfid,'remove_X: %s\n',num2str(gistic_params_struct.remove_X));
fprintf(logfid,'\n');


%%%%%%%%%%%%%%%%%%%%%%%%
    calling_gistic (seg_file,array_list_file,markers_file,ref_gene_file, cnv_file, gistic_params_struct,outdir_run,gistic_type,center_tag,logfid);
%%%%%%%%%%%%%%%%%%%%%%%%%%



end


function run_gistic_on_submitted_data_work (center,center_tag,batches,outdir,gistic_type,cnv_file,logfid)

sd = seg_data();
%fetch all batches to run
num_batches = length(batches);
for i = 1:num_batches
    data_perbatch{i} = sd.get_batch_data(center,batches(i));
end  

%concat all batches
data_all = vertcat(data_perbatch{:});


%scrub out uninteresting data
data_filt1 = seg_data.filter_keep_tumors_only_ds(data_all);
data_filt2 = ...
    seg_data.filter_remove_excessively_segmented_samples_ds(data_filt1);
scrubbed_data = data_filt2;

%fetch correct markers file
markers_file = seg_data.get_probe_file(center,batches(1));
fprintf(logfid,'markers_file: %s\n',markers_file);

%Select GISTIC parameters
platform = sd.get_platform(center,batches(1));
fprintf(logfid,'platform: %s\n',platform);

switch platform
    case {'HG-CGH-244A' 'CGH-1x1M_G4447A' 'Human1MDuo'}%harvard, mskcc1; mskcc2; stanford.
        gistic_params_struct.t_amp = 0.12;
        gistic_params_struct.t_del = 0.12;
        gistic_params_struct.join_segment_size = 3;
    case 'Genome_Wide_SNP_6' %broad
        gistic_params_struct.t_amp = 0.1699; %determined to be good for snp 6.0 platform by Gaddy.
        gistic_params_struct.t_del = 0.1926;
        gistic_params_struct.join_segment_size = 10;
    otherwise
        throw(MException('run_gistic_on_adhoc_data:unrecognizedPlatform','unsupported platform type'));
end

gistic_params_struct.qv_thresh=0.25;
gistic_params_struct.save_seg_data= 1;
gistic_params_struct.res=0.0001; % was 0.001 (may change in the future)
gistic_params_struct.remove_X=1;


fprintf(logfid,'t_amp: %s\n',num2str(gistic_params_struct.t_amp));
fprintf(logfid,'t_del: %s\n',num2str(gistic_params_struct.t_del));
fprintf(logfid,'join_segment_size: %s\n',num2str(gistic_params_struct.join_segment_size));
fprintf(logfid,'qv_thresh: %s\n',num2str(gistic_params_struct.qv_thresh));
fprintf(logfid,'save_seg_data: %s\n',num2str(gistic_params_struct.save_seg_data));
fprintf(logfid,'res: %s\n',num2str(gistic_params_struct.res));
fprintf(logfid,'remove_X: %s\n',num2str(gistic_params_struct.remove_X));

%reference input files
ref_gene_file = '/xchip/tcga/Annotation_Files/UCSC_hg18/hg18_with_miR_20080407.mat';
%cnv_file = '/xchip/tcga/gbm/analysis/mokelly/080429_convert_CNV_to_BED/CNV.verified_080606.combined.Ovarian.081202.txt';
%cnv_file = '/xchip/tcga/Annotation_Files/CNVs/CNV_coordinates.combined.Ovarian.090414.txt';

fprintf(logfid,'ref_gene_file: %s\n',ref_gene_file);
fprintf(logfid,'cnv_file: %s\n',cnv_file);

%%%%%%%%%%%%%%%%%%%%%%%%
     calling_gistic (seg_file,array_list_file,markers_file,ref_gene_file, cnv_file, gistic_params_struct,outdir_run,gistic_type,center_tag,logfid)
%%%%%%%%%%%%%%%%%%%%%%%%%%


end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calling_gistic (seg_file,array_list_file,markers_file,ref_gene_file, cnv_file, gistic_params_struct,outdir_run,gistic_type,center_tag,logfid)

dbstop if error
set_verbose_level(30)
%add paths
cd ~/CancerGenomeAnalysis/trunk/matlab
startup


cd(outdir_run);


t_amp = gistic_params_struct.t_amp;
t_del = gistic_params_struct.t_del;
join_segment_size = gistic_params_struct.join_segment_size;
qv_thresh = gistic_params_struct.qv_thresh;
save_seg_data = gistic_params_struct.save_seg_data;
res = gistic_params_struct.res;
remove_X = gistic_params_struct.remove_X;

switch gistic_type
    case 'standard'
        ext=['.' center_tag '.standard'];
        focal = 0;
        run_gistic_from_seg(outdir_run,seg_file,markers_file,ref_gene_file,array_list_file,cnv_file,t_amp,t_del, ...
            join_segment_size,ext,qv_thresh,save_seg_data,res,remove_X,focal,[-1.5 1.5]);
        
    case 'focal'
        ext=['.' center_tag '.focal'];
        focal = 1;
        run_gistic_from_seg(outdir_run,seg_file,markers_file,ref_gene_file,array_list_file,cnv_file,t_amp,t_del, ...
            join_segment_size,ext,qv_thresh,save_seg_data,res,remove_X,focal,[-1.5 1.5]);
        
    case 'peeloff.focal'
        ext=['.' center_tag '.peeloff.focal'];
        focal = 1;
        peeloff_focal = struct('k',1,'peel_off_method','regions');
        run_gistic_from_seg(outdir_run,seg_file,markers_file,ref_gene_file,array_list_file,cnv_file,t_amp,t_del, ...
            join_segment_size,ext,qv_thresh,save_seg_data,res,remove_X,focal,[-1.5 1.5],...
            peeloff_focal);
        
    otherwise
        throw(MException('run_gistic_on_adhoc_data:unknownGisticType',['unexpected value for Gistic_type: ' gistic_type]));
        
end
end
