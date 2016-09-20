function fh_CopyNumberQCReport(individual_name,tumor_rcl,normal_rcl,tumor_lanelist,normal_lanelist,...
        lane_blacklist,region_list,normals_db,tumor_seg,normal_seg,params)
% fh_CopyNumberQCReport(individual_name,tumor_rcl,normal_rcl,tumor_lanelist,normal_lanelist,
%       lane_blacklist,region_list,normals_db,tumor_seg,normal_seg,params)
%
%  Inputs:
%
%  individual_name = e.g. "GBM-0188"
%
%  *_rcl = location of RCL file (= concatenated output of all RegionCovPerLane.java jobs)
%            (no header line)
%            <targname/genename> <chr> <start> <end> <lane0_count> <lane1_count> ... <laneN_count>
%
%  *.lanelist = location of output from MakeLanelist.java (Firehose version; i.e. Lanetable)
%
%  lane_blacklist = location of lane blacklist, containing list of blacklisted lanes
%
%  region_list = location of textfile containing the region definitions
%            (no header line)
%            <targname/genename> <chr> <start> <end> (<membership>)
%            (optional column <membership> is used for distinguishing C2K, C6K, WE datasets)
%
%  normals_db = location of file containing list of normal RCL files for doing median-of-normals normalization
%          
%  tumor_seg = location of tumor seg file (ground truth)
%           (THIS PARAMETER IS OPTIONAL: If omitted, will only be able to diagnose certain cases)
%               (a header line may be present or not)
%               <ID> <chr> <start> <end> <nprobes> <segmean>
%               (<ID> is ignored)
%  Outputs:
%
%  writes the following files to the cwd:
%    *_CopyNumberQC.png
%    *_CopyNumberQC_report.txt     * = <individual_name>
%    *_CopyNumberQC.seg.txt          = "segfile" for displaying on CIRCOS plot
%    report.html
%    num_mixups.txt    = single line telling count of possible mixups
%
% Mike Lawrence 2009-10-16

if nargin<8
  fprintf(['input parameters required:\n'...
        '   individual_name, tumor_rcl, normal_rcl, tumor_lanelist, normal_lanelist,\n'...
        '   lane_blacklist, region_list, normals_db, (tumor_seg), (normal_seg)\n']);
  error('Please supply all required parameters');
end
if ~exist('tumor_seg','var'), tumor_seg = []; end
if ~exist('normal_seg','var'), normal_seg = []; end
if ~exist('params','var'), params = []; end

params = impose_default_value(params,'ignore_XY',true);

fprintf('CopyNumberQC: %s\n',individual_name);

% load and process data
[Lt Ln R T N Z] = load_and_process_CNQC_data(tumor_rcl,normal_rcl,tumor_lanelist,normal_lanelist,...
  lane_blacklist,region_list,normals_db,tumor_seg,normal_seg,params);

% visualize
draw_CNQC_plot(Lt,Ln,R,T,N,Z,individual_name)
print('-dpng','-r150',[individual_name '_CopyNumberQC.png']);

% write "segfile"
S = [];
S.sample = repmat({individual_name},slength(R),1);
S.chromosome = R.chr;
S.start = R.start;
S.end = R.end;
S.notused = zeros(slength(R),1);
S.log2copyratio = log2(median(T(:,Lt.use),2));
outfile = [individual_name '_CopyNumberQC.seg.txt'];
fprintf('Writing %s\n',outfile);
save_struct(S,outfile);

% QC report: textfile
L = generate_CNQC_report(Lt,Ln,individual_name);
save_struct(L,[individual_name '_CopyNumberQC_report.txt']);

% QC report: HTML
H = convert_CNQC_report_to_html(L,Z,individual_name);
save_textfile(H,'report.html');

% n_mixups.txt
n_mixups = sum(L.is_mixup & ~L.is_blacklisted);
save_textfile(sprintf('%d\n',n_mixups),'num.mixups.txt');

% done
fprintf('Report complete.\n');
close all;



