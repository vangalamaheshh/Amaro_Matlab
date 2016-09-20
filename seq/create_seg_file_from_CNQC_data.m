function create_seg_file_from_CNQC_data(individual_name, tumor_rcl,normal_rcl,tumor_lanelist,normal_lanelist,...
        lane_blacklist,region_list,normals_db,outfile)

if nargin<9
  error(['input parameters required:\n'...
        '   individual_name, tumor_rcl, normal_rcl, tumor_lanelist, normal_lanelist,\n'...
        '   lane_blacklist, region_list, normals_db, outfile\n']);
end

% load and process data
fprintf('Loading data\n');
[Lt Ln R T N Z] = load_and_process_CNQC_data(tumor_rcl,normal_rcl,tumor_lanelist,normal_lanelist,...
  lane_blacklist,region_list,normals_db);

% write segfile
S = [];
S.sample = repmat({individual_name},slength(R),1);
S.chromosome = R.chr;
S.start = R.start;
S.end = R.end;
S.notused = zeros(slength(R),1);
S.log2copyratio = log2(median(T(:,Lt.enoughreads),2));

fprintf('Writing %s\n',outfile);
save_struct(S,outfile);
