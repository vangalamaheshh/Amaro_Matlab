function V = load_sequenom_results(fname)

demand_file(fname);

V = load_struct(fname);
V = parse_in(V,'sequenom_assay_name','^(.*)\|(.*)',{'assay_set','assay_name'});
V = parse_in(V,'assay_name','^(chr.+)_(\d+)_(gD|gI|)_?(\d+)$',...
  {'assay_chr','assay_pos','assay_type','assay_num'},[2 4]);



