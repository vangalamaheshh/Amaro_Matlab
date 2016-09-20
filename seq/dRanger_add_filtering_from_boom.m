function dRanger_add_filtering_from_boom(infile,t_boomdir,n_boomdir,outfile,params)
% Mike Lawrence 2009-09

if ~exist('params','var'), params=[]; end

fprintf('Loading input file\n');
X = load_struct(infile);
X = make_numeric(X,{'chr1','min1','max1','chr2','min2','max2'});

X = dRanger_add_filtering_from_boom_part2(X,t_boomdir,n_boomdir,params);

fprintf('Saving output file\n');
save_struct(X,outfile);

fprintf('Done\n');
