function run_gistic_on_all(srcpath,outpath,rg_file,params,ncpu,min_samples)
% TODOC

if ~exist('ncpu','var') || isempty(ncpu)
    ncpu = 6;
end
if ~exist('min_samples','var') || isempty(min_samples)
    min_samples = 40;
end

% get all files in the source path that match the canonical D.*.mat pattern
files = dir(srcpath);
[~,~,~,dfiles] = regexp({files.name},'D\.(.*)\.mat');
dfiles = [dfiles{~cellfun(@isempty,dfiles)}];
% extract TCGA cancer type from file name
dnames = regexprep(dfiles,'^D\.','');
dnames = regexprep(dnames,'\.mat','');

ng = length(dfiles);
dstructs = cell(1,ng);
dsizes = zeros(1,ng);
for i=1:ng
    dstructs{i} = load_D([srcpath,dfiles{i}]);
    dsizes(i) = size(dstructs{i}.dat,2);
end

% eliminates types with fewer than 40 samples
keepers = dsizes >= min_samples;
dstructs = dstructs(keepers);
dnames = dnames(keepers);
dsizes = dsizes(keepers);
ng = sum(keepers);

%% Create batches of gistic runs 
% The algorithm used to get roughly even numbers of cumulative samples is,
% starting from the largest sample and working down to the smallest, put
% each sample in the batch having the least number of samples.

batch = cell(1,ncpu);         % D structs for each {CPU}{batch member}
paramses = cell(1,ncpu);      % params for each {CPU}(batch member)
outdirs = cell(1,ncpu);       % output directory for each {CPU}(batch member)
est_times = zeros(1,ncpu);    % estimated time of batch (arbitrary units)

% order gistic runs by descending number of samples
[~,order] = sort(dsizes,'descend');
for i=1:ng
    % put each (D struct, run parameters) set in the batch
    % with the least cumulative samples
    [~,next] = min(est_times);
    batch{next} = [batch{next},dstructs(order(i))];
    tcgatype = dnames{order(i)};
    params.ext = ['.',tcgatype];
    paramses{next} = [paramses{next}, params];
    outdirs{next} = [outdirs{next},{[outpath,tcgatype,'/']}];
    % adjust batch sizes
    nsamp = dsizes(order(i));
    est_times(next) = est_times(next) + nsamp*nsamp; %! sad but closer to true
%!  est_times(next) = est_times(next) + nsamp*log(nsamp);
end

%% start the batches of gistic runs
matlabpool('open',num2str(ncpu));
parfor i = 1:ncpu
    for j = 1:length(batch{i})
        D = batch{i}{j};
        params = paramses{i}(j);
        new_dir = outdirs{i}{j};
        run_gistic20(new_dir,D,rg_file,params);
    end
end
matlabpool('close');



