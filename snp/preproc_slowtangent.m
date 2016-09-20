function [C,components] = preproc_slowtangent(C,atatime)
%PREPROC_SLOWTANGENT - slower tangent for normalizing normals
%
%   [D,COMPONENTS] = preproc_slowtangent(C,ATATIME)
%
%  D is a copy number data structure containing sample data to normalize
%  in D.dat. D.supdat must contain an 'N' field marking the normals. D is
%  returned with normalized tumor and normal data. The returned COMPONENTS 
%  is a matrix of projection weights.
%
%  ATATIME is a group size for the number of normal samples to normalize
%  against their complement. The larger the group size the faster it will
%  run at the expense of a lower quality normalization
%

%% normalize normals

normalidxs = find(C.supdat(strmatch('N',C.supacc),:));
N_normals = length(normalidxs);

% initialize array to hold normalized normal data
normdata = nan(size(C.dat,1),N_normals);

if ~exist('atatime','var')
    atatime = max(1,floor(N_normals/15));
end
% calculate how many passes will be required
nbunches = max(2,floor(N_normals/atatime)); %!

% Rameen's "ideal male" trick
ideal_male = 0.0 - (C.chrn==23); %!!!HUMGEN specific to human species

% array that keeps track of weights of each sample in the projection
% each column has the weights of the sample vector. The first row is 
% the weight of the ideal male vector.
components = nan(N_normals+1,getsize(C,'dat',2));

counter = 0;
tic;
% loop across bunches of normals to normalize against the rest
for i = 1:nbunches
    % attempt to spread the samples we're normalizing across plates
    leave_out = i:nbunches:N_normals;
    leave_in = setdiff(1:N_normals,leave_out);
    
    % extract the normal samples from the data
    N = double(C.dat(:,normalidxs));
    % pick the bunch of normals to normalize
    X = N(:,leave_out);
    % build a normalizing plane from the remaining normals
    N = N(:,leave_in);
    all_center = mean(N,2);
    N = [ ideal_male, N-repmat(all_center,1,size(N,2)) ];

    q = []; % uninitialized inverse plane
    for j = 1:length(leave_out)
        % normalize using tangent method
        [x,q,w] = tangent_normalization(X(:,j)-all_center,N,q);
        % save normalized vector
        normdata(:,leave_out(j)) = x;
        % save weights (ideal male is first)
        components([1,leave_in+1],normalidxs(leave_out(j))) = w;
     end
    counter = counter + length(leave_out);
    verbose('normalized %d out of %d normals (%0.2f sec)',20,counter,N_normals,toc);
end

%% normalize tumors

% index the tumors
tumoridxs = setdiff(1:size(C.dat,2),normalidxs);
N_tumors = length(tumoridxs);

% grab a copy of all unnormalized data from C.dat
N = double(C.dat);
% store normalized normals in C.dat
C.dat(:,normalidxs) = normdata;

% split unnormalized data into normals and tumors
X = N(:,tumoridxs);
N = N(:,normalidxs);
% prepare the plane
all_center = mean(N,2);
N = [ideal_male N];

% re-initialize array to hold normalized tumor data
normdata = nan(size(C.dat,1),N_tumors);
verbose('Normalizing Tumors (at %0.2f sec)',10,toc)
q = [];
for k = 1:N_tumors;
    [x,q,w] = tangent_normalization(X(:,k)-all_center,N,q);
    % save normalized vector
    normdata(:,k) = x;
    % save weights
    components(:,tumoridxs(k)) = w;
    if ~mod(k,atatime)
        verbose('normalized %d out of %d tumors (%0.2f sec)',20,k,N_tumors,toc);
    end
end
verbose('Finished slow tangent normalization in %0.1f seconds.',20,toc)
% store the normalized tumor data
C.dat(:,tumoridxs) = normdata;

%!!! temporary !!!
save('components_111111.mat','components');
