function [M,P,hlarge,snps]=correct_batch_effect(M,batch_effect_params)
%CORRECT_BATCH_EFFECT corrects the batch effect for data in structure M.
%
%       [M,P,HLARGE,SNPS] = CORRECT_BATCH_EFFECT(M,BATCH_EFFECT_PARAMS)
%       inputs data structure, M, and BATCH_EFFECT_PARAMS with optional
%       fields:
%                   'min_sz' -- gives the minimum number of samples a batch
%                   needs in order to be corrected  (default = 5)
%
%                   'bonf_pv_threshold' -- gives the p-value threshold for
%                   bonferoni corrected batch data (default = 0.05)
%                   ***Bonferroni is never done with current code structure***
%
%                   'absolute_pv' -- gives the p-value threshold when
%                   bonferoni correction is not used  (default = 0.001)
%
%                   'conserve_memory' -- for large data sets, sets the
%                   conserve memory field  (not operational at  the moment)
%
%                   'cbediary_filename' -- location to write stdout
%                   from batch effect correction report.  If field exists,
%                   data will be written to filename specified.  If field
%                   does not exist, data will not be written.
%                    
%
%       Revisions:
%               9 Nov 07: Added support for inc_batch supdat.
%               (Documentation added, too).  (jendobson@broad.mit.edu)
%                   
%               19 Nov 07: Added option to write line snp correct data to
%               file.
%
%               3 Dec 07:  In last part of function (when changing M.dat)
%               load all of dat into variable dat so that whole array is
%               written to HDF5 rather than elements of array.  (this is
%               faster -- need to figure out how to fix this)           
%---
% $Id$ 
% $Date$ 
% $LastChangedBy$
% $Rev$

%% Initialize parameters


if ~exist('maxsnps','var')
    maxsnps = [];
end


%set minimum batch size
if isfield(batch_effect_params,'min_sz')
  min_sz = batch_effect_params.min_sz;
else
  min_sz = 5;
end

%set p-value threshold for bonferoni corrected data
if isfield(batch_effect_params,'bonf_pv_threshold')
  bonf_pv_thresh = batch_effect_params.bonf_pv_threshold;
else
  bonf_pv_thresh = 0.05;
end

%set p-value threshold non-corrected data
if isfield(batch_effect_params,'absolute_pv')
  absolute_pv = batch_effect_params.absolute_pv;
else
  absolute_pv = .001;
end

%set conserve memory flag
if isfield(batch_effect_params,'conserve_memory') && batch_effect_params.conserve_memory
  conserve_memory=1;
end

% if iscell(M) && ~conserve_memory
%   clear M;
%   global M;
% end  
%   
%% Find batches with size >= min_sz; if using includes, get include samples

batch_supid=strmatch('BATCH',M.supacc,'exact');
%control_supid=strmatch('CTRL',M.supacc,'exact');  %no longer used

nbatch=max(M.supdat(batch_supid,:));
h=histc(M.supdat(batch_supid,:),1:nbatch);
hlarge=find(h>=min_sz);

M = check_array_includes(M,'batch');
include_supid = strmatch('inc_batch',M.supacc,'exact');

%% Loop through batches returning p-values using differential_analysis.m

P=zeros(size(M.dat,1),length(hlarge));
%DC=zeros(size(M.dat,1),length(hlarge));  %DC is no longer used!
b=zeros(length(hlarge),size(M.dat,2));
for i=1:length(hlarge)
    
  %%% FIXME: add inner loop on tissue type (include normals in each tissue type)
  
  
  
  if ~isempty(include_supid)  %if array list file specified include samples, only use those in differential analysis
      thisbatch = intersect(find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(include_supid,:)));
      otherbatches = intersect(find(M.supdat(batch_supid,:)~=hlarge(i)),find(M.supdat(include_supid,:)));
  else
      thisbatch = find(M.supdat(batch_supid,:)==hlarge(i));
      otherbatches = find(M.supdat(batch_supid,:)~=hlarge(i));
  end
  
  verbose(['(large) BATCH: ' num2str(i)])
  [p,s]=differential_analysis(M,thisbatch,otherbatches,struct('method','ttest_minvar','minvar',(0.4)^2),0);
                     
% if the number of control samples in a batch is greater than the min batch
% size, get the differential mean between controls from the other batches
%  if nnz(M.supdat(control_supid,M.supdat(batch_supid,:)==hlarge(i)))>0
%    DC(:,i)=mean(M.dat(:,find(M.supdat(control_supid,M.supdat(batch_supid,:)==hlarge(i)))),2)-...
%            mean(M.dat(:,find(M.supdat(control_supid,M.supdat(batch_supid,:)~=hlarge(i)))),2);
%  else
%    DC(:,i)=NaN;
%  end
%             

%   [p,s]=differential_analysis(M,find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(batch_supid,:)~=hlarge(i)),...
%                               struct('method','ttest_minvar', ...
%                                      'minvar',(0.4)^2),0);
    
  P(:,i)=p;  %P gives p-values down the snps by batch (P is #snps X #batches)
  
  if ~isempty(include_supid)
      b(i,intersect(find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(include_supid,:))))=1;  
      bghosts(i,find(M.supdat(batch_supid,:)==hlarge(i)))=1;  %These are all the samples in the batch (even ones who will be corrected but won't be included in the analysis)
  else
      b(i,find(M.supdat(batch_supid,:)==hlarge(i)))=1;  % b is #batches X #samples (1 if sample in batch, 0 else)
             
      bghosts = b;% 
   
                                                        %                                                     
  end                                               % 
end


%% Figure out which snps to correct

% if diary file is specified, write correction data to file
if isfield(batch_effect_params,'CBEdiary_filename')
    diary batcheffect_params.CBEdiary_filename
end

if exist('absolute_pv','var') && absolute_pv~=0
  verbose(['Finding SNPs with P-value < ' num2str(absolute_pv) ' for at least one batch'],10);
  [minv,mini]=min(P,[],2);
  snps=find(minv<absolute_pv); %these are the snps that will need to be corrected
%  keyboard
else
  verbose(['Finding SNPs with Bonferroni corrected P-value < ' num2str(bonf_pv_thresh) ' for at least one batch'],10);
  P2=min(P*size(P,1),1);
  [minv,mini]=min(P2,[],2);
  snps=find(minv<bonf_pv_thresh);  %these are the snps that will need to be corrected
end

verbose(['correcting overall ' num2str(length(snps)) ' SNPs'],10);
mini=mini(snps);  %for each snp, the number of the batch that needs correction
verbose('number of SNPs corrected in each (large) batch :',10);
verbose('\t \t %7.0f \t   %7.0f \n',10, [hlarge ;histc(mini,1:length(hlarge))']);

%% Correct the snps

%mini is an NX1 vector (N is number of snps).  Each element of mini gives
%the batch number to be corrected for that snp

%correct data in snp chunks
memchunks = getmemchunkdims(M,'dat',2);
chunksize = memchunks(1);

k = 0;

%while k < length(snps)
while k < getsize(M,'dat',1)
    %maxidx = min(length(snps),k+chunksize);
    maxidx = min(getsize(M,'dat',1),k+chunksize);
    dat = M.dat((k+1):maxidx,:);

    thisloopsnps = (intersect((k+1):maxidx,snps(find(snps<=maxidx)))-k)';
    thismini = mini(1:length(thisloopsnps));  %each index of mini gives the batch to use for that row (snp)
    mini = mini(setdiff(1:length(mini),1:length(thismini)));
    
    B=sparse(b(thismini,:));  % dim1:number of corrected snps; dim2:total samples;  element = 1 if snp for sample needs correcting
    nB=sum(B,2);  %number in Batch
    
    x = dat(thisloopsnps,:);  %subset of data needing correction
    m1=sum(double(x).*B,2)./nB;    %column vector: each row is mean of data in the bad batch for that snp
    m2=sum(double(x).*(1-B),2)./(size(B,2)-nB);  % mean of data not in that batch for that snp
    delta=(m1-m2)./m2;  %
    verbose('Chunk %.0f: \t %f \t %f',10,ceil((k+1)/chunksize),mean(abs(delta)),std(abs(delta)));

    Bghosts = sparse(bghosts(thismini,:));
    %hist(100*delta,50);

    x=single(double(x)-Bghosts.*repmat((m1-m2),1,size(B,2)));   %adjust mean s.t. bad snps get mean of not in batch
    dat(thisloopsnps,:) = x;
    M.dat((k+1):maxidx,:) = dat;
    k = k+chunksize;
   % keyboard
end


    % Turn off diary write 
if isfield(batch_effect_params,'CBEdiary_filename')
    diary off
end


            
