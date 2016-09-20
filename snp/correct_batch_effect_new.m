function [M,P,hlarge]=correct_batch_effect_new(M,batch_effect_params)
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
% $Date: 2008-02-28 17:15:19 -0500 (Thu, 28 Feb 2008) $ 
% $LastChangedBy: jdobson $
% $Rev$

%% Initialize parameters


%set minimum batch size
if isfield(batch_effect_params,'min_sz')
  min_sz = batch_effect_params.min_sz;
else
  min_sz = 5;
end


%set p-value threshold non-corrected data
if isfield(batch_effect_params,'absolute_pv')
  absolute_pv = batch_effect_params.absolute_pv;
else
  absolute_pv = .001;
end



%% Find batches with size >= min_sz; if using includes, get include samples

batch_supid=strmatch('BATCH',M.supacc,'exact');
%control_supid=strmatch('CTRL',M.supacc,'exact');  %no longer used
BN = M.supdat(batch_supid,:);  %batch numbers for all samples
nbatch=max(BN);
h=histc(M.supdat(batch_supid,:),1:nbatch);
hlarge=find(h>=min_sz);

M = check_array_includes(M,'batch');
include_supid = strmatch('inc_batch',M.supacc,'exact');

%% Loop through batches returning p-values using differential_analysis.m

P=ones(size(M.dat,1),nbatch);   %Initialize P-value array
MN = zeros(size(M.dat,1),nbatch);
b=zeros(length(hlarge),size(M.dat,2));

for i=1:length(hlarge)
      
  if ~isempty(include_supid)  %if array list file specified include samples, only use those in differential analysis
      thisbatch = intersect(find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(include_supid,:)));
      otherbatches = find(M.supdat(include_supid,:));
  else
      thisbatch = find(M.supdat(batch_supid,:)==hlarge(i));
      otherbatches = find(M.supdat(batch_supid,:)~=hlarge(i));
  end
  
  verbose(['(large) BATCH: ' num2str(i)])
  [p]=differential_analysis(M,thisbatch,otherbatches,struct('method','ttest_minvar','minvar',(0.4)^2),0);
                    
  P(:,hlarge(i))=p;  %P gives p-values down the snps by batch (P is #snps X #batches)
  MN(:,hlarge(i)) = mean(M.dat(:,thisbatch),2); %#ok 
  %MN gives the mean for each snp for each batch
  
  
  if ~isempty(include_supid)
      b(i,intersect(find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(include_supid,:))))=1;  
      bghosts(i,find(M.supdat(batch_supid,:)==hlarge(i)))=1; %#ok  %These are all the samples in the batch (even ones who will be corrected but won't be included in the analysis)
  else
      b(i,find(M.supdat(batch_supid,:)==hlarge(i)))=1; %#ok  % b is #batches X #samples (1 if sample in batch, 0 else)           
      bghosts = b;%                                                       %                                                     
  end                                               % 
  
  
end


%% Find snps and batches that need correcting

T = sparse(P <= absolute_pv.*ones(size(P)));  % Threshold matrix T is 1/0 where 1 means p-value was below threshold






verbose(['Finding SNPs with P-value < ' num2str(absolute_pv) ' for at least one batch'],10);
 

verbose(['correcting overall ' num2str(length(find(T))) ' SNPs'],10);

verbose('number of SNPs corrected in each batch :',10);
verbose('\t \t %7.0f \t   %7.0f \n',10, [1:nbatch ;full(sum(T,1))]);

%% Correct the snps

%mini is an NX1 vector (N is number of snps).  Each element of mini gives
%the batch number to be corrected for that snp

%correct data in snp chunks
memchunks = getmemchunkdims(M,'dat',2);
chunksize = memchunks(1);

k = 0;

%while k < length(snps)
while k < getsize(M,'dat',1)

    maxidx = min(getsize(M,'dat',1),k+chunksize);
    thisloopsnps = (k+1):maxidx;
    dat = M.dat(thisloopsnps,:);



    snpMN = mean(dat,2);  % the mean across each snp for this group of data

    delMN = full(T(thisloopsnps,:)).*(MN(thisloopsnps,:) - repmat(snpMN,1,size(T,2)));  % the amount the mean needs to change in each batch

    delMNexp = delMN(:,BN); %  amount to add or subtract from each data element
    dat = dat - delMNexp;


    M.dat(thisloopsnps,:) = dat;
    k = k+chunksize;
    % keyboard
end



            
