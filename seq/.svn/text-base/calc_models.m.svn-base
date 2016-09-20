function [LL,refref_model,ref,coverage,R,B,ri,bi,refi,LLmu,LLvar]=calc_models(fname,chr,start,stop,priors,params)
%
% fname='/xchip/tcga_scratch/ng/GBM-0188/wgs/bam/normal.bam';
% chr=17;
% start=7417451;
% stop=7617851;
%                ref/ref    ref/alt alt1/alt1 alt1/alt2 hom_N het_N
% priors=log10([ 0.999-1e-8 0.0007  0.0003    1e-8      1/3   2/3]);

if ~exist('params','var'), params=[]; end
params=impose_default_value(params,'filter_mapq0',1);
params=impose_default_value(params,'oldstyle',false);
params=impose_default_value(params,'quiet',0);
params=impose_default_value(params,'trim_offtarget_bases',true);
params=impose_default_value(params,'calc_mu_var',false);
params=impose_default_value(params,'non_ref_sum_q',Inf);
params=impose_default_value(params,'readgroup_blacklist',[]);

tt0 = tic;

params.include_unmapped_reads = false;
if length(fname)>4 && strcmp(fname(end-3:end),'.bam')
  [R B ref ri bi refi] = pull_from_bam(fname,chr,start,stop,params);
elseif length(fname)>5 && strcmp(fname(end-4:end),'.boom')
  [R B ref ri bi refi] = pull_from_boom(fname,chr,start,stop,params);
else
  error('Unknown file format');
end

% remove the +64 from the nonreference bases
nonref = B(:,1)>=63;
B(nonref,1) = B(nonref,1) - 64;

if ~isempty(params.readgroup_blacklist)
  subfprintf('WARNING:  blacklist functionality untested\n');
  subfprintf('Checking reads against lane blacklist\n');
  blacklisted_reads = find(ismember(R(:,1),P.readgroup_blacklist));
  do_not_use_idx=find(ismember(B(:,3),blacklisted_reads));
  subfprintf('   Removed %d blacklisted reads',length(blacklisted_reads));
  B(do_not_use_idx,2)=0; % make base quality = 0
  B(do_not_use_idx,1)=1; % replace with A
end

conv = [];
conv('ACGTN')=1:5;
ref=conv(ref);

if isempty(R) | isempty(B)
  fprintf('No data!\n');
  LL=[]; refref_model = []; coverage = [];
  return
end

subfprintf('Model computation');
subfprintf('\n   Remove dels, Ns, [MAPQ=0''s]...');

if params.filter_mapq0
  do_not_use_idx=find(B(:,1)<0 | B(:,4)<0 | R(B(:,3),8)<1);     % remove dels, Ns, and MAPQ=0
else
  do_not_use_idx=find(B(:,1)<0 | B(:,4)<0);     % remove only dels, Ns
end
B(do_not_use_idx,2)=0; % make base quality = 0
B(do_not_use_idx,1)=1; % replace with A

if ~isinf(params.non_ref_sum_q)
  if size(R,2)<13
    fprintf('\n***   R lacks nonrefsumq column (so far only available via boom)');
    fprintf('\n***   Can''t impose non_ref_sum_q filter\n');
  else
    noisyreads = find(R(:,13)>params.non_ref_sum_q); 
    do_not_use_idx=find(ismember(B(:,3),noisyreads));
    B(do_not_use_idx,2)=0; % make base quality = 0
    B(do_not_use_idx,1)=1; % replace with A
    subfprintf(['\n   Remove noisy reads [' num2str(length(noisyreads)) ']']);
  end

%  do_not_use=zeros(size(B,1),1);
%  for reg=1:length(start)
%    if reg==1
%      firstbase = 1;
%      firstref = 1;
%    else
%      firstbase = bi(reg-1)+1;
%      firstref = refi(reg-1)+1;
%    end
%    lastbase = bi(reg);
%    lastref = refi(reg);
%
%    idx=firstbase:lastbase;
%    not_ref=(B(idx,1)~=ref(B(idx,4)));
%    tmp1=sparse(B(idx,4),B(idx,3),idx);
%    tmp=sparse(B(idx,4),B(idx,3),B(idx,2).*not_ref);
%    tot_q=full(sum(tmp,1));
%    noisy_reads=find(tot_q>params.non_ref_sum_q);
%    tmp1=tmp1(:,noisyreads);
    
  %do_not_use_idx=full(tmp1(:));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  added 2009/10/08 ML
%%%%%%%%  trim bases that lie outside the targeted region
%%%%%%%%  --> yields huge time and memory savings, especially for N phase of somcall

if params.trim_offtarget_bases
  keep = false(size(B,1),1);
  pos1 = 1;
  for i=1:length(start)
    pos2 = bi(i);
    in_bounds = find(B(pos1:pos2,4)>=start(i) & B(pos1:pos2,4)<=stop(i));
    keep(pos1-1+in_bounds) = true;
    pos1 = pos2+1;
  end

  % might as well also get rid of do-not-use bases and MAPQ=0 bases
  keep(B(:,2)==0) = false;

  origlen = size(B,1);
  % B_orig = B;
  % bi_orig = bi;

  B = B(keep,:);
  kcs = cumsum(keep);
  bi = kcs(bi);

  subfprintf('\n   Kept %d/%d bases',size(B,1),origlen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

subfprintf('\n   Compute LP table...');

uq = 0:50;
LP=zeros(10,4,max(uq)+1); % qualities from 0 to 50
if params.calc_mu_var
  mu_table=zeros(max(uq)+1,1);
  var_table=zeros(max(uq)+1,1);
end

for i=1:length(uq)
  e=10^(-uq(i)/10);
  e=min(e,0.75); % make sure error is not >0.75
% AA AC AG AT CC GG TT CG CT GT
% 1  2  3  4  5  6  7  8  9  10
% Prob(base|model)
%        A                        C                    G                    T
  LP(:,:,uq(i)+1)=log10( ...
      [ (1-e)                    (e/3)                (e/3)                (e/3); ...
        (0.5*(1-e)+0.5*(e/3))    (0.5*(1-e)+0.5*e/3)  (0.5*e/3+0.5*e/3)    (0.5*e/3+0.5*e/3); ...
        (0.5*(1-e)+0.5*(e/3))    (0.5*e/3+0.5*e/3)    (0.5*(1-e)+0.5*e/3)  (0.5*e/3+0.5*e/3); ...
        (0.5*(1-e)+0.5*(e/3))    (0.5*e/3+0.5*e/3)    (0.5*e/3+0.5*e/3)    (0.5*(1-e)+0.5*e/3); ...
        (e/3)                    (1-e)                (e/3)                (e/3); ...
        (e/3)                    (e/3)                (1-e)                (e/3); ...
        (e/3)                    (e/3)                (e/3)                (1-e); ...
        (0.5*(e/3)+0.5*(e/3))    (0.5*(1-e)+0.5*e/3)  (0.5*(1-e)+0.5*e/3)  (0.5*e/3+0.5*e/3); ...
        (0.5*(e/3)+0.5*(e/3))    (0.5*(1-e)+0.5*e/3)  (0.5*e/3+0.5*e/3)    (0.5*(1-e)+0.5*e/3); ...
        (0.5*(e/3)+0.5*(e/3))    (0.5*e/3+0.5*e/3)    (0.5*(1-e)+0.5*e/3)  (0.5*(1-e)+0.5*e/3); ...
      ]);
  
  if params.calc_mu_var
    mu1=(1-e)*log10(1-e)+e*log10(e/3);
    tmp1=(1-e)*(log10(1-e).^2)+e*(log10(e/3).^2);
    var1=tmp1-mu1.^2;
        
    mu_table(i)=mu1;
    var_table(i)=var1;
  end
end

subfprintf('\n   Assign L...');

if params.oldstyle
  I=repmat((1:10)',1,size(B,1))+repmat(((B(:,1)-1)*10+(B(:,2)*4*10))',10,1);
else
  I=bsxfun(@plus,(1:10)',((B(:,1)-1)*10+(B(:,2)*40))');                            % 4x faster
end

L=LP(I);

%fprintf('[PAUSE AFTER L]\n');
%keyboard

%              ref/ref    ref/alt alt1/alt1 alt1/alt2 hom_N het_N
%priors=log10([ 0.999-1e-8 0.0007  0.0003    1e-8      1/3   2/3]);

% AA AC AG AT CC GG TT CG CT GT
RI=[ 1 2 2 2 3 3 3 4 4 4;   ... % A
     3 2 4 4 1 3 3 2 2 4;   ... % C
     3 4 2 4 3 1 3 2 4 2;   ... % G
     3 4 4 2 3 3 1 4 2 2;   ... % T
     5 6 6 6 5 5 5 6 6 6;   ... % N
   ]';

subfprintf('\n   Assign LR...');

RP=priors(RI);

if params.oldstyle
  LR=RP(repmat((1:10)',1,size(ref,2))+repmat((ref-1)*10,10,1));
else
  LR=RP(bsxfun(@plus,(1:10)',(ref-1)*10));
end

%% original code

if (0)

  % method 1 (sparse matrix)
  tic
    offset = min(B(:,4))-1;
    K = sparse(1:size(B,1),B(:,4)-offset,1);
    CL = L*K;
    CL = CL(:,start-offset:stop-offset);
  toc   % 2.95 sec (for me 3.34)
  CL1 = CL;

  % method 2 (sort + cumsum)
  % [correct version]
  tic
  [tmp ord] = sort(B(:,4));
  idx = find(tmp(2:end)~=tmp(1:end-1));
  L2 = cumsum(L(:,ord),2);
  CL = [];
  CL(:,u-offset) = L2(:,[idx;end]) - [zeros(10,1),L2(:,idx)];
  CL = CL(start-offset:stop-offset,:);
  toc   % 3.6 sec
  CL2 = CL;
end


% for each requested region,
% (1) collapse L by position -> CL
% (2) LL = LR + CL
% (3) calculate Mut_type, refref_model, and coverage

% allocate memory for concatenated results

subfprintf('\n   Collapse L; Calculate LL, refref_model, and coverage...');

% VECTORIZED VERSION
%if params.vectorized == true

LL = LR;
refref_model = 5*ones(1,refi(end));    % Ns for no-coverage areas
coverage = zeros(refi(end),1);

if params.calc_mu_var
  LLmu=zeros(refi(end),1);
  LLvar=zeros(refi(end),1);
else
  LLmu=[];
  LLvar=[];
end

for reg=1:length(start)

  if reg==1
    firstbase = 1;
    firstref = 1;
  else
    firstbase = bi(reg-1)+1;
    firstref = refi(reg-1)+1;
  end
  lastbase = bi(reg);
  lastref = refi(reg);

  if firstbase>lastbase, continue; end    % skip regions with no coverage

  % may be able to make even faster by trimming bases outside the requested-region

  % collapse L by position -> CL

  offset = B(firstbase,4)-1;
  Kt = sparse(B(firstbase:lastbase,4)-offset,1:(lastbase-firstbase+1),1);
  Lsub = L(:,firstbase:lastbase);
  CL = (Kt*Lsub')';

  if params.calc_mu_var
    base_mu=mu_table(B(firstbase:lastbase,2)+1); 
    base_var=var_table(B(firstbase:lastbase,2)+1);

    reg_LLmu = Kt*base_mu;
    reg_LLvar = Kt*base_var;
    
    reg_LLmu=reg_LLmu+priors(1)-priors(2); % add priors
  end
  
  % LL = LR + CL

  if start(reg)>offset
    CLstart=start(reg)-offset;
    LLstart=1;
  else
    CLstart=1;
    LLstart=offset-start(reg)+2;
  end
  n_read=size(Kt,1);
  n_asked=lastref-firstref+1;
  if stop(reg)<(offset+n_read)
    CLstop=stop(reg)-offset;
    LLstop=n_asked;
  else
    CLstop=n_read;
    LLstop=offset+n_read-start(reg)+1;
  end

  if (LLstop-LLstart+1) ~= (CLstop-CLstart+1), error('sizes do not match'); end

  % perhaps could be made faster if we check whether range covers the whole matrix
  llo = firstref-1;
  LL(:,llo+LLstart:llo+LLstop)=LL(:,llo+LLstart:llo+LLstop)+CL(:,CLstart:CLstop);
  if params.calc_mu_var
    LLmu(llo+LLstart:llo+LLstop)=reg_LLmu(CLstart:CLstop);
    LLvar(llo+LLstart:llo+LLstop)=reg_LLvar(CLstart:CLstop);
  end

  % extract refref_model and coverage
  refsub = ref(firstref:lastref);
  Mut_type=RI(bsxfun(@plus,(1:10)',(refsub-1)*10));
  [mn_Mut_type,mn_MTi]=min(Mut_type,[],1); 
  refref_model(firstref:lastref)=mn_MTi(1,:);

  coverage=zeros(n_asked,1);
  tmp=sum(Kt,2);
  coverage(llo+LLstart:llo+LLstop)=tmp(CLstart:CLstop);

%disp(reg)
%keyboard

end % next reg
%end


if 0 %params.vectorized==false
% NON-VECTORIZED VERSION
offset = B(1,4)-1;
Kt = sparse(B(:,4)-offset,1:size(B,1),1);
CL = Kt*L';
CL=CL';
LL = LR;
if start>offset
  CLstart=start-offset;
  LLstart=1;
else
  CLstart=1;
  LLstart=offset-start+2;
end
n_read=size(Kt,1);
n_asked=size(LL,2);
if stop<(offset+n_read)
  CLstop=stop-offset;
  LLstop=n_asked;
else
  CLstop=n_read;
  LLstop=offset+n_read-start+1;
end
if (LLstop-LLstart+1) ~=  (CLstop-CLstart+1)
  error('sizes do not match');
end
LL(:,LLstart:LLstop)=LL(:,LLstart:LLstop)+CL(:,CLstart:CLstop);
Mut_type=RI(bsxfun(@plus,(1:10)',(ref-1)*10));
[mn_Mut_type,mn_MTi]=min(Mut_type,[],1);
refref_model=mn_MTi(1,:);
coverage=zeros(n_asked,1);
tmp=sum(Kt,2);
coverage(LLstart:LLstop)=tmp(CLstart:CLstop);
end

subfprintf('\n');

  function subfprintf(str,varargin)
    if ~(params.quiet), fprintf(str,varargin{:}); end
  end     
  
end
  
