function [M n_R n_B] = germcall(n_bamfile,chr,start,stop,params)
% M = germcall(n_bamfile,chr,start,stop)
%
% columns of M are:
%   (1) chr
%   (2) position (1-based)
%   (3) reference base
%   (4) variant model
%   (5) normal LOD: Best vs. next best
%
%  [filtering statistics:]
%   (6) fraction of overlapping reads that are mapping-quality-zero (in normal)
%   (7) average number of mismatches in overlapping reads (in normal)
%   (8) number of distinct weird-pair destinations among overlapping reads (in normal)
%   (9) estimated SNP realness probability (0=artifact -->  1=real)
%
% print_mutation_list(M) to convert to human-readable table
%
% Gaddy Getz and Mike Lawrence 2009-08


if ~exist('params','var') 
  params=[];
end
params=impose_default_value(params,'lod_type','best_vs_next_best');
params=impose_default_value(params,'use_fixed_lod_cutoff',false);
params=impose_default_value(params,'lod_cutoff',5);  % (no longer used)
params=impose_default_value(params,'oldstyle',false);
params=impose_default_value(params,'quiet',0);
params=impose_default_value(params,'K',[0.2 8 5]);  % filtering: inflection point
params=impose_default_value(params,'Z',[5 5 5]);    % filtering: sigmoidicity

subfprintf('GERMCALL starting at %s\n',datestr(now));
tt0=tic;

het_r=1e-3;
tri_r=1e-5;
%                hom         het              hom         het       hom           het 
%                ref/ref     ref/alt          alt1/alt1   alt1/alt2 hom_N         het_N
n_priors=log10([ 1-3*het_r/2 (het_r-tri_r)/3  het_r/2/3   tri_r/3   (1-het_r)/4   het_r/6 ]); % no priors 
n_params = params;
n_params.filter_mapq0 = 0;
n_params.trim_offtarget_bases = false;

subfprintf('Calc_models for N:');

[n_LL,n_refref_model,ref,n_coverage,n_R,n_B,n_ri,n_bi,n_refi]=calc_models(n_bamfile,chr,start,stop,n_priors,n_params);

tt1=toc(tt0);
ref_N_idx=find(ref>4);
ref_no_N=ref;
ref_no_N(ref>4)=1;

if isempty(n_R) || isempty(n_B) 
  subfprintf('No data!\n');
  M = [];
  return;
end

n_refref_LL=n_LL(n_refref_model+10*(0:(size(n_LL,2)-1)));

[sLL,sLLi]=sort(n_LL,'descend');

switch params.lod_type
 case 'best_vs_next_best',
  n_LOD=sLL(1,:)-sLL(2,:); 
 case 'best_vs_ref',
  n_LOD=sLL(1,:)-n_refref_LL; 
 case 'other_vs_ref'
  mxLOD=max(n_LL,[],1);
  tmp=n_LL-repmat(mxLOD,10,1);
  tmp(n_refref_model+10*(0:(size(tmp,2)-1)))=-Inf;
  tmp=sum(10.^(tmp),1);
  n_LOD=log10(tmp)-(n_refref_LL-mxLOD);
 otherwise
  error('no such lod type');
end

n_LOD(ref_N_idx)=-Inf;

% AA AC AG AT CC GG TT CG CT GT
RI=[ 1 2 2 2 3 3 3 4 4 4;   ... % A
     3 2 4 4 1 3 3 2 2 4;   ... % C
     3 4 2 4 3 1 3 2 4 2;   ... % G
     3 4 4 2 3 3 1 4 2 2;   ... % T
     5 6 6 6 5 5 5 6 6 6;   ... % N
   ]';


if params.oldstyle
  mut_type=RI(repmat((1:10)',1,size(ref,2))+repmat((ref-1)*10,10,1));
else
  mut_type=RI(bsxfun(@plus,(1:10)',(ref-1)*10));
end

[smut_type,sMTi]=sort(mut_type,'ascend');

tt2=toc(tt0);

%% LOD statistics
uq=0:50;
mu_table=zeros(length(uq)+1,1);
var_table=zeros(length(uq)+1,1);

for i=1:length(uq)
  e=10^(-uq(i)/10);
  e=min(e,0.75); % make sure error is not >0.75
  
  mu1=(1-e)*log10(1-e)+e*log10(e/3);
  tmp1=(1-e)*(log10(1-e).^2)+e*(log10(e/3).^2);
  var1=tmp1-mu1.^2;
  
  mu2=log10(0.5*(1-e)+0.5*e/3);
  tmp2=(log10(0.5*(1-e)+0.5*e/3)).^2;
  var2=tmp2-mu2.^2;
  
  mu_table(i)=mu1-mu2;
  var_table(i)=var1+var2;
end

%if length(start)>1, error('not yet finished vectorizing to this point'); end

% allocate n_LOD_cutoff
n_LOD_cutoff = nan(1,n_refi(end));

% VECTORIZED VERSION
for reg = 1:length(start)

  if reg==1
    firstbase = 1;
    firstref = 1;
  else
    firstbase = n_bi(reg-1)+1;
    firstref = n_refi(reg-1)+1;
  end
  lastbase = n_bi(reg);
  lastref = n_refi(reg);

  if firstbase>lastbase, continue; end    % skip regions with no coverage

  offset = n_B(firstbase,4)-1;
  Kt = sparse(n_B(firstbase:lastbase,4)-offset,1:(lastbase-firstbase+1),1);

  n_mu=mu_table(n_B(firstbase:lastbase,2)+1);
  n_var=var_table(n_B(firstbase:lastbase,2)+1);

  n_MU = Kt*n_mu;
  n_VAR = Kt*n_var;

  n_MU=n_MU+n_priors(1)-n_priors(2); % add priors

  sub_n_LOD_cutoff=n_MU+norminv(0.05/3e9,0,1)*sqrt(n_VAR);

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

  n_LOD_cutoff((firstref+LLstart-1):(firstref+LLstop-1)) = sub_n_LOD_cutoff(CLstart:CLstop)';
end

if 0
% NON-VECTORIZED VERSION
offset = n_B(1,4)-1; Kt = sparse(n_B(:,4)-offset,1:size(n_B,1),1);
n_mu=mu_table(n_B(:,2)+1); n_var=var_table(n_B(:,2)+1);
n_MU = Kt*n_mu; n_VAR = Kt*n_var;
n_MU=n_MU+n_priors(1)-n_priors(2); % add priors
n_LOD_cutoff=n_MU+norminv(0.05/3e9,0,1)*sqrt(n_VAR);
if start>offset, CLstart=start-offset; else CLstart=1; end
n_read=size(Kt,1); n_asked=size(n_LL,2);
if stop<(offset+n_read),CLstop=stop-offset;else CLstop=n_read;end
n_LOD_cutoff=n_LOD_cutoff(CLstart:CLstop)';
end

tt3=toc(tt0);


%---------------------------------------------------------------------------------

if params.use_fixed_lod_cutoff
  germline=find(sLLi(1,:)~=sMTi(1,:) & n_LOD>params.lod_cutoff);
%germline=find(sLLi(1,:)~=sMTi(1,:) & (n_LOD>params.lod_cutoff | (n_LOD>params.lod_cutoff_2 & n_coverage'> ...
%                                                  params.coverage_cutoff)) );
else
  germline=find(sLLi(1,:)~=sMTi(1,:) & (n_LOD>-n_LOD_cutoff));
end
% germline=find(n_LOD>-n_LOD_cutoff);

% germline=find(sLLi(1,:)~=sMTi(1,:) & n_LOD>0 & n_LOD<=params.lod_cutoff & n_coverage'>30 );

% determine region and pos of each mutation

[region,pos]=map_to_intervals(germline,start,n_refi);

% FILTERING: STATISTICS

subfprintf('\nFiltering...');

haveboom = ~isempty(grep('\.boom$',n_bamfile,1));

if haveboom   % faster method (pre-computed stats in boom)

  n_stem = [n_bamfile '/chr' num2str(chr) '.'];
  n_fmq0 = get_block([n_stem 'fmapqz'],'byte',somatic_pos-1)/100;
  n_avgnmm = get_block([n_stem 'avgnmm'],'byte',somatic_pos-1)/10;
  n_nuwp = get_block([n_stem 'nuwp'],'byte',somatic_pos-1);

else   % old, slower method (for bam)

  z = nan(length(germline),1);
  n_fmq0 = z;
  n_avgnmm = z;
  n_nuwp = z;
  n_ri2 = [0;n_ri];
  for i=1:length(germline)
    r1 = n_ri2(region(i))+1; r2 = n_ri2(region(i)+1);
    nidx = (r1-1)+find(n_R(r1:r2,4)<=pos(i) & n_R(r1:r2,5)>=pos(i)); 
    n_fmq0(i) = mean(n_R(nidx,8)==0); % fraction MAPQZ=0 at this position
    n_avgnmm(i) = mean(n_R(nidx,7)); % average number of mismatches
    n_nuwp(i) = count_uwp(n_R,nidx,chr,pos(i)); % number of unique weird-pairs destinations
  end

end

% FILTERING: SCORE

% tests=
%    (1)   fraction of overlapping reads with MAPQ=0
%    (2)   average number of mismatches in overlapping reads
%    (3)   number of unique weird-pair destinations among overlapping reads

% scoring=
%    (1) calculate probability that the SNP is real, based on each test
%        -- from sigmoidal curve with two test-specific parameters
%    (2) multiply probabilities from each test

n_S = [n_fmq0 n_avgnmm n_nuwp];

tmp =  bsxfun(@power,bsxfun(@rdivide,n_S,params.K),params.Z);

n_P = 1 - (tmp ./ (1+tmp));
filter = prod(n_P,2);

M = [...
    repmat(chr,length(germline),1),...
    pos,...
    ref(germline)',...
    sLLi(1,germline)',...
    n_LOD(germline)',...
    n_fmq0,...
    n_avgnmm,...
    n_nuwp,...
    filter...
];

tt4=toc(tt0);

subfprintf('\n\nGERMCALL finished at %s\n',datestr(now));
subfprintf('calc_models %.2f   n_LOD %.2f   n_LOD_cutoff %.2f   filter %.2f   total %.2f s\n',tt1,tt2-tt1,tt3-tt2,tt4-tt3,tt4);

return



% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% models_types: ref/ref ref/alt alt1/alt1 alt1/alt2 hom_N het_N
% AA AC AG AT CC GG TT CG CT GT
RI=[ 1 2 2 2 3 3 3 4 4 4;   ... % A
     3 2 4 4 1 3 3 2 2 4;   ... % C
     3 4 2 4 3 1 3 2 4 2;   ... % G
     3 4 4 2 3 3 1 4 2 2;   ... % T
     5 6 6 6 5 5 5 6 6 6;   ... % N
   ]';

model_type=RI(sLLi(1,:)+(ref-1)*10);


%---- visualization


n_C=sparse(n_B(:,4),n_B(:,3),n_B(:,1));
n_Q=sparse(n_B(:,4),n_B(:,3),n_B(:,2));

[tmp,j]=min(n_LOD(germline)); 
%[tmp,j]=max(n_LOD(germline));
models={'AA','AC','AG','AT','CC','GG','TT','CG','CT','GT'};
ref_base={'A','C','G','T','N'};
% j=2;
germline(j)
pos=start+germline(j)-1;
disp(pos);
disp(histc(full(n_C(pos,:)),1:4));
disp(['Ref=' ref_base{ref(germline(j))}]);
for k=1:10
  fprintf(1,'%s\t%f\n',models{k},n_LL(k,germline(j)));
end



% non-ref sites in genomic order
for j=1:length(germline)
  fprintf('chr%d:%d   %s -> %s   LOD = %f\n',...
    chr, start+germline(j)-1, ref_base{ref(germline(j))},...
    models{sLLi(1,germline(j))}, n_LOD(germline(j)));
end

pos = 10946283-start+1;
for k=1:10
  fprintf(1,'%s\t%f\n',models{k},n_LL(k,pos));
end
disp(n_LOD(pos));   % 4.72


figure(1); clf;
cols='brg';
for i=2:4
  idx=find(model_type(germline)==i);
  plot(n_coverage(germline(idx)),n_LOD(germline(idx)),[cols(i-1) '.']); hold on;
end
xlabel('Coverage','FontSize',14);
ylabel('LOD','FontSize',14);
title(['chr' num2str(chr) ':' num2str(start) '-' num2str(stop)],'FontSize',18);

prc=zeros(50,1);
stats=zeros(50,3);
for i=5:50
  cov_idx=find(n_coverage'==i);
  idx=cov_idx(find(model_type(cov_idx)==1));
  prc(i)=length(find(n_LOD(idx)<5))/length(idx); % prctile(n_LOD(idx),99);
  stats(i,1)=mean(n_LOD(cov_idx));
  stats(i,2)=std(n_LOD(cov_idx));
  stats(i,3)=length(find(n_LOD(cov_idx)<5))/length(cov_idx); 
  stats(i,4)=length(find(n_LOD(cov_idx)<5 & model_type(cov_idx)==1))/length(find(n_LOD(cov_idx)<5));
end


figure(1); clf;
errorbar(stats(:,1),stats(:,2));

figure(2); clf;
semilogy(stats(:,3));

figure(3); clf;
plot(stats(:,4));


cov1=n_coverage;
cov1(cov1>50)=50;
A=hist2d(n_LOD(model_type~=1),cov1(model_type~=1)',-10:0.1:20,4:51);
figure(1);clf;
imagesc(log(A+1));

cov1=n_coverage;
cov1(cov1>50)=50;
A=hist2d(n_LOD,cov1',-10:0.1:20,4:51);
figure(1);clf;
imagesc(log(A+1));

find(n_coverage'==30 & model_type~=1 & abs(n_LOD-3)<0.4)
find(n_coverage'>30 & model_type~=1 & n_LOD<5 & n_LOD>0)


%
%
%  function nuwp = subfunction_nuwp(R,idx,pos)
%    wp = R(idx,10:12);                                    % chr start strand
%    wp = wp(wp(:,1)>0,:);                                 % remove unpaired and unmapped-pair
%    nw = find(wp(:,1)==chr & abs(wp(:,2)-pos)<10000);     % remove non-weird
%    wp = wp(setdiff(1:size(wp,1),nw),:);
%    wp(:,2) = round(wp(:,2) / 2000);                      % smooth to windows of 2kB
%    nuwp = size(unique(wp,'rows'),1);
%  end
%

  function subfprintf(str,varargin)
    if ~(params.quiet), fprintf(str,varargin{:}); end
  end     
  

end % of main function
