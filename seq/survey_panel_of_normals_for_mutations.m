function [M C] = survey_panel_of_normals_for_mutations(mutfile,bamfiles,outdir,P)
% survey_panel_of_normals_for_mutations(mutfile,bamfiles,outdir,P)
% note: mutfile can be either an M.mut struct itself or the filename of a maf
%       bamfiles can be either a list of filenames or the filename of a list of a filenames
%
% in returned M, "metric" field columns are as follows:
%     1 = NUMBER OF SAMPLES with high level of the alternate allele
%     2 = "       "       "      low  "            "              "
%     3 = NUMBER OF SAMPLES with high level of noise
%     4 = "       "       "      low  "            "
%     5 = total PERCENT OF READS bearing the alternate allele

if ischar(P), error('parameters have changed'); end
if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'force_report',false);
P = impose_default_value(P,'lane_blacklist','none');
P = impose_default_value(P,'basequalcutoff',20);
P = impose_default_value(P,'bad_normal_nmad_threshold',4);
P = impose_default_value(P,'exclude_indel_counts_from_noise_metric',false); % NO LONGER NEEDED
P = impose_default_value(P,'min_tot_nreads', 20);
P = impose_default_value(P,'hi_min_nreads', 10);
P = impose_default_value(P,'hi_min_pct', 10);
P = impose_default_value(P,'lo_min_nreads', 2);
P = impose_default_value(P,'lo_min_pct', 10);
P = impose_default_value(P,'min_pct_normals',0.5);
P = impose_default_value(P,'min_pct_reads',0.2);
if isfield(P,'suppress_bamidx') && isfield(P,'suppress_bamnames'), error('please use only one'); end
P = impose_default_value(P,'suppress_bamidx',[]);

% BAMS 
if exist([outdir '/bamfile_list.txt'],'file')
  % bamfile list is already finalized: do not alter
  bamfiles = load_lines([outdir '/bamfile_list.txt']);
else
  if P.force_report, error('can''t force_report: apparently hasn''t been run yet!'); end
  orig_len = length(bamfiles);
  % preprocess bamfiles list
  if ischar(bamfiles), bamfiles = load_lines(bamfiles); end
  if ~iscell(bamfiles), error('bamfiles should be a list of filenames'); end
  bamfiles(strcmp('.',bamfiles) | strcmp('..',bamfiles) | strcmp('/',bamfiles))=[];
  bamfiles = bamfiles(demand_files(bamfiles));  % keep only the ones that exist
  % unique bamfiles list
  bamfiles = unique(bamfiles);
  if length(bamfiles)<orig_len, fprintf('Keeping %d/%d unique existent BAMs in list\n',length(bamfiles),orig_len); end
end
if isfield(P,'suppress_bamidx'), P.suppress_bamnames = nansub(bamfiles,P.suppress_bamidx); end
if isempty(P.suppress_bamnames), P.suppress_bamnames = {}; end % % because ismember will choke on []
B=[]; B.bam = bamfiles; clear bamfiles;
B.name = regexprep(B.bam,'^.*/([^/]+)$','$1');
B.name = regexprep(B.name,'\.bam$','');
B.name = regexprep(B.name,'-Normal$','');
B.name = regexprep(B.name,'^(TCGA-..-....).*$','$1');

% MUTATIONS 
if ischar(mutfile)
  M = load_struct(mutfile);
elseif isstruct(mutfile)
  M = mutfile;
  clear mutfile
else
  error('unknown format for mutfile');
end
Morig = M;
M = add_simple_fieldnames(M);
try
  M = add_helper_is_fields(M);
catch me
end
demand_fields(M,{'chr','start','ref_allele','newbase'});
M.chr = convert_chr(M.chr);
M = make_numeric(M,{'start'});

bad_idx = find(isnan(M.chr) | isnan(M.start));
if ~isempty(bad_idx)
  fprintf('%d/%d positions had NaN for chr and/or start:  ignoring these positions\n',length(bad_idx),slength(M));
  M = reorder_struct_exclude(M,bad_idx);
end

ede(outdir);
if ~exist([outdir '/muts_to_screen.maf'],'file')
  save_struct(M,[outdir '/muts_to_screen.maf']);
else
  fprintf('Mutation file already written.\n');
end

% determine which column of C will be the alt
M.alt_idx = listmap(regexprep(M.newbase,'^(.).*$','$1'),{'A','C','G','T'});
if isfield(M,'classification')
  M.alt_idx(strcmp('DEL',M.classification)) = 5;
  M.alt_idx(strcmp('INS',M.classification)) = 6;
end
if isfield(M,'Variant_Type')
  M.alt_idx(strcmp('DEL',M.Variant_Type)) = 5;
  M.alt_idx(strcmp('INS',M.Variant_Type)) = 6;
end
if isfield(M,'ref_allele') && isfield(M,'newbase')
  M.alt_idx(strcmp('-',M.newbase)) = 5; % del
  M.alt_idx(strcmp('-',M.ref_allele)) = 6; % ins
end

bad = isnan(M.alt_idx);
if all(bad)
  error('Don''t know how to find INS/DEL status');
end
if any(bad)
  fprintf('Don''t know how to find INS/DEL status for the following %d mutations:\n',sum(bad));
  if sum(bad)<20
    look(M,isnan(M.alt_idx));
  else
    fprintf('[too many to show]\n');
  end
  keyboard
  error('Please fix this problem');
end

% submit jobs and wait for them to finish
C = get_calls_at_positions(B.bam,M.chr,M.start,outdir,P.lane_blacklist,P.basequalcutoff,P);

% get rid of "suppressed" BAMs
if ~isempty(P.suppress_bamidx)
  C(:,:,P.suppress_bamidx) = [];
  B = reorder_struct_exclude(B,P.suppress_bamidx);
end

% get rid of BAMs that yielded little/no data
B.sum = squeeze(sum(sum(C,1),2));
bad = find(B.sum<0.05*median(B.sum));
if ~isempty(bad)
  C(:,:,bad) = [];
  B = reorder_struct_exclude(B,bad);
end

% convergently analyze results

max_rounds = 1;
round_num = 0;
while(true)
  round_num = round_num + 1;

  tot = sum(C,2);
  alt = nan(size(C,1),size(C,3)); for i=1:slength(M), alt(i,:) = C(i,M.alt_idx(i),:); end
  pct = 100*bsxfun(@rdivide,C,tot);
  covok = repmat(tot,[1 size(C,2) 1])>=P.min_tot_nreads;
  hi = C>=P.hi_min_nreads & pct>=P.hi_min_pct & covok;
  lo = ~hi & C>=P.lo_min_nreads & pct>=P.lo_min_pct & covok;
  sumHi = squeeze(sum(hi,2));
  sumLo = squeeze(sum(lo,2));
  pctAlt = 100*alt./squeeze(tot);
  altHi = (alt>=P.hi_min_nreads & pctAlt>=P.hi_min_pct);
  altLo = ~altHi & (alt>=P.lo_min_nreads & pctAlt>=P.lo_min_pct);
  hi_speci = altHi & sumHi==1;
  lo_speci = altLo & sumLo==1;
  if P.exclude_indel_counts_from_noise_metric
    % OBSOLETE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This was a temporary measure to combat the problem of inflated INS calls
    % due to the Java bamgrasp class returning only an approximate reading about insertions.
    % Now that bamgrasp has been modified to handle insertions exactly, this option is no longer needed.
    sumHiSNP = squeeze(sum(hi(:,1:4,:),2));
    sumLoSNP = squeeze(sum(lo(:,1:4,:),2));
    altHiSNP = altHi; altHiSNP(M.alt_idx>4)=0;
    altLoSNP = altLo; altLoSNP(M.alt_idx>4)=0;
    hi_noise = sumHiSNP-altHiSNP>=2;
    lo_noise = sumLoSNP-altLoSNP>=2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
    hi_noise = sumHi-altHi>=2;
    lo_noise = sumLo-altLo>=2;
  end
  totPctAlt = 100*sum(alt,2)./sum(squeeze(tot),2);

  if round_num==max_rounds, break; end

  % IDENTIFY problematic normals:
  % get rid of samples with bright columns

  B.metric = [sum(hi_speci,1)' sum(lo_speci,1)' sum(hi_noise,1)' ...
              sum(lo_noise)' squeeze(sum(C(:,5,:),1)) squeeze(sum(C(:,6,:),1))];
  me = median(B.metric,1);
  ma = mad(B.metric,1);
  thresh = me + P.bad_normal_nmad_threshold*ma;
  thresh = max(thresh,2); % don't ever throw out a sample because of a single event
  
  B.outlier = bsxfun(@ge,B.metric,thresh);
  B.bad = any(B.outlier,2);
  if any(B.bad)
    B = reorder_struct(B,~B.bad);
    C = C(:,:,~B.bad);
    continue
  else
    break
  end
end

M.metric = [sum(hi_speci,2) sum(lo_speci,2) sum(hi_noise,2) sum(lo_noise,2) totPctAlt];

% determine which mutations are bad
nb = slength(B);
min_nsamps_hi_speci = ceil(nb*P.min_pct_normals/100);
min_nsamps_lo_speci = ceil(nb*P.min_pct_normals/100);
min_nsamps_hi_noise = ceil(nb*P.min_pct_normals/100);
min_nsamps_lo_noise = ceil(nb*P.min_pct_normals/100);
min_totPctAlt = P.min_pct_reads;
mthresh = [min_nsamps_hi_speci min_nsamps_lo_speci min_nsamps_hi_noise min_nsamps_lo_noise min_totPctAlt];
M.bad = any(bsxfun(@ge,M.metric,mthresh),2);

% make list of reasons for calling mutations bad
explained = false(slength(M),1);
M.reason = repmat({'OK'},slength(M),1);
for i=1:size(M.metric,2)
  idx = find(~explained & M.metric(:,i)>=mthresh(i));
  if i<5
    for j=1:length(idx), k=idx(j);
      if i==1, samps = find(hi_speci(k,:)); txt = 'alt allele strongly present in';
      elseif i==2, samps = find(lo_speci(k,:)); txt = 'alt allele weakly present in';
      elseif i==3, samps = find(hi_noise(k,:)); txt = 'high noise at this locus in';
      elseif i==4, samps = find(lo_noise(k,:)); txt = 'weak noise at this locus in';
      end
      M.reason{k} = [txt ' '];
      if length(samps)<10
        for z=1:length(samps)
          M.reason{k} = [M.reason{k} B.name{samps(z)}];
          if z<length(samps), M.reason{k} = [M.reason{k} ', ']; end
        end
      else
        M.reason{k} = [M.reason{k} num2str(length(samps)) ' normals'];
      end
      explained(k) = true;
    end
  elseif i==5
    for j=1:length(idx), k=idx(j);
      M.reason{k} = sprintf('alt allele present in %0.2f%% of normal reads',M.metric(k,i));
      explained(k)=true;
    end
  end 
end

%%%%
%%%% REPORTING 
%%%%

% print list of fates of recurrent sites
fprintf('\nList of recurrent sites and their fates:\n');
[u ui uj] = unique_combos(M.chr,M.start);
h = histc(uj,1:length(ui));
M.nsites = h(uj);
Mr = reorder_struct(M,M.nsites>=2);
Mr.removal_reason = regexprep(Mr.reason,'^(.{1,80}).*$','$1');
[u ui uj] = unique_combos(Mr.chr,Mr.start,Mr.removal_reason);
Mr = reorder_struct(Mr,ui);
Mr = sort_struct(Mr,{'nsites','chr','start'});
pr(Mr,{'nsites','gene','type','chr','start','Protein_Change','removal_reason'});

% print list of all filtered mutations
fprintf('\nList of all mutations flagged for removal:\n');
Mb = reorder_struct(M,M.bad);
[tmp ord] = sort(sum(Mb.metric,2));
Mb.removal_reason = regexprep(Mb.reason,'^(.{1,80}).*$','$1');
pr(Mb,{'gene','type','chr','start','Protein_Change','removal_reason'},ord);

% print summary of bad mutations removed
fprintf('\nSummary of mutations flagged for removal:\n');
status = nansub({'noncoding';'nonsilent';'silent'},1+M.is_coding+M.is_silent);
[a b c] = xcount(M.gene(M.bad),status(M.bad),1);
xcount(M.gene(M.bad),status(M.bad),1);
if length(c)==3
  fprintf('                    %s  %s  %s\n\n',c{1},c{2},c{3});
end

% look for cancer gene coding mutations removed
cgenes = {'TP53','NRAS','KRAS','PTEN','NF1','RB1','CDKN2A','PIK3CA','PIK3R1','EGFR','SF3B1','MYD88','VHL','PBRM1',...
   'SETD2','APC','BRAF','ATM','ARID1A','FBXW7','IDH1','IDH2','CTNNB1','NFE2L2','KEAP1','SMAD4','STK11','SPOP','HRAS',...
   'NOTCH1','MAP3K1','ALK'};
idx = find(ismember(M.gene,cgenes) & M.bad & M.is_coding & ~M.is_silent);
if ~isempty(idx)
  flds = {'patient','gene','chr','start','type','Protein_Change','dbsnp','reason'};
  [tmp ord] = sort(M.start(idx));
  fprintf('\nWARNING: Cancer gene nonsilent mutations flagged for removal:\n');
  pr(M,flds,idx(ord));
end

% write output files
% (0) M,B to results.mat
% (1) filtered MAF
% (2) list of removed mutations

fprintf('Removed %d/%d mutations on the basis of Panel of Normals with %d members.\n',sum(M.bad),slength(Morig),nb);

Mfilt = reorder_struct(Morig,~M.bad);
Mbad = reorder_struct(Morig,M.bad);

outname0 = [outdir '/results.mat'];
outname1 = [outdir '/PoN_filtered.maf'];
outname2 = [outdir '/bad_mutations_only.maf'];

save(outname0,'M','B');

fprintf('Saving filtered MAF to %s\n',outname1);
save_struct(Mfilt,outname1);
fprintf('Saving blacklisted mutations to %s\n',outname2);
save_struct(Mbad,outname2);

if nargout==0
  clear M C
end


