function compare_dRanger_and_SegSeq(sample,tolerance,pvalue_cutoff,bkpt_uncertainty_cutoff)
if ~exist('tolerance','var'), tolerance = 1000; end
if ~exist('pvalue_cutoff','var'), pvalue_cutoff=1e-100; end
if ~exist('bkpt_uncertainty_cutoff','var'), bkpt_uncertainty_cutoff=2000; end

fprintf('dRanger_dRanger_and_SegSeq\n\tsample = %s\n',sample);
fprintf('\ttolerance = %d\n\tpvalue_cutoff = %d\n\tbkpt_uncertainty_cutoff = %d\n',...
   tolerance,pvalue_cutoff,bkpt_uncertainty_cutoff)

fprintf('\nLoading dRanger data\n');
direc = ['/xchip/tcga_scratch/lawrence/' sample];
name2 = upper(regexprep(sample,'/','-'));
fname = [direc '/' name2 '_dRanger_results_bfiltered.txt'];
if ~exist(fname,'file'), error('Can''t find file %s',fname);end
D = load_struct(fname);
D = reorder_struct(D,strcmp(D.filterB,'0'));
D = reorder_struct(D,strcmp(D.filterHCL,'0'));
D = reorder_struct(D,strcmp(D.normreads,'0'));
D = make_numeric(D,{'chr1','chr2','pos1','pos2'});
fprintf('%d filtered somatic rearrangements\n',slength(D));

fprintf('\nLoading SegSeq data\n');
ssr = [direc '/SegSeq_results.txt'];       % may also be SegSeq_results.seg.txt 
if ~exist(ssr,'file'), error('Can''t find file %s',ssr); end
S = load_struct(ssr,[],1,'lowercase_fieldnames');

if isfield(S,'segseqpval')
  S.p = S.segseqpval;
elseif isfield(S,'segseqpvalue')
  S.p = S.segseqpvalue;
elseif isfield(S,'pvalue')
  S.p = S.pvalue;
elseif isfield(S,'pval')
  S.p = S.pval;
elseif isfield(S,'p')
  S.p = S.p;
else
  error('Can''t find p-value column in %s\n',ssr);
end

S = make_numeric(S,{'chromosome','start','end','p'});

% convert to breakpoints

idx = find([0;~diff(S.chromosome)]);

B.chr = S.chromosome(idx);
x1 = S.end(idx-1); x2 = S.start(idx);
B.start = min([x1 x2],[],2);
B.end = max([x1 x2],[],2);
B.p1 = S.p(idx-1);
B.p2 = S.p(idx);

B = reorder_struct(B,B.p1<=pvalue_cutoff | B.p2<=pvalue_cutoff);

idx = find(B.end-B.start>bkpt_uncertainty_cutoff);
if ~isempty(idx)
  fprintf('Removing %d breakpoints with uncertainty greater than %d bp\n',....
    length(idx),bkpt_uncertainty_cutoff);
  B = reorder_struct(B,setdiff(1:slength(B),idx));
end

fprintf('%d segments --> %d filtered copy-number breakpoints\n',...
  slength(S),slength(B));

% compare positions

nb = slength(B);
nd = slength(D);

bchr = repmat(B.chr,1,nd);
bstart = repmat(B.start,1,nd);
bend = repmat(B.end,1,nd);
dpos1 = repmat(D.pos1',nb,1);
dpos2 = repmat(D.pos2',nb,1);
dchr1 = repmat(D.chr1',nb,1);
dchr2 = repmat(D.chr2',nb,1);

hit1 = bstart <= dpos1+tolerance & bend >= dpos1-tolerance;
hit2 = bstart <= dpos2+tolerance & bend >= dpos2-tolerance;
hit1(bchr~=dchr1) = false;
hit2(bchr~=dchr2) = false;

% print report

fprintf('\nOf %d filtered copy-number breakpoints,\n',nb);
fprintf('\t%d have at least one matching rearrangement breakpoint.\n',sum(sum(hit1,2)+sum(hit2,2)>0));
fprintf('\t%d have no matching rearrangement breakpoint.\n',sum(sum(hit1,2)+sum(hit2,2)==0));

fprintf('\nOf %d filtered somatic rearrangements,\n',nd);
fprintf('\t%d have at least one breakpoint matching a copy-number breakpoint\n',...
   sum(sum(hit1,1)>0 | sum(hit2,1)>0));
fprintf('\t%d have both breakpoints matching copy-number breakpoints\n',...
   sum(sum(hit1,1)>0 & sum(hit2,1)>0));
fprintf('\t%d have neither breakpoint matching a copy-number breakpoint\n',...
   sum(sum(hit1,1)==0 & sum(hit2,1)==0));

