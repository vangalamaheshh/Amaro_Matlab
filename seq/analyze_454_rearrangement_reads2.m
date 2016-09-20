function [hit,dir,pct] = analyze_454_rearrangement_reads(X,R,P)
% [hit,dir,pct] = analyze_454_rearrangement_reads(X,R,P)
%
% X = struct of 454 reads from load_and_trim_454_data
% R = dRanger validation design file
%
% returns:
%
% hit = which element of R is the best hit 
% dir = which direction is the match (0=left 1=right)
% pct = pct identity of the alignment to amplicon ri
%
% Mike Lawrence 2010-06-08

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'pause',false);
P=impose_default_value(P,'verbose',false);

require_fields(X,{'seq','len'});
require_fields(R,{'amplicon'});
if ~isfield(R,'leftprimer') || ~isfield(R,'rightprimer')
  fprintf('No primer info supplied: will try matching all amplicons (VERY SLOW!)\n');
  blindmode = true;
else
  blindmode = false;
end

nx = slength(X);
nr = slength(R);

X.seq = upper(X.seq);
R.amplicon = upper(R.amplicon);
if ~blindmode
  R.leftprimer = upper(R.leftprimer);
  R.rightprimer = upper(R.rightprimer);
  nlp = cellfun('length',R.leftprimer);
  nrp = cellfun('length',R.rightprimer);
  mnp = min([nlp;nrp]);
  if min(X.len)<mnp, error('some read lengths are shorter than shortest primer'); end
end

% for each read:
%    (1) make list of amplicons for which this read is an exact match to the first N=mnp bases of its leftprimer or rightprimer
%    (2) for each amplicon, do an alignment of the read to the amplicon, measure the pct identity
%    (3) the chosen amplicon is the one with the highest pct identity

% smith-waterman parameters
scmat = -2*ones(size(nuc44)); for z=1:4, scmat(z,z)=4; end
swparams = {'gapopen',32,'extendgap',1,'alphabet','nt','scoringmatrix',scmat};

hit = nan(nx,1); dir = nan(nx,1); pct = nan(nx,1);

tic
fprintf('Assigning to amplicons...\n');
if blindmode, ss=1; else ss=1000; end
for x=1:nx
  if ~mod(x,ss), q=toc;
    fprintf('%d/%d\t%.2f sec\t%.2f min remaining\t%0.2f%% matched\t%0.2f%% >=50pct\t%0.2f%% >=80pct\t%0.2f%% >=90pct \n',...
            x,nx,q,(nx-x)*((q/60)/x),100*mean(~isnan(hit(1:x-1))),...
            100*mean(pct(1:x-1)>=50),100*mean(pct(1:x-1)>=80),100*mean(pct(1:x-1)>=90));
    if P.pause, keyboard, end;
  end
  % find candidates
  if ~blindmode
    lmatches = find(strncmp(X.seq{x}(1:mnp),R.leftprimer,mnp));
    rmatches = find(strncmp(X.seq{x}(1:mnp),R.rightprimer,mnp));
  else
    lmatches = (1:slength(R))';
    rmatches = (1:slength(R))';
  end
  nl = length(lmatches); nr = length(rmatches);
  % try each left candidate
  if nl>0
    lpct = nan(nl,1);
    for i=1:nl
      amp = R.amplicon{lmatches(i)}; maxsc = sum(amp~='N');
      [a b c] = swalign(X.seq{x},amp,swparams{:}); lpct(i) = 100*sum(b(2,:)=='|')/maxsc;
    end
    if P.verbose
      fprintf('  histogram of leftmatch%\n');
      rr = (0:10:100)';
      [rr histc(lpct,rr)]
    end
    [lpct best] = max(lpct); lhit = lmatches(best);
  else
    lhit = nan; lpct = 0;
  end
  % try each right candidate
  if nr>0
    rpct = nan(nr,1); rcseq = rc(X.seq{x});
    for i=1:nr
      amp = R.amplicon{rmatches(i)}; maxsc = sum(amp~='N');
      [a b c] = swalign(rcseq,amp,swparams{:}); rpct(i) = 100*sum(b(2,:)=='|')/maxsc;
    end
    if P.verbose
      fprintf('  histogram of rightmatch%\n');
      rr = (0:10:100)';
      [rr histc(rpct,rr)]
    end
    [rpct best] = max(rpct); rhit = rmatches(best);
  else
    rhit = nan; rpct = 0;
  end
  % choose best candidate
  if nl==0 & nr==0
    hit(x) = nan; pct(x) = 0; dir(x) = nan;
  elseif lpct>rpct
    hit(x) = lhit; dir(x) = 0; pct(x) = lpct;
  else
    hit(x) = rhit; dir(x) = 1; pct(x) = rpct;
  end
  if P.verbose
    fprintf('  best hit: %d pct match\n', pct(x));
  end
end, fprintf('\n');

