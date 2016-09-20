function x = load_multifind_results_and_sequences(dirname,freeze)
% load_multifind_results_and_sequences(dirname,freeze)
%
% given multifind output directory <dirname>,
% retrieves the alignments.
%
% returns struct x with fields:
%   class (A,A2,B,C,D,E,F)
%   fle (000-999)
%   targ (which bpzoom target was matched)
%   id
%   chr1,start1,end1,strand1,editstring1
%   chr2,start2,end2,strand2,editstring2
%   nerr1,nerr2 -> number of errors (gaps+mismatches)
%   is_tumor, is_normal
%   seq1,seq2 (or blank if none)
%
% freeze (default value = 4)
%
% Mike Lawrence 2009-03-03

if ~exist('freeze','var'), freeze=4; end

x = load_multifind_results(dirname);
nx = slength(x);
lanes = load_lanelist;

x.seq1 = repmat({''},nx,1);
x.seq2 = x.seq1;

[u ui uj] = unique(x.fle);
fprintf('Retrieving sequences... ');
for i=1:length(u)
  fprintf('%d/%d ',i,length(u));
  fle1 = u(i); fle2 = find_pairmate(fle1);
  idx1 = find(lanes.FLE==fle1); idx2 = find(lanes.FLE==fle2);
  if strcmp(lanes.end(idx1),'2'), tmp=idx1; idx1=idx2; idx2=tmp; end
  fname1 = generate_fastb_fname(lanes.FLE(idx1),freeze);
  fname2 = generate_fastb_fname(lanes.FLE(idx2),freeze);
  xidx = find(uj==i);
  if ~strcmpi(fname1,'none')
    try
      s = extract_from_fastb(fname1,x.id(xidx));
      if slength(s) == length(xidx), x.seq1(xidx) = s.seq;
      else fprintf('Problem\n'); end
    catch err, fprintf('Error: %s\n',err.message);
    end
  end
  if ~strcmpi(fname2,'none')
    try
      s = extract_from_fastb(fname2,x.id(xidx));
      if slength(s) == length(xidx), x.seq2(xidx) = s.seq;
      else fprintf('Problem\n'); end
    catch err, fprintf('Error: %s\n',err.message);
    end
  end
end
fprintf('Done.\n');
