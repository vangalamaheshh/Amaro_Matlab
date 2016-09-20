function bpzoom(B,window_size,freeze,temp_dir,output_dir,cohort_range,min_cohort)
% 
%  bpzoom
%
%  Input:
%     B = list of breakpoints.
%         fields: chr (numeric), start, end
%
%     window_size = size of window for matching reads to breakpoints.
%         -->finds reads overlapping with a window of this size centered on each breakpoint
%         -->if window_size==0, then B.start and B.end are taken to specify windows exactly.
%
%     freeze = number of data freeze to use
%
%     temp_dir = directory for writing intermediate results
%
%     output_dir = directory for writing final output
%
%     cohort_range, min_cohort = parameters for filtering "weird" reads
%
%   For each breakpoint:
%     1. Retrieves all pairs with at least one read that was mapped to within <window_size> of the breakpoint.
%     2. Filters out "noise"
%     3. Outputs a "MOKl-o-gram" input file listing all the read pairs.
%
%  Mike Lawrence 2009-02-04

if ~exist('cohort_range','var'), cohort_range = 2000; end;
if ~exist('min_cohort','var'), min_cohort = 4; end;

require_fields(B,{'chr','start','end'});
nb = slength(B);

X=[];
X.chr = B.chr;

if window_size>0
  % check sizes of breakpoints
  B.length = B.end-B.start;
  idx = find(B.length>window_size);
  if ~isempty(idx)
    fprintf('Warning: some breakpoints are larger than the window size:\n');
    look(B,idx);
    fprintf('The query window will not completely contain these breakpoint(s).\n');
    fprintf('Recommend splitting each of them into two queries, or expanding window size.\n');
  end
  B.pos = round((B.start + B.end) / 2);
  X.start = B.pos - round(window_size/2);
  X.end = B.pos + round(window_size/2);
else
  fprintf('Using B.start and B.end to specify windows exactly\n');
  X.start = B.start;
  X.end = B.end;
end

% CALL MULTIFIND

multifind(X,temp_dir,freeze);

x = load_multifind_results(temp_dir);

% MARK IN-TARGET vs. NOT-IN-TARGET

nx = slength(x);
x.is1intarget = false(nx,1);
x.is2intarget = false(nx,1);
w = round(window_size/2);
for i=1:nx
  if ~mod(i,10000), fprintf('%d/%d ',i,nx); end
  b = x.targ(i);
  if x.chr1(i)==B.chr(b) && x.start1(i)<B.pos(b)+w && x.end1(i)>B.pos(b)-w, x.is1intarget(i)=true; end
  if x.chr2(i)==B.chr(b) && x.start2(i)<B.pos(b)+w && x.end2(i)>B.pos(b)-w, x.is2intarget(i)=true; end
end
fprintf('\n');


if 0
% print tallies
xcount([148;x.targ(x.is_normal)],[{'X'};x.class(x.is_normal)]);
xcount([148;x.targ(x.is_tumor)],[{'X'};x.class(x.is_tumor)]);

save('x20090204-2','x');

% survey variety
nb = slength(B);
xD = reorder_struct(x,strcmp('D',x.class));
for b=1:nb
  idx = find(xD.targ==b);
  fprintf('Breakpoint %d\n',b);
  xcount(xD.chr1(idx),xD.chr2(idx));
end
end

% FILTERING RULES
%
%    --> A2/B/C/D events with fewer than min_cohort events within cohort_range are filtered out

xN = reorder_struct(x,grep('^(A|E|F)$',x.class,1));      % Normal
xW = reorder_struct(x,grep('^(A2|B|C|D)$',x.class,1));   % Weird: 1588 pairs
nW = slength(xW);
tum = xW.is_tumor;
bp = xW.targ;
chr = xW.chr1; st = xW.start1; en = xW.end1;
idx = find(~xW.is2intarget);
chr(idx) = xW.chr2(idx); st(idx) = xW.start1(idx); en(idx) = xW.end1(idx);
pos = round((st+en)/2);

win1 = round(pos/cohort_range);
win2 = round((pos/cohort_range)+0.5);
cohort_size = nan(nW,1);
for i=1:nW, cohort_size(i) = sum(tum==tum(i)&bp==bp(i)&chr==chr(i)&(win1==win1(i)|win2==win2(i))); end

xWf = reorder_struct(xW,cohort_size>=min_cohort);    % 366 pairs remain (min_cohort=4)
xf = concat_structs({xN,xWf});

if 0
xcount([148;xf.targ(xf.is_normal)],[{'X'};xf.class(xf.is_normal)]);
xcount([148;xf.targ(xf.is_tumor)],[{'X'};xf.class(xf.is_tumor)]);

% survey variety
nb = slength(B);
xD = reorder_struct(xf,strcmp('D',xf.class));
for b=1:nb
  idx = find(xD.targ==b);
  fprintf('Breakpoint %d\n',b);
  xcount(xD.chr1(idx),xD.chr2(idx));
end
end

%
% CREATE MKLGRAM INPUTFILES
% one for each breakpoint

if ~exist(output_dir,'dir'), mkdir(output_dir), end;
for b=1:nb
  for t={'tum','norm'}
    output_file = sprintf('%s/bp%03d_%s.txt',output_dir,b,t{1});
    y = reorder_struct(xf,xf.targ==b & xf.is_tumor==strcmp(t{1},'tum'));
    ny = slength(y);
    z = [];
    z.pair_no = repmat((1:ny)',2,1);
    z.pairmate = [ones(ny,1);2*ones(ny,1)];
    z.chrA = [y.chr1; y.chr2];
    z.chrA_strand = [y.strand1; y.strand2];
    z.chrA_start = [y.start1; y.start2];
    z.chrA_end = [y.end1; y.end2];
    z.chrB = -ones(ny*2,1); z.chrB_strand = z.chrB;
    z.chrB_start = z.chrB; z.chrB_end = z.chrB;
    fprintf('Outputting %s\n',output_file);
    save_struct(z,output_file);
  end
end

fprintf('\nDone!\n');
fprintf('Type "return" to exit.\n');
keyboard

return

% 2009-02-11
% CREATE MKLGRAM INPUTFILES
% one for each candidate somatic translocation
% TL1=1+2, TL2=3+4, etc.

if ~exist(output_dir,'dir'), mkdir(output_dir), end;
for b=2:2:nb
  for t={'tum','norm'}
    output_file = sprintf('%s/TL%03d_%s.txt',output_dir,b/2,t{1});
    y = reorder_struct(xf,(xf.targ==b-1 | xf.targ==b) & xf.is_tumor==strcmp(t{1},'tum'));
    ny = slength(y);
    z = [];
    z.pair_no = repmat((1:ny)',2,1);
    z.pairmate = [ones(ny,1);2*ones(ny,1)];
    z.chrA = [y.chr1; y.chr2];
    z.chrA_strand = [y.strand1; y.strand2];
    z.chrA_start = [y.start1; y.start2];
    z.chrA_end = [y.end1; y.end2];
    z.chrB = -ones(ny*2,1); z.chrB_strand = z.chrB;
    z.chrB_start = z.chrB; z.chrB_end = z.chrB;
    fprintf('Outputting %s\n',output_file);
    save_struct(z,output_file);
  end
end
