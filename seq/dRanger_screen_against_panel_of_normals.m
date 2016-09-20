function X = dRanger_screen_against_panel_of_normals(X,P)
% dRanger_screen_against_panel_of_normals(X,P)
%
% Mike Lawrence 2010-03-24

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'panel_of_normals_database',[]);
P = impose_default_value(P,'minimum_normal_window',2000);
P = impose_default_value(P,'minimum_normal_span_frac',0.5);

if ~exist(P.panel_of_normals_database,'file')
  fprintf('Could not open panel_of_normals_database\n');
  fprintf('Returning without doing anything\n')'
  return
end

list = load_lines(P.panel_of_normals_database);

X = make_numeric(X,{'chr1','chr2','min1','min2','max1','max2'});

id_col=1;chr1_col=2;strand1_col=3;start1_col=4;end1_col=5;
chr2_col=6;strand2_col=7;start2_col=8;end2_col=9;
good1_col=10;good2_col=11;fle_col=12;switch_col=13;

E = zeros(slength(X),length(list));

fprintf('SCREENING AGAINST PANEL OF NORMALS\n');

used = false(length(list),1);
for n=1:length(list)
  fprintf('\n%d/%d  *************************************\n',n,length(list));
  try
    N = load_dRanger_input(list{n},P);   % (performs dedup)
  catch me
    fprintf('Error loading\n');
    continue;
  end

  % make chromosome index
  idx1 = cell(24,1);
  for i=1:24, idx1{i} = find(N(:,chr1_col)==i); end
  idx2 = cell(24,1);
  for i=1:24, idx2{i} = find(N(:,chr2_col)==i); end
  nidx = cell(24,24);
  for i=1:24, for j=i:24
    nidx{i,j} = intersect(idx1{i},idx2{j});
  end, end

  % test each rearrangement for evidence in this panel normal

  for x=1:slength(X)
    if ~mod(x,1000), fprintf('%d/%d ',x,slength(X)); end
    % match chromosomes
    idx = nidx{X.chr1(x),X.chr2(x)};
    % match start1-end1
    if (X.max1(x)-X.min1(x)+1)<P.minimum_normal_window
      win_cen = (X.max1(x)+X.min1(x))/2;
      win_st = round(win_cen-P.minimum_normal_window/2);
      win_en = win_st + P.minimum_normal_window;
    else
      win_st = X.min1(x);
      win_en = X.max1(x);
    end
    idx = idx(N(idx,start1_col)>=win_st);
    idx = idx(N(idx,end1_col)<=win_en);
    % match start2-end2
    if (X.max2(x)-X.min2(x)+1)<P.minimum_normal_window
      win_cen = (X.max2(x)+X.min2(x))/2;
      win_st = round(win_cen-P.minimum_normal_window/2);
      win_en = win_st + P.minimum_normal_window;
    else
      win_st = X.min2(x);
      win_en = X.max2(x);
    end
    idx = idx(N(idx,start2_col)>=win_st);
    idx = idx(N(idx,end2_col)<=win_en);
    % match requirement for minimum fraction of span
    %    to prevent edge-of-distribution pairs from disqualifying authentic local events
    %    this refinement is due to the SPATS2 example in PR-2832 @ chr12:48,195,500 / 48,198,200
    %    where the matched normal has a ~600-bp-insert pair @ chr12:48,192,700
    if P.minimum_normal_span_frac>0
      span = abs(N(idx,end2_col) - N(idx,end1_col));
      tspan = abs(X.min2(x) - X.min1(x));
      idx = idx(span >= tspan * P.minimum_normal_span_frac);
    end
    % done matching normal
    E(x,n) = length(idx);
  end % next rearrangement

  used(n) = true;

end % next normal in panel

% add fields to X

X.normpanelreads = sum(E,2);
X.normpanelsamps = sum(E>0,2);

names = regexprep(list,'.*/(\S+)-Normal.*','$1');

for i=1:slength(X)
  txt = '';
  z = 0;
  for j=1:length(list)
    if ~used(j), continue; end
    if E(i,j)>0
      if ~isempty(txt), txt = [txt ',']; end
      txt = [txt names{j} ':' num2str(E(i,j))];
    else
      z=z+1;
    end
  end
  if z>0
    if ~isempty(txt), txt = [txt ',']; end
    txt = [txt num2str(z) '_samples:0'];
  end
  X.normpaneldetails{i,1} = txt;
end


