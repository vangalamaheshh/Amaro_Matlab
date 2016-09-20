function [X,tidx,nidx] = dRanger_core_algorithmXXX(T,N,P)
% dRanger_core_algorithm(T,N,P)
%
% Gaddy Getz and Mike Lawrence 2008-2009

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'minpairs',2);
P = impose_default_value(P,'discovery_window_size',2000);   % should be >4x the largest insert size
P = impose_default_value(P,'minimum_normal_window',2000);   % look for normal evidence over at least this range
P = impose_default_value(P,'minimum_normal_span_frac',0.5); % normal weirdpair must be at least this frac span
P = impose_default_value(P,'build',[]);
P = impose_default_value(P,'refdir',[]);

id_col=1;chr1_col=2;strand1_col=3;start1_col=4;end1_col=5;
chr2_col=6;strand2_col=7;start2_col=8;end2_col=9;
good1_col=10;good2_col=11;fle_col=12;switch_col=13;

Tmask = zeros(size(T,1),1);
tidx = nan(size(T,1),1);
nidx = nan(size(N,1),1);
trans = {};
nr=0;
NC=24;
if ((~isempty(P.build)) &  (~isempty(P.refdir)))
      filename=[P.refdir '/' P.build '_info.txt'];
      if ~exist(filename,'file')
        error(['missing ' filename])  
      end
      CC = load_table(filename);
      CC = trimStruct(CC,strfindk(CC.use,'D'));
      NC=max(CC.num);
end

for i=1:NC
  for j=i:NC
    for hi=1:2
      disp([i j hi]);
      pos=find(T(:,chr1_col)==i & T(:,chr2_col)==j);
      rx=[round(T(pos,[start1_col start2_col])/P.discovery_window_size-hi/2) T(pos,[ strand1_col strand2_col]) ];
      
      if exist('N','var') && ~isempty(N)
        posN=find(N(:,chr1_col)==i & N(:,chr2_col)==j);
      end
      
      [urx,ui,uj]=unique(rx,'rows');
      hc=histc(uj,1:max(uj));
      [shc,shci]=sort(hc,'descend');

      if ~isempty(shc) && shc(1)>=P.minpairs
        nk = length(find(shc>=P.minpairs));
        for k=1:nk
          if ~mod(k,100), fprintf('i=%d j=%d hi=%d k=%d/%d\n',i,j,hi,k,nk); end

          pk=find(uj==shci(k));
          if length(find(Tmask(pos(pk))==0))<P.minpairs % if less than P.minpairs continue to next peak
            continue;   % used to be "break;"  --> fixed bug 2009-05-27 ML
          end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if abs(T(pos(pk(1)),4)-57.4785e6)>2000, continue; end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          nr=nr+1;
          cur_trans=[];
          cur_trans.end1.chr=i;
          cur_trans.end1.strand=T(pos(pk(1)),strand1_col);
          cur_trans.end2.chr=j;
          cur_trans.end2.strand=T(pos(pk(1)),strand2_col);
          cur_trans.support1=pos(pk);
          
          % first guess based on current set of pairs
          cur_trans.end1.min = min(T(pos(pk),start1_col),[],1);
          cur_trans.end1.max = max(T(pos(pk),end1_col),[],1);
          if cur_trans.end1.strand==0
            cur_trans.end1.bp1=cur_trans.end1.max;
          else
            cur_trans.end1.bp1=cur_trans.end1.min;
          end
          
          cur_trans.end2.min = min(T(pos(pk),start2_col),[],1);
          cur_trans.end2.max = max(T(pos(pk),end2_col),[],1);
          if cur_trans.end2.strand==0
            cur_trans.end2.bp1=cur_trans.end2.max;
          else
            cur_trans.end2.bp1=cur_trans.end2.min;
          end
          
          % look for better definition of the boundary
          % (also, find evidence of a larger field of supporting pairs, i.e. homology region)

          idx = pos;
          idx = idx(T(idx,strand1_col)==cur_trans.end1.strand);
          idx = idx(T(idx,strand2_col)==cur_trans.end2.strand);
          idx = idx(Tmask(idx)==0);
          rad = round(P.discovery_window_size/2);
          if cur_trans.end1.strand==0
            idx = idx(T(idx,end1_col)>cur_trans.end1.bp1-rad);
            idx = idx(T(idx,end1_col)<cur_trans.end1.bp1+rad);
          else
            idx = idx(T(idx,start1_col)>cur_trans.end1.bp1-rad);
            idx = idx(T(idx,start1_col)<cur_trans.end1.bp1+rad);
          end
          if cur_trans.end2.strand==0
            idx = idx(T(idx,end2_col)>cur_trans.end2.bp1-rad);
            idx = idx(T(idx,end2_col)<cur_trans.end2.bp1+rad);
          else
            idx = idx(T(idx,start2_col)>cur_trans.end2.bp1-rad);
            idx = idx(T(idx,start2_col)<cur_trans.end2.bp1+rad);
          end
          if ~isempty(idx)
            allsup = unique([pos(pk);idx]);
          else
            allsup = pos(pk);
          end

          cur_trans.support2 = allsup;
          Tmask(allsup)=1;
          tidx(allsup)=nr;

          cur_trans.end1.min = min(T(allsup,start1_col),[],1);
          cur_trans.end1.max = max(T(allsup,end1_col),[],1);
          if cur_trans.end1.strand==0
            cur_trans.end1.bp2=cur_trans.end1.max;
          else
            cur_trans.end1.bp2=cur_trans.end1.min;
          end
          
          cur_trans.end2.min = min(T(allsup,start2_col),[],1);
          cur_trans.end2.max = max(T(allsup,end2_col),[],1);
          if cur_trans.end2.strand==0
            cur_trans.end2.bp2=cur_trans.end2.max;
          else
            cur_trans.end2.bp2=cur_trans.end2.min;
          end

% cur_trans.end1,cur_trans.end2,keyboard   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if (abs(cur_trans.end1.bp2-16)<1000 && abs(cur_trans.end2.bp2-46944206)<1000) keyboard; end  %%%%%%%%%%%%%%
           
          % look for evidence in normal
          if exist('N','var') && ~isempty(N)
            idx = posN;
            % match start1-end1
            if (cur_trans.end1.max-cur_trans.end1.min+1)<P.minimum_normal_window
              win_cen = (cur_trans.end1.max+cur_trans.end1.min)/2;
              win_st = round(win_cen-P.minimum_normal_window/2);
              win_en = win_st + P.minimum_normal_window;
            else
              win_st = cur_trans.end1.min;
              win_en = cur_trans.end1.max;
            end
            idx = idx(N(idx,start1_col)>=win_st);
            idx = idx(N(idx,end1_col)<=win_en);
            % match start2-end2
            if (cur_trans.end2.max-cur_trans.end2.min+1)<P.minimum_normal_window
              win_cen = (cur_trans.end2.max+cur_trans.end2.min)/2;
              win_st = round(win_cen-P.minimum_normal_window/2);
              win_en = win_st + P.minimum_normal_window;
            else
              win_st = cur_trans.end2.min;
              win_en = cur_trans.end2.max;
            end
            idx = idx(N(idx,start2_col)>=win_st);
            idx = idx(N(idx,end2_col)<=win_en);
            % match requirement for minimum fraction of span
            %    to prevent edge-of-distribution pairs from disqualifying authentic local events
            %    this refinement is due to the SPATS2 example in PR-2832 @ chr12:48,195,500 / 48,198,200
            %    where the matched normal has a ~600-bp-insert pair @ chr12:48,192,700
            if P.minimum_normal_span_frac>0
              span = abs(N(idx,end2_col) - N(idx,end1_col));
              tspan = abs(cur_trans.end2.min - cur_trans.end1.min);
              idx = idx(span >= tspan * P.minimum_normal_span_frac);
            end
            % done matching normal
            cur_trans.normal_support = idx;
            nidx(idx)=nr;
            % check how many match strand1 and strand2
%            idx = idx(N(idx,strand1_col)==cur_trans.end1.strand);
%            idx = idx(N(idx,strand2_col)==cur_trans.end2.strand);
%            cur_trans.normal_support_strict = idx;
          else
            cur_trans.normal_support = [];
          end

          trans{end+1}=cur_trans;

        end % next k
      end % endif
    end % next hi
  end % next j
end % next i

% reformat
X=[];
z = nan(length(trans),1);
X.num=z;
X.chr1=z;X.str1=z;X.pos1=z;
X.chr2=z;X.str2=z;X.pos2=z;
X.span=z;X.tumreads=z;X.normreads=z;%X.normreads_strict=z;
X.min1=z;X.max1=z;X.range1=z;X.stdev1=z;
X.min2=z;X.max2=z;X.range2=z;X.stdev2=z;
for i=1:length(trans), if ~mod(i,1000), fprintf('%d/%d ',i,length(trans)); end
  X.chr1(i) = trans{i}.end1.chr;
  X.str1(i) = trans{i}.end1.strand;
  X.pos1(i) = trans{i}.end1.bp2;
  X.chr2(i) = trans{i}.end2.chr;
  X.str2(i) = trans{i}.end2.strand;
  X.pos2(i) = trans{i}.end2.bp2;
  if (X.chr1(i)==X.chr2(i)), X.span(i) = X.pos2(i)-X.pos1(i); end
  X.tumreads(i)=length(trans{i}.support2);
  X.normreads(i)=length(trans{i}.normal_support);
%  X.normreads_strict(i)=length(trans{i}.normal_support_strict);
  X.min1(i) = trans{i}.end1.min;
  X.max1(i) = trans{i}.end1.max;
  X.range1(i) = trans{i}.end1.max-trans{i}.end1.min+1;
  X.stdev1(i) = round(std(T(trans{i}.support2,start1_col)));
  X.min2(i) = trans{i}.end2.min;
  X.max2(i) = trans{i}.end2.max;
  X.range2(i) = trans{i}.end2.max-trans{i}.end2.min+1;
  X.stdev2(i) = round(std(T(trans{i}.support2,start2_col)));
end, fprintf('\n');

% get rid of candidate rearrangements with "negative span"
% (this is a problem first encountered when trying to run dRanger on capture data)
idx = find(X.span<0);
keep = setdiff(1:slength(X),idx);
a = nan(slength(X),1);
X = reorder_struct(X,keep);
a(keep) = 1:slength(X);
a = as_column(a);
tidx = nansub(a,tidx);
nidx = nansub(a,nidx);

% sort
a = nan(slength(X),1);
[X ord] = sort_struct(X,{'normreads','tumreads'},[1 -1]);
a(ord) = 1:slength(X);
a = as_column(a);
tidx = nansub(a,tidx);
nidx = nansub(a,nidx);
X.num = (1:slength(X))';
