function [trans,F,E,filtered_E]=dRanger2(T,N,min_pairs,window_size,DB,bugfix1_flag)
% dRanger2(T,N,min_pairs,window_size,DB,bugfix1_flag)
%
% Gaddy Getz and Mike Lawrence 2008-2009

if ~exist('bugfix1_flag','var'), bugfix1_flag = true; end

id_col=1;chr1_col=2;strand1_col=3;start1_col=4;end1_col=5;
chr2_col=6;strand2_col=7;start2_col=8;end2_col=9;
good1_col=10;good2_col=11;fle_col=12;switch_col=13;

F=[];
E=[];
filtered_E=[];
Tmask=zeros(size(T,1),1);
trans={};
for i=1:24
  for j=i:24
    for hi=1:2
      disp([i j hi]);
      pos=find(T(:,chr1_col)==i & T(:,chr2_col)==j);
      rx=[round(T(pos,[start1_col start2_col])/window_size-hi/2) T(pos,[ strand1_col strand2_col]) ];
      
      if exist('N','var') && ~isempty(N)
        posN=find(N(:,chr1_col)==i & N(:,chr2_col)==j);
        rn=[round(N(posN,[start1_col start2_col])/window_size-hi/2) N(posN,[ strand1_col strand2_col]) ];
      end
      
      [urx,ui,uj]=unique(rx,'rows');
      hc=histc(uj,1:max(uj));
      [shc,shci]=sort(hc,'descend');

      if ~isempty(shc) && shc(1)>=min_pairs
        for k=1:length(find(shc>=min_pairs))
          pk=find(uj==shci(k));
          if length(find(Tmask(pos(pk))==0))<min_pairs % if less than min_pairs continue to next peak
            continue;   % used to be "break;"  --> fixed bug 2009-05-27 ML
          end
          
          cur_trans=[];
          cur_trans.end1.chr=i;
          cur_trans.end1.strand=T(pos(pk(1)),strand1_col);
          cur_trans.end2.chr=j;
          cur_trans.end2.strand=T(pos(pk(1)),strand2_col);
          cur_trans.support1=pos(pk);
          
          % first guess based on current set of pairs
          if cur_trans.end1.strand==0
            cur_trans.end1.bp1=max(T(pos(pk),end1_col),[],1);
          else
            cur_trans.end1.bp1=min(T(pos(pk),start1_col),[],1);
          end
          
          if cur_trans.end2.strand==0
            cur_trans.end2.bp1=max(T(pos(pk),end2_col),[],1);
          else
            cur_trans.end2.bp1=min(T(pos(pk),start2_col),[],1);
          end
          
          % look for better defnition of the boundary
         
          fprintf('i=%d j=%d hi=%d k=%d/%d chr1=%d strand1=%d pos1=%d chr2=%d strand2=%d pos2=%d\n',...
                  i,j,hi,k,length(find(shc>=min_pairs)),...
                           cur_trans.end1.chr,cur_trans.end1.strand,cur_trans.end1.bp1,...
                           cur_trans.end2.chr,cur_trans.end2.strand,cur_trans.end2.bp1);

          r2 = window_size/2; if bugfix1_flag, r1 = r2; else r1 = 0; end    % bugfix 2009-05-26 ML

          if cur_trans.end1.strand==0 && cur_trans.end2.strand==0
            additional_pk=find(T(pos,end1_col)>cur_trans.end1.bp1-r1 & T(pos,end1_col)<cur_trans.end1.bp1+r2 &...
                               T(pos,strand1_col)==cur_trans.end1.strand  &...
                               T(pos,end2_col)>cur_trans.end2.bp1-r1 & T(pos,end2_col)<cur_trans.end2.bp1+r2 &...
                               T(pos,strand2_col)==cur_trans.end2.strand & Tmask(pos)==0);
          elseif cur_trans.end1.strand==0 && cur_trans.end2.strand==1
            additional_pk=find(T(pos,end1_col)>cur_trans.end1.bp1-r1 & T(pos,end1_col)<=cur_trans.end1.bp1+r2 &...
                               T(pos,strand1_col)==cur_trans.end1.strand  &...
                               T(pos,start2_col)<cur_trans.end2.bp1+r1 & T(pos,start2_col)>cur_trans.end2.bp1-r2 &...
                               T(pos,strand2_col)==cur_trans.end2.strand & Tmask(pos)==0);
          elseif cur_trans.end1.strand==1 && cur_trans.end2.strand==0
            additional_pk=find(T(pos,start1_col)<cur_trans.end1.bp1+r1 & T(pos,start1_col)>cur_trans.end1.bp1-r2 &...
                               T(pos,strand1_col)==cur_trans.end1.strand  &...
                               T(pos,end2_col)>cur_trans.end2.bp1-r1 & T(pos,end2_col)<=cur_trans.end2.bp1+r2 &...
                               T(pos,strand2_col)==cur_trans.end2.strand & Tmask(pos)==0);
          else % 1, 1
            additional_pk=find(T(pos,start1_col)<cur_trans.end1.bp1+r1 & T(pos,start1_col)>cur_trans.end1.bp1-r2 &...
                               T(pos,strand1_col)==cur_trans.end1.strand  &...
                               T(pos,start2_col)<cur_trans.end2.bp1+r1 & T(pos,start2_col)>cur_trans.end2.bp1-r2 &...
                               T(pos,strand2_col)==cur_trans.end2.strand & Tmask(pos)==0);
          end
   
          if ~isempty(additional_pk)
            pk=unique([pk; additional_pk]);
          end
          cur_trans.support2=pos(pk);
          
          if cur_trans.end1.strand==0
            cur_trans.end1.bp2=max(T(pos(pk),end1_col),[],1);
          else
            cur_trans.end1.bp2=min(T(pos(pk),start1_col),[],1);
          end
          
          if cur_trans.end2.strand==0
            cur_trans.end2.bp2=max(T(pos(pk),end2_col),[],1);
          else
            cur_trans.end2.bp2=min(T(pos(pk),start2_col),[],1);
          end          
          Tmask(pos(pk))=1;
          
          
          %% look for related pairs
          %                               T(pos,end1_col)<=cur_trans.end1.bp2 & T(pos,end1_col)>cur_trans.end1.bp2-window_size & ...
          related_pk{1,1}=find(abs(T(pos,end1_col)-cur_trans.end1.bp2)<window_size & ...
                               T(pos,strand1_col)==0  & ...
                               abs(T(pos,end2_col)-cur_trans.end2.bp2)<window_size & ...
                               T(pos,strand2_col)==0);
          related_pk{1,2}=find(...
%abs(T(pos,end1_col)-window_size) & ...
                               abs(T(pos,end1_col)-cur_trans.end1.bp2)<window_size & ...
                               T(pos,strand1_col)==0  & ...
                               abs(T(pos,start2_col)-cur_trans.end2.bp2)<window_size & ...
                               T(pos,strand2_col)==1);
          related_pk{2,1}=find(abs(T(pos,start1_col)-cur_trans.end1.bp2)<window_size & ...
                               T(pos,strand1_col)==1  & ...
                               abs(T(pos,end2_col)-cur_trans.end2.bp2)<window_size & ...
                               T(pos,strand2_col)==0);
          related_pk{2,2}=find(abs(T(pos,start1_col)-cur_trans.end1.bp2)<window_size & ...
                               T(pos,strand1_col)==1  & ...
                               abs(T(pos,start2_col)-cur_trans.end2.bp2)<window_size & ...
                               T(pos,strand2_col)==1);
          cur_trans.related_pk=related_pk;
          
          if exist('N','var') && ~isempty(N)
            normal_related_pk{1,1}=find(abs(N(posN,end1_col)-cur_trans.end1.bp2)<window_size & ...
                                        N(posN,strand1_col)==0  & ...
                                        abs(N(posN,end2_col)-cur_trans.end2.bp2)<window_size & ...
                                        N(posN,strand2_col)==0);
            normal_related_pk{1,2}=find(...
%abs(N(posN,end1_col)-window_size) & ...
                                        abs(N(posN,end1_col)-cur_trans.end1.bp2)<window_size & ...   
%                                        abs(N(posN,end1_col)-cur_trans.end1.bp2)-window_size & ...     % was BUG!
%                                       abs(N(posN,end1_col)-cur_trans.end1.bp2)<window_size & ...   <-- should be
                                        N(posN,strand1_col)==0  & ...
                                        abs(N(posN,start2_col)-cur_trans.end2.bp2)<window_size & ...
                                        N(posN,strand2_col)==1);
            normal_related_pk{2,1}=find(abs(N(posN,start1_col)-cur_trans.end1.bp2)<window_size & ...
                                        N(posN,strand1_col)==1  & ...
                                        abs(N(posN,end2_col)-cur_trans.end2.bp2)<window_size & ...
                                        N(posN,strand2_col)==0);
            normal_related_pk{2,2}=find(abs(N(posN,start1_col)-cur_trans.end1.bp2)<window_size & ...
                                        N(posN,strand1_col)==1  & ...
                                        abs(N(posN,start2_col)-cur_trans.end2.bp2)<window_size & ...
                                        N(posN,strand2_col)==1);
            
            cur_trans.normal_related_pk=normal_related_pk;
            cur_trans.normal_support=normal_related_pk{cur_trans.end1.strand+1,cur_trans.end2.strand+1};
          else
            cur_trans.normal_related_pk={[],[];[],[]};
            cur_trans.normal_support=[];
          end

          out_dat=[];
          for kk=1:length(pk)
            % filter based on DB each entry
            dbe=nan(2,length(DB));
            for dbi=1:length(DB)
              tmp=find(T(pos(pk(kk)),end1_col)>=DB{dbi}{i}.st & T(pos(pk(kk)),start1_col)<=DB{dbi}{i}.en,1);
              if ~isempty(tmp)
                dbe(1,dbi)=tmp(1); % first that it found
              end
              tmp=find(T(pos(pk(kk)),end2_col)>=DB{dbi}{j}.st & T(pos(pk(kk)),start2_col)<=DB{dbi}{j}.en,1);
              if ~isempty(tmp)
                dbe(2,dbi)=tmp(1); % first that it found
              end
            end
            out_dat=[ out_dat; i j k shc(k) kk T(pos(pk(kk)),:) NaN dbe(1,:) dbe(2,:)];
          end
          cur_trans.annotated_support=out_dat;
          
          F=[F; out_dat];
          
          tmp=floor(nanmedian(out_dat,1));
          tmp([ 10 14])=NaN;
          tmp(9)=cur_trans.end1.bp2;
          tmp(13)=cur_trans.end2.bp2;
          tmp(19)=length(cur_trans.normal_support);
         
          E=[E; tmp];
          if all(isnan(out_dat(:,(end-length(DB)*2+1):end)))
            filtered_E=[filtered_E; tmp];
          end

          trans{end+1}=cur_trans;

        end % next k
      end % endif

%    fprintf('Reached point 1\n');  keyboard

    end % next hi
  end % next j
end % next i

% sort trans and E

nt = length(trans);
tsupp = nan(nt,1); nsupp = tsupp;
for i=1:nt
  tsupp(i) = length(trans{i}.support2);
  nsupp(i) = length(trans{i}.normal_related_pk{1,1}) +...
             length(trans{i}.normal_related_pk{2,1}) +...
             length(trans{i}.normal_related_pk{1,2}) +...
             length(trans{i}.normal_related_pk{2,2});
end

[tmp ord1] = sort(tsupp,'descend');
[tmp ord2] = sort(nsupp(ord1));
ord = ord1(ord2);

trans = trans(ord);
E = E(ord,:);
