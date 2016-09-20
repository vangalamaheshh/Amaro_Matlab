function H = agglomerate_blast_hits(H,proximity_cutoff,allow_overlaps)
% agglomerates blast hits and marks whether they're shuffled (as in circular hits)
%
% Mike Lawrence 2008-06-11

if ~exist('proximity_cutoff','var')
  proximity_cutoff = 20;
end
if ~exist('allow_overlaps','var')
  allow_overlaps = false;
end

nhits = length(H.hit);

if ~isfield(H, 'is_shuffled')
  H.is_shuffled = false(nhits,1);
end

to_be_discarded = false(nhits,1);
changed = true;
while(changed)
  changed = false;
  for i1=1:nhits
    if to_be_discarded(i1), continue; end
    for i2=1:nhits
      if to_be_discarded(i2), continue; end
      if i1==i2, continue; end
      if ~strcmp(H.query{i1}, H.query{i2}), continue; end
      if ~strcmp(H.hit{i1}, H.hit{i2}), continue; end
      if ~strcmp(H.qstrand{i1}, H.qstrand{i2}), continue; end
      if ~strcmp(H.hstrand{i1}, H.hstrand{i2}), continue; end

      % determine relative positions of hit fragments and query fragments
      CLOSE12 = 1;
      CLOSE21 = 2;
      OVERLAP = 3;
      FAR = 0;

      if H.qstart(i1)<=H.qend(i2) && H.qend(i1)>=H.qstart(i2)
        qgap = 0;
        qrel = OVERLAP;
      elseif H.qend(i1)<H.qstart(i2) && (H.qstart(i2)-H.qend(i1)-1)<=proximity_cutoff
        qgap = (H.qstart(i2)-H.qend(i1)-1);
        qrel = CLOSE12;
      elseif H.qend(i2)<H.qstart(i1) && (H.qstart(i1)-H.qend(i2)-1)<=proximity_cutoff
        qgap = (H.qstart(i1)-H.qend(i2)-1);
        qrel = CLOSE21;
      else
        qrel = FAR;
      end

      if H.hstart(i1)<=H.hend(i2) && H.hend(i1)>=H.hstart(i2)
        hgap = 0;
        hrel = OVERLAP;
      elseif H.hend(i1)<H.hstart(i2) && (H.hstart(i2)-H.hend(i1)-1)<=proximity_cutoff
        hgap = (H.hstart(i2)-H.hend(i1)-1);
        hrel = CLOSE12;
      elseif H.hend(i2)<H.hstart(i1) && (H.hstart(i1)-H.hend(i2)-1)<=proximity_cutoff
        hgap = (H.hstart(i1)-H.hend(i2)-1);
        hrel = CLOSE21;
      else
        hrel = FAR;
      end

      if (qrel==FAR || hrel==FAR), continue; end
      if (qrel==OVERLAP || hrel==OVERLAP) && ~allow_overlaps, continue; end

      % it's a match!

      % first resolve overlaps

      if qrel==OVERLAP
        if H.qstart(i1)<H.qstart(i2), qrel = CLOSE12;
        else qrel = CLOSE21; end
      end
      if hrel==OVERLAP
        if H.hstart(i1)<H.hstart(i2), hrel = CLOSE12;
        else hrel = CLOSE21; end
      end

      % then agglomerate these two hits

      H.qstart(i1) = min(H.qstart(i1),H.qstart(i2));
      H.qend(i1) = max(H.qend(i1),H.qend(i2));
      H.hstart(i1) = min(H.hstart(i1),H.hstart(i2));
      H.hend(i1) = max(H.hend(i1),H.hend(i2));
      gap_length = max(hgap,qgap);
      new_len = H.qend(i1)-H.qstart(i1)+1;
      H.pctid(i1) =100*((H.matchlen(i1)*H.pctid(i1)/100)+(H.matchlen(i2)*H.pctid(i2)/100))/new_len;
      H.matchlen(i1) = new_len;
      H.n_mm(i1) = H.n_mm(i1) + H.n_mm(i2);
      H.n_gaps(i1) = H.n_gaps(i1) + H.n_gaps(i2) + gap_length;
      H.score(i1) = min(H.score(i1),H.score(i2));
      H.E(i1) = max(H.E(i1),H.E(i2));
      if (strcmp(H.qstrand{i1},H.hstrand{i1} && qrel~=hrel)) || ...
         (~strcmp(H.qstrand{i1},H.hstrand{i1} && qrel==hrel)), ...
         % if it's shuffled already
         if H.is_shuffled(i1) || H.is_shuffled(i2)
           % make special note of three-way shuffles!
           fff = fopen('three.txt','at');
           fprintf(fff,'%s %s %d-%d %d-%d\n',H.hit{i1},H.hit{i2},H.hstart(i1),H.hend(i1),H.hstart(i2),H.hend(i2));
           fclose(fff);               
         end
         H.is_shuffled(i1) = true;
      end
      to_be_discarded(i2) = true;
      changed = true;
      break;  % start over

    end % next i2
  end % next i1
end % keep working until no more changes

% now discard the redundant records

to_keep = find(~to_be_discarded);
H = reorder_struct(H, to_keep);

end
