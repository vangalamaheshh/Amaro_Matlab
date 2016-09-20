function R = dRanger_guess_junction(dRanger_results_file,bridgepairs_dir)
% to be renamed ChooseTags

read_length = 76;
insert_length = 450;

if ~exist(bridgepairs_dir,'dir') error('directory %s does not exist',bridgepairs_dir); end

R = load_struct(dRanger_results_file);
R = reorder_struct(R,strcmp(R.filterB,'0'));
R = reorder_struct(R,strcmp(R.filterHCL,'0'));
R = make_numeric(R,{'chr1','chr2','pos1','pos2','min1','max1','min2','max2','str1','str2'});
nr = slength(R);

for i=1:nr
  X1 = load_15column([bridgepairs_dir '/reg' num2str(i) '.sj.txt']);
  X2 = load_15column([bridgepairs_dir '/reg' num2str(i+nr) '.sj.txt']);
  X = concat_structs({X1,X2});
  [tmp ord] = unique([X.rgrp X.namenumber],'rows');
  X = reorder_struct(X,ord);  % remove duplicates
  nx = slength(X);
  Xw = reorder_struct(X, X.chr1==R.chr1(i) & X.chr2==R.chr2(i) &...
     X.start1<R.max1(i) & X.end1>R.min1(i) &...
     X.start2<R.max2(i) & X.end2>R.min2(i));
  nxw = slength(Xw);
  if nxw==0, fprintf('Rearrangement %d: No supporting weirdpairs!\n',i); continue; end

  % (1) determine neighborhoods
  chr1 = R.chr1(i); str1 = R.str1(i);
  chr2 = R.chr2(i); str2 = R.str2(i);
  pos1 = [X.start1(X.chr1==chr1);X.start2(X.chr2==chr1)];
  pos2 = [X.start1(X.chr1==chr2);X.start2(X.chr2==chr2)];
  center1 = round(median(pos1)+read_length/2);
  center2 = round(median(pos2)+read_length/2);
  span = 2000;
  min1 = center1-round(span/2); max1 = min1+span-1;
  min2 = center2-round(span/2); max2 = min2+span-1;

  % Maq-aligned seqs
  seqs = [X.seq1; X.seq2];

  % Merlin-aligned seqs
  F = load_fasta(['/xchip/tcga/gbm/analysis/lawrence/wgs/jump/gg/cp/' R.name{i} '/ClosePairs.both.fasta']);
  seqs = F.seq;

  find_split_reads(seqs,chr1,str1,min1,max1,chr2,str2,min2,max2);
end
