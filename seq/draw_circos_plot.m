function draw_circos_plot(individual,dRangerfile,segfile)
% draw_circos_plot(individual,dRangerfile,segfile)
%
% individual = name for output files
%
% draws a CIRCOS plot: outputs files to current working directory
%
% dRangerfile = should have chr1, pos1, chr2, pos2, score, [somatic_score]
%
% all rearrangements in file are displayed.
% link strengths are drawn proportional to score
%
% if dRangerfile has "somatic_score" AND "score", then both somatic and germline are drawn.
% if dRangerfile has only "score", then only somatic are drawn.
%
% segfile (optional)
%
% Mike Lawrence 2010-01-28

fprintf('draw_circos_plot: %s\n', individual);

fprintf('Loading dRanger data\n');
X = load_struct(dRangerfile);
X = make_numeric(X,{'chr1','chr2','pos1','pos2','score'});
if isfield(X,'somatic_score')
  X = make_numeric(X,'somatic_score');
  X = rename_field(X,'score','germline_score');
else
  X = rename_field(X,'score','somatic_score');
  X.germline_score = nan(slength(X),1);
end
pngname = [individual '.circos.png'];

% load copy-number data if it exists

have_cn = false;
if exist('segfile','var')
  fprintf('Loading copy-number data\n');
  if exist(segfile,'file')
    prep_SegSeq_for_circos(ssr,segfile);
    have_cn = true;
  else
    fprintf('NOT FOUND: %s\n',segfile);
  end
end

% % link strengths
% nreads = [
% somatic inter   1000    inf
% somatic inter   201     1000
% somatic inter   51      200
% somatic inter   21      50
% somatic inter   12      20
% somatic inter   8       11
% somatic inter   4       7
% somatic inter   0       3];




% output germline inter links
germline_inter_file = [individual '.germline_inter.links'];
x='';
for i=1:slength(X)
  if X.normreads(i)==0, continue; end
  if X.chr1(i)==X.chr2(i), continue; end
    x = [x sprintf('link%05d hs%d %d %d\nlink%05d hs%d %d %d\n',...
       i,X.chr1(i),X.min1(i),X.max1(i),i,X.chr2(i),X.min2(i),X.max2(i))];
end
x = regexprep(x,'hs23','hsX'); x = regexprep(x,'hs24','hsY');
save_textfile(x,germline_inter_file);

% output germline intra links
germline_intra_file = [individual '.germline_intra.links'];
x = '';
for i=1:slength(X)
  if X.normreads(i)==0, continue; end
  if X.chr1(i)~=X.chr2(i), continue; end
  x = [x sprintf('link%05d hs%d %d %d\nlink%05d hs%d %d %d\n',...
     i,X.chr1(i),X.min1(i),X.max1(i),i,X.chr2(i),X.min2(i),X.max2(i))];
end
x = regexprep(x,'hs23','hsX'); x = regexprep(x,'hs24','hsY');
save_textfile(x,germline_intra_file);

% output somatic inter links
somatic_inter_file = [individual '.somatic_inter.links'];
x = '';
for i=1:slength(X)
  if X.normreads(i)>0, continue; end
  if X.chr1(i)==X.chr2(i), continue; end
  x = [x sprintf('link%05d hs%d %d %d\nlink%05d hs%d %d %d\n',...
     i,X.chr1(i),X.min1(i),X.max1(i),i,X.chr2(i),X.min2(i),X.max2(i))];
end
x = regexprep(x,'hs23','hsX'); x = regexprep(x,'hs24','hsY');
save_textfile(x,somatic_inter_file);

% output somatic intra links
somatic_intra_file = [individual '.somatic_intra.links'];
x = '';
for i=1:slength(X)
  if X.normreads(i)>0, continue; end
  if X.chr1(i)~=X.chr2(i), continue; end
  x = [x sprintf('link%05d hs%d %d %d\nlink%05d hs%d %d %d\n',...
     i,X.chr1(i),X.min1(i),X.max1(i),i,X.chr2(i),X.min2(i),X.max2(i))];
end
x = regexprep(x,'hs23','hsX'); x = regexprep(x,'hs24','hsY');
save_textfile(x,somatic_intra_file);

% output CIRCOS configuration file

conffile = [individual '.circos.conf'];

confdir = '/xchip/tcga/gbm/analysis/lawrence/dRanger/circos';

if have_cn
  x = load_textfile([confdir '/template_with_cn.conf']);
  x = regexprep(x,'{copynumber}',cnfile);
else
  x = load_textfile([confdir '/template_without_cn.conf']);
end

pngname = [individual '.circos.png'];

if P.omit_germlines
  x = regexprep(x,'{show_germlines\?}','no');
else
  x = regexprep(x,'{show_germlines\?}','yes');
end

x = regexprep(x,'{ideogram}',[confdir '/ideogram.conf']);
x = regexprep(x,'{ticks}',[confdir '/ticks.conf']);
x = regexprep(x,'{colors}',[confdir '/colors.conf']);
x = regexprep(x,'{fonts}',[confdir '/fonts.conf']);
x = regexprep(x,'{karyotype}',[confdir '/karyotype_custom.txt']);
x = regexprep(x,'{germline_inter_links}',germline_inter_file);
x = regexprep(x,'{germline_intra_links}',germline_intra_file);
x = regexprep(x,'{somatic_inter_links}',somatic_inter_file);
x = regexprep(x,'{somatic_intra_links}',somatic_intra_file);
x = regexprep(x,'{germline_inter_color}',P.germline_inter_color);
x = regexprep(x,'{germline_intra_color}',P.germline_intra_color);
x = regexprep(x,'{somatic_inter_color}',P.somatic_inter_color);
x = regexprep(x,'{somatic_intra_color}',P.somatic_intra_color);
x = regexprep(x,'{dir}',outputdir);
x = regexprep(x,'{file}',pngname);

save_textfile(x,conffile);

% run CIRCOS

fprintf('Calling CIRCOS\n');
cmd = ['circos -conf ' conffile];
[result output] = system(cmd);
if result~=0
  disp(output);
  error('CIRCOS failed');
end
