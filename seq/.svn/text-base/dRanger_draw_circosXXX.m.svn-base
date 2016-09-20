function dRanger_draw_circosXXX(sample,P)

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'dRanger_results_filespec',[]);
P=impose_default_value(P,'confdir','/xchip/tcga/gbm/analysis/lawrence/dRanger/circos');
P=impose_default_value(P,'conffile',[]);
P=impose_default_value(P,'executable_name','circos-0.49');
P=impose_default_value(P,'perl_modules_library','');
P=impose_default_value(P,'segfile_filespec',[]);
P=impose_default_value(P,'omit_chrY_cn',false);
P=impose_default_value(P,'results_name','dRanger_results');
P=impose_default_value(P,'omit_germlines',true);
P=impose_default_value(P,'filter_by_filter',false);
P=impose_default_value(P,'filter_by_tumreads',false);
P=impose_default_value(P,'tumreads_cutoff',[]);
P=impose_default_value(P,'filter_by_score',false);
P=impose_default_value(P,'score_cutoff',[]);
P=impose_default_value(P,'filter_by_somatic_score',false);
P=impose_default_value(P,'somatic_score_cutoff',[]);
P=impose_default_value(P,'germline_intra_color','vvdgrey');
P=impose_default_value(P,'germline_inter_color','vvdgrey');
P=impose_default_value(P,'somatic_intra_color','green');  % old style = "blue"
P=impose_default_value(P,'somatic_inter_color','purple');  % old sdtyle = "red"
P=impose_default_value(P,'CIRCOS_outname',[]);
P=impose_default_value(P,'CIRCOS_outdir',[]);
P=impose_default_value(P,'ideogram_file',[P.confdir '/ideogram.conf']);
P=impose_default_value(P,'ticks_file',[P.confdir '/ticks.conf']);
P=impose_default_value(P,'colors_file',[P.confdir '/colors.conf']);
P=impose_default_value(P,'fonts_file',[P.confdir '/fonts.conf']);
P=impose_default_value(P,'karyotype_file',[P.confdir '/karyotype_custom.txt']);

fprintf('dRanger_draw_circos\n\tsample = %s\n\n',sample);

fprintf('Loading dRanger data\n');
if ~isempty(P.dRanger_results_filespec)
  fname = P.dRanger_results_filespec;
  direc = '.';
else
  basedir = '/xchip/cga1/lawrence';
  direc = [basedir '/' sample];
  name2 = upper(regexprep(sample,'/','-'));
  fname = [direc '/' P.results_name '.txt'];
end
if ~exist(fname,'file'), error('Can''t find file %s',fname);end
X = load_struct(fname);

% FILTERING

if P.filter_by_filter
  fprintf('Filtered by filter==0\n');
  X = reorder_struct(X,strcmp(X.filter,'0'));
end

if P.filter_by_tumreads
  if isempty(P.tumreads_cutoff)
    error('Must specify P.tumreads_cutoff');
  else
   fprintf('Filtered by tumreads >= %d\n',P.tumreads_cutoff);
    X = make_numeric(X,{'tumreads'});
    X = reorder_struct(X,X.tumreads>=P.tumreads_cutoff);
  end
end

if P.filter_by_score
  if isempty(P.score_cutoff)
    error('Must specify P.score_cutoff');
  else
   fprintf('Filtered by score >= %d\n',P.score_cutoff);
    X = make_numeric(X,{'score'});
    X = reorder_struct(X,X.score>=P.score_cutoff);
  end
end

if P.filter_by_somatic_score
  if isempty(P.somatic_score_cutoff)
    error('Must specify P.somatic_score_cutoff');
  else
   fprintf('Filtered by somatic_score >= %d\n',P.somatic_score_cutoff);
    X = make_numeric(X,{'somatic_score'});
    X = reorder_struct(X,X.somatic_score>=P.somatic_score_cutoff);
  end
end

if ~isfield(X,'min1'), X.min1 = X.pos1; end
if ~isfield(X,'max1'), X.max1 = X.pos1; end
if ~isfield(X,'min2'), X.min2 = X.pos2; end
if ~isfield(X,'max2'), X.max2 = X.pos2; end

if ~isempty(grep('chr',X.chr1)), X.chr1 = convert_chr(X.chr1,P); end
if ~isempty(grep('chr',X.chr2)), X.chr2 = convert_chr(X.chr2,P); end

X = make_numeric(X,{'chr1','chr2','normreads','min1','max1','min2','max2'});

% make sure min1 and max2 are all in range for chromosomes
K = load_lines(P.karyotype_file);
K = grep('^chr - ',K);
% species (hs human, mm mouse)
spec=regexp(K{1},' ','split'); spec=spec{3}(1:2);
fmt=['chr - ' spec '(.*) (.*) (.*) (.\d+) .*'];
k = parse(K,fmt,{spec,'chr','min','max'},3:4);

if length(unique(k.chr))<10, error('Error processing karytype file'); end
k.chr = convert_chr(k.chr,P);
for i=1:slength(X)
  ci = find(X.chr1(i)==k.chr);
  if isempty(ci)
    fprintf('WARNING: chr not found\n');
  else
    X.min1(i) = min(max(X.min1(i),k.min(ci)),k.max(ci));
    X.max1(i) = min(max(X.max1(i),k.min(ci)),k.max(ci));
  end
  ci = find(X.chr2(i)==k.chr);
  if isempty(ci)
    fprintf('WARNING: chr not found\n');
  else
    X.min2(i) = min(max(X.min2(i),k.min(ci)),k.max(ci));
    X.max2(i) = min(max(X.max2(i),k.min(ci)),k.max(ci));
  end
end


if isempty(P.CIRCOS_outdir)
  outputdir = direc;
else
  outputdir = P.CIRCOS_outdir;
  ensure_dir_exists(P.CIRCOS_outdir);
end

if isempty(P.CIRCOS_outname)
  if ~P.omit_germlines
   pngname = [name2 '_circos.png'];
  else
    pngname = [name2 '_circos_no_germlines.png'];
  end
else
  if P.CIRCOS_outname(1)=='/'
    error('P.CIRCOS_outname should not be a full path.  Use P.CIRCOS_outdir for that.');
  end
  pngname = P.CIRCOS_outname;
end

auxfilestem = regexprep([outputdir '/' pngname],'^(.*)\.png$','$1');

% load copy-number data if it exists

if ~isempty(P.segfile_filespec)
  ssr = P.segfile_filespec;
  if ~exist(ssr,'file'), fprintf('Not found: %s\n',ssr); end
else
  ssr = [direc '/SegSeq_results.seg.txt'];
  if ~exist(ssr,'file')
    ssr = [direc '/SegSeq_results.txt'];
  end
end

if exist(ssr,'file')
  fprintf('Loading copy-number data\n');
  cnfile = [auxfilestem '.circos_copynumber.txt'];
  prep_SegSeq_for_circos(ssr,cnfile,P);
  have_cn = true;
else
  have_cn = false;
end

% output germline inter links
germline_inter_file = [auxfilestem '.germline_inter.links'];
x='';

X.chromomsome1 = convert_chr_backXXX(X.chr1,P);
X.chromomsome2 = convert_chr_backXXX(X.chr2,P);

for i=1:slength(X)
  if X.normreads(i)==0, continue; end
  if X.chr1(i)==X.chr2(i), continue; end
    %x = [x sprintf('link%05d %s%d %d %d\nlink%05d %s%d %d %d\n',...
    %   i,spec,X.chr1(i),X.min1(i),X.max1(i),i,spec,X.chr2(i),X.min2(i),X.max2(i))];
    x = [x sprintf('link%05d %s%s %d %d\nlink%05d %s%s %d %d\n',...
       i,spec,X.chromosome1(i),X.min1(i),X.max1(i),i,spec,X.chromosome2(i),X.min2(i),X.max2(i))];
end
%x = regexprep(x,[spec '23'],[spec 'X']); x = regexprep(x,[spec '24'],[spec 'Y']);
x = regexprep(x,[spec 'chr'],spec); 
save_textfile(x,germline_inter_file);

% output germline intra links
germline_intra_file = [auxfilestem '.germline_intra.links'];
x = '';
for i=1:slength(X)
  if X.normreads(i)==0, continue; end
  if X.chr1(i)~=X.chr2(i), continue; end
  %x = [x sprintf('link%05d %s%d %d %d\nlink%05d %s%d %d %d\n',...
  %   i,spec,X.chr1(i),X.min1(i),X.max1(i),i,spec,X.chr2(i),X.min2(i),X.max2(i))];
  x = [x sprintf('link%05d %s%s %d %d\nlink%05d %s%s %d %d\n',...
     i,spec,X.chromosome1(i),X.min1(i),X.max1(i),i,spec,X.chromosome2(i),X.min2(i),X.max2(i))];
end
%x = regexprep(x,[spec '23'],[spec 'X']); x = regexprep(x,[spec '24'],[spec 'Y']);
x = regexprep(x,[spec 'chr'],spec); 

save_textfile(x,germline_intra_file);

% output somatic inter links
somatic_inter_file = [auxfilestem '.somatic_inter.links'];
x = '';
for i=1:slength(X)
  if X.normreads(i)>0, continue; end
  if X.chr1(i)==X.chr2(i), continue; end
  %x = [x sprintf('link%05d %s%d %d %d\nlink%05d %s%d %d %d\n',...
  %   i,spec,X.chr1(i),X.min1(i),X.max1(i),i,spec,X.chr2(i),X.min2(i),X.max2(i))];
  x = [x sprintf('link%05d %s%s %d %d\nlink%05d %s%s %d %d\n',...
     i,spec,X.chromosome1(i),X.min1(i),X.max1(i),i,spec,X.chromosome2(i),X.min2(i),X.max2(i))];
end
%x = regexprep(x,[spec '23'],[spec 'X']); x = regexprep(x,[spec '24'],[spec 'Y']);
x = regexprep(x,[spec 'chr'],spec);
save_textfile(x,somatic_inter_file);

% output somatic intra links
somatic_intra_file = [auxfilestem '.somatic_intra.links'];
x = '';
for i=1:slength(X)
  if X.normreads(i)>0, continue; end
  if X.chr1(i)~=X.chr2(i), continue; end
  %x = [x sprintf('link%05d %s%d %d %d\nlink%05d %s%d %d %d\n',...
  %   i,spec,X.chr1(i),X.min1(i),X.max1(i),i,spec,X.chr2(i),X.min2(i),X.max2(i))];
  x = [x sprintf('link%05d %s%s %d %d\nlink%05d %s%s %d %d\n',...
     i,spec,X.chromosome1(i),X.min1(i),X.max1(i),i,spec,X.chromosome2(i),X.min2(i),X.max2(i))];
end
%x = regexprep(x,[spec '23'],[spec 'X']); x = regexprep(x,[spec '24'],[spec 'Y']);
x = regexprep(x,[spec 'chr'],spec);

save_textfile(x,somatic_intra_file);

% output CIRCOS configuration file

conffile = [auxfilestem '.circos.conf'];

confdir = P.confdir;

if ~isempty(P.conffile)
  fname = P.conffile;
else
  if have_cn
    fname = [confdir '/template_with_cn.conf'];
  else
    fname = [confdir '/template_without_cn.conf'];
  end
end

fprintf('Loading %s\n',fname);
x = load_textfile(fname);
if exist('cnfile','var')
  x = regexprep(x,'{copynumber}',cnfile);
end

if P.omit_germlines
  x = regexprep(x,'{show_germlines\?}','no');
else
  x = regexprep(x,'{show_germlines\?}','yes');
end

x = regexprep(x,'{ideogram}',P.ideogram_file);
x = regexprep(x,'{ticks}',P.ticks_file);
x = regexprep(x,'{colors}',P.colors_file);
x = regexprep(x,'{fonts}',P.fonts_file);
x = regexprep(x,'{karyotype}',P.karyotype_file);
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

fprintf('Output name will be: %s\n',pngname);

save_textfile(x,conffile);

% run CIRCOS

fprintf('Calling CIRCOS\n');
cmd = 'perl';
if ~isempty(P.perl_modules_library) cmd = [cmd ' -I' P.perl_modules_library]; end
cmd = [cmd ' ' P.executable_name ' -conf "' conffile '"'];
fprintf('Calling command: %s\n',cmd);
[result output] = system(cmd);
fprintf('Result: '); disp(result);
disp(output);

