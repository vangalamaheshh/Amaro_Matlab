function generate_gene_coverage_plots(M,P)

require_fields(M,{'np','mut','cov'});

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'number_of_top_genes_to_generate_coverage_plots_for',20);
P=impose_default_value(P,'siggene_report_filename','*required*');
P=impose_default_value(P,'per_gene_coverage_plots_outdir','*required*');
P=impose_default_value(P,'power_per_gene_outstem',[]);

fprintf('Generating gene coverage plots...\n');
ngenes = P.number_of_top_genes_to_generate_coverage_plots_for;

% load significance calculations

S = load_struct(P.siggene_report_filename);
require_fields(S,{'gene','p','q','n','N'});

% load power calculations

pow = [3 5 10 15 20];
if ~isempty(P.power_per_gene_outstem)
  pow_outfile = [P.power_per_gene_outstem '.power_stats.txt'];
  powfull_outfile = [P.power_per_gene_outstem '.power_stats_full.txt'];
  x1 = load_struct(pow_outfile);
  x2 = load_struct(powfull_outfile);
  x1 = make_numeric(x1,{'x3','x5','x10','x15','x20'});
  x2 = make_numeric(x2,{'x3','x5','x10','x15','x20'});
  p1 = [x1.x3 x1.x5 x1.x10 x1.x15 x1.x20];
  p2 = [x2.x3 x2.x5 x2.x10 x2.x15 x2.x20];
end

% generate plots

outdir = P.per_gene_coverage_plots_outdir;
ensure_dir_exists(outdir);

close all;figure(3); clf
PP=[];
PP.mainxpos = 0.15;
PP.omit_sample_names = (M.np > 30);
PP.mutations = M.mut;
if isnumeric(PP.mutations.gene), PP.mutations.gene = PP.mutations.gene_name; end
if isnumeric(PP.mutations.patient), PP.mutations.patient = PP.mutations.patient_name; end
PP.mutation_sample_names_to_match = M.cov.sample.name;

for i=1:ngenes
  gname = S.gene{i};
  PP.genelist = {gname};
  PP.mainwidth = 0.68;
  axes('position',[0 0 1 1],'visible','off');
  capture_covplot(M.cov,PP);
  axes('position',[0 0 1 1],'visible','off');
  % print power statistics
  if ~isempty(P.power_per_gene_outstem)
    idx = find(strcmp(gname,x1.Gene));
  end
  for p=-9:length(pow)+3
    txt = ''; fmt = {};
    if p==-9, txt = gname; fmt = {'fontweight','bold'};
    elseif p==-8, txt = sprintf('rank = %d',i);
    elseif p==-6, txt = sprintf('n = %s',S.n{i});
    elseif p==-5, txt = sprintf('N = %s',S.N{i});
    elseif p==-4, txt = sprintf('p = %s',S.p{i});
    elseif p==-3, txt = sprintf('q = %s',S.q{i});
    elseif p==-1, txt = 'Detection Power'; fmt = {'fontweight','bold'};
    elseif p==0, txt = 'freq      power_A    power_F';
    else
      if ~isempty(P.power_per_gene_outstem)
        if p>=1 && p<=length(pow), txt = sprintf('%2d%%       %0.2f        %0.2f%t',pow(p),p1(idx,p),p2(idx,p));
        elseif p==length(pow)+2, txt = 'power_A = actual coverage';
        elseif p==length(pow)+3, txt = 'power_F = full coverage';
        end
      else
        txt='';
      end
    end
    text(0.8,0.5-0.05*p,txt,fmt{:});
  end
  print_to_file([outdir '/' num2str(i) '.png']);
end

fprintf('generate_gene_coverage_plots finished\n');

close all
