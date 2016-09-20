function generate_coverage_plot_and_bargraphs(M,outstem,P)

if ~exist('P','var'),P=[]; end
P = impose_default_value(P,'suppress_targets_with_membership_equals_zero',true);
P = impose_default_value(P,'output_per_exon_coverage',false);
P = impose_default_value(P,'output_per_exon_mutations',false);
P = impose_default_value(P,'superimpose_mutations',false);
P = impose_default_value(P,'output_pdfs',false);
P = impose_default_value(P,'impute_full_coverage',false);

if ~exist('outstem','var'), outstem=[]; end

if ~isempty(outstem) && outstem(end)~='/' && outstem(end)~='.', outstem = [outstem '.']; end

if iscell(M)
  error('please collapse M before submitting to this function');
end

close all

% DRAW COVERAGE PLOT

% make sure M has the bare minimum fields needed to allow capture_covplot to run
if ~isfield(M,'cov')
  M.cov = [];
end
if ~isfield(M.cov,'sample')
  M.cov.sample.name = M.patient.name;
  M.cov.ns = length(M.cov.sample.name);
end
if ~isfield(M.cov,'targ')
  M.cov.targ = [];
  M.cov.targ.gene = {'---'};
  M.cov.targ.gc = [1];
  M.cov.nt = 1;
end

PP=[];
if P.suppress_targets_with_membership_equals_zero
  if isfield(M.cov,'targ')
    if isfield(M.cov.targ,'col7'), M.cov.targ = rename_field(M.cov.targ,'col7','mem'); end
    if isfield(M.cov.targ,'mem')    % suppress targets not in the baitset
      PP.omit_these_targets = find(M.cov.targ.mem==0);   
    end
  end
end
PP=impose_default_value(PP,'sample_name_xadj',0);
PP=impose_default_value(PP,'sample_name_fontsize',8);
PP=impose_default_value(PP,'colorbar_fontsize',7);
PP=impose_default_value(PP,'GC_plot_fontsize',7);
if P.impute_full_coverage
  PP.omit_these_samples = [];
  PP.impute_full_coverage = true;
else
  PP.omit_these_samples = find(~ismember(M.cov.sample.name,M.patient.name));
  PP.impute_full_coverage = false;
end

if P.output_per_exon_coverage
  if ~isempty(outstem)
    PP.dump_cov_filename = [outstem 'per_exon.coverage.txt'];
  end
end
if P.output_per_exon_mutations
  PP.mutations = M.mut;
  PP.mutations.patient = M.mut.patient_name;
  PP.mutations_match_margin = 2;
  PP.superimpose_mutations = P.superimpose_mutations;
  if ~isempty(outstem)
    PP.dump_mut_filename = [outstem 'per_exon.mutation_counts.txt'];
  end
end

% multiM-mode: lastsort = iset
if isfield(M.cov,'sample') && isfield(M.cov.sample,'isetname')
  [PP.x_axis_lastsort_labels tmp PP.x_lastsort] = unique(M.cov.sample.isetname);
  tmp = distinct_colors(length(PP.x_axis_lastsort_labels)+1);
  PP.x_axis_lastsort_colors = tmp(2:end,:); % (no black)
end

figure(1);clf;
snameidx = capture_covplot(M.cov,PP);

if ~isempty(outstem)
  print('-dpng','-r120',[outstem 'covplot.png']);
  print('-dpng','-r360',[outstem 'covplot_hires.png']);
  if P.output_pdfs
    print_to_file([outstem 'covplot.pdf']);
  end
end

% DRAW BARGRAPHS

PP=[];
PP=impose_default_value(PP,'plot_ypos',[0.70 0.40 0.10]);
PP=impose_default_value(PP,'plot_width',0.91);
PP=impose_default_value(PP,'plot_left',0.07);
PP=impose_default_value(PP,'legend','best');
PP=impose_default_value(PP,'y_axis_fontsize',7);
PP=impose_default_value(PP,'ylabel_fontsize',12);
PP=impose_default_value(PP,'sample_name_fontsize',8);
PP=impose_default_value(PP,'off_scale_fontsize',5);
PP=impose_default_value(PP,'text_xadj',-0.1);
PP=impose_default_value(PP,'barwidth',0.7);

if P.impute_full_coverage
  D = M.patient;
else
  D = reorder_struct(M.cov.sample,snameidx);
end

% if "dataset" is available, use it to set background colors
if isfield(D,'dataset') && ~isfield(D,'label')
  [u ui uj] = unique(D.dataset);
  if length(u)>1
    D.label = uj;
    PP.label_names = u;
  end
end

idx = listmap(D.name,M.patient.name);
if any(isnan(idx)), fprintf('unexpected problem with sample name mapping'); keyboard; error('error'); end
D = merge_structs({reorder_struct(M.patient,idx),D}); % (mostly to get "type")
D = rmfield_if_exist(D,{'cidx','pidx','covfile','fwbfile'});
if isfield(M,'N_cov')
  D.N_tot = sum(squeeze(M.N_cov(:,M.TOT,idx)))';
else
  fprintf('mutation rate bargraphs: imputing full coverage\n');
  D.N_tot = repmat(sum(M.N_terr(:,M.TOT)),length(D.name),1);
end
D.nsil_tot = sum(squeeze(M.n_silent(:,M.TOT,idx)))';
D.nnon_tot = sum(squeeze(M.n_nonsilent(:,M.TOT,idx)))';

% if information is available, count how many of the mutations are at dbSNP sites
if isfield(M.mut,'dbSNP_RS')
  try
    dbsnp_sil = M.use_silent(grep('rs',M.mut.dbSNP_RS(M.use_silent),1));
    dbsnp_non = M.use_nonsilent(grep('rs',M.mut.dbSNP_RS(M.use_nonsilent),1));
    D.ndbsnp_tot = zeros(slength(D),1);
    if isfield(M.mut,'pat_idx')
      pidx = M.mut.pat_idx;
    elseif isnumeric(M.mut.patient)
      pidx = M.mut.patient;
    else
      error('couldn''t find pat_idx field');
    end
    for i=1:length(dbsnp_sil)
      p = find(snameidx==pidx(dbsnp_sil(i)));
      if ~isempty(p)
        D.nsil_tot(p) = D.nsil_tot(p) - 1;
        D.ndbsnp_tot(p) = D.ndbsnp_tot(p) + 1;
      end
    end
    for i=1:length(dbsnp_non)
      p = find(snameidx==pidx(dbsnp_non(i)));
      if ~isempty(p)
        D.nnon_tot(p) = D.nnon_tot(p) - 1;
        D.ndbsnp_tot(p) = D.ndbsnp_tot(p) + 1;
      end
    end
    D.rate_dbsnp = D.ndbsnp_tot ./ D.N_tot;
  catch me
    fprintf('Error processing dbSNP information\n');
    D = rmfield_if_exist(D,{'ndbsnp_tot','rate_dbsnp'});
  end
else
  fprintf('(No dbSNP information available for mutations\n');
end 

D.rate_sil = D.nsil_tot ./ D.N_tot;
D.rate_non = D.nnon_tot ./ D.N_tot;
if ~isempty(outstem)
  save_struct(D,[outstem 'patients.counts_and_rates.txt']);
end

figure(2);clf;
display_mutation_bargraph2(D,PP);

if ~isempty(outstem)
  print('-dpng','-r120',[outstem 'bargraphs.png']);
  print('-dpng','-r360',[outstem 'bargraphs_hires.png']);
  if P.output_pdfs
    print_to_file([outstem 'bargraphs.pdf']);
  end
  close all   % (so FH job terminates properly)
end


