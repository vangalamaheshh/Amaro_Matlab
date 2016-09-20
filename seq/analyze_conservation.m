function [pval pval_cum] = analyze_conservation(P)

if ~exist('P','var'), P=[]; end

P=impose_default_value(P,'summed_cov_fwb','*required*');
P=impose_default_value(P,'cons_fwb','*required*');
P=impose_default_value(P,'cons_missing_value_token',[]);
P=impose_default_value(P,'genes_to_analyze','*required*');
P=impose_default_value(P,'maf_file','*required*');
P=impose_default_value(P,'exclude_silent_mutations',true);
P=impose_default_value(P,'target_list','*required*');
P=impose_default_value(P,'outdir','*required*');
P=impose_default_value(P,'generate_figures',true);

import('org.broadinstitute.cga.tools.seq.FixedWidthBinary');

% open coverage and conservation tracks
COV = FixedWidthBinary(P.summed_cov_fwb);
COV_nullval = COV.getNullVal;
CONS = FixedWidthBinary(P.cons_fwb);
CONS_nullval = CONS.getNullVal;

% genes to analyze
G=[];
G.name = P.genes_to_analyze;
ng = slength(G);

% load targets
Ta = load_target_file(P.target_list);
Ta.gidx = listmap(Ta.gene,G.name);
Ta = reorder_struct(Ta,~isnan(Ta.gidx));
Ta.chr = convert_chr(Ta.chr);
Ta = make_numeric(Ta,{'start','end'});

% load mutations
fprintf('Loading MAF file... ');
Ma = load_struct(P.maf_file);
Ma = add_simple_fieldnames(Ma);
Ma.chr = convert_chr(Ma.chr);
Ma = make_numeric(Ma,{'start','end'});
if isfield(Ma,'gene_name'), Ma.gene = Ma.gene_name; end
Ma.gidx = listmap(Ma.gene,G.name);
Ma = reorder_struct(Ma,~isnan(Ma.gidx));
if P.exclude_silent_mutations
  Ma = reorder_struct(Ma,grepv('synon|silent',Ma.type,1));
end
fprintf('\n');

% ANALYSIS

if P.generate_figures
  figure(1)
  set(gcf,'position',[440 302 1236 498],'visible','off');
end

G.pval_conservation = nan(ng,1);
G.pval_conservation_cumulative = nan(ng,1);
ensure_dir_exists(P.outdir);

ccons_tot = [];
mcons_tot = [];

for g=1:ng
  fprintf('Gene %s: ',G.name{g});

  T = reorder_struct(Ta,Ta.gidx==g);
  nt = slength(T);
  if nt==0
    fprintf('no targets in target list!\n');
    continue;
  end
  T = sort_struct(T,'start');
  T.start_pos = find_cDNA_pos(T,T);

  M = reorder_struct(Ma,Ma.gidx==g);
  nm = slength(M);
%  if nm==0
%    fprintf('no mutations!\n');
%    continue;
%  end
  M.cDNA_pos = find_cDNA_pos(M,T);

  plusstrand = true;
  if (isfield(M,'Transcript_Strand') && any(strcmp(M.Transcript_Strand,'-')))
    plusstrand = false;
  end

  cons = double(CONS.get(T.chr,T.start,T.end));
  cons(cons==CONS_nullval) = NaN;
  if ~isempty(P.cons_missing_value_token)
    cons(cons==P.cons_missing_value_token) = NaN;
  end
  cov = double(COV.get(T.chr,T.start,T.end));
  cov(cov==COV_nullval) = 0;

  len = length(cov);

  [ucons ui uj] = unique(cons);
  mincons = min(ucons);
  maxcons = max(ucons);
  ch = zeros(length(ucons),1);
  C_cons = cell(len,1);
  for i=1:len
    ch(uj(i)) = ch(uj(i)) + cov(i);
    C_cons{i} = repmat(cons(i),cov(i),1); 
  end
  C_cons = cat(1,C_cons{:});

  M.cons = nansub(cons,M.cDNA_pos);
  idx = listmap(M.cons,ucons);
  mh = histc(idx,1:length(ucons));

  [h G.pval_conservation(g)] = ttest2(M.cons,C_cons,[],'right');

  mcons_tot = [mcons_tot; M.cons];
  ccons_tot = [ccons_tot; C_cons];
  [h G.pval_conservation_cumulative(g)] = ttest2(mcons_tot,ccons_tot,[],'right');

  if P.generate_figures
    clf
    subplot(3,1,1);
      hold on
      if plusstrand, range = 1:len; mutpos = M.cDNA_pos;
      else range = len:-1:1; mutpos = len-M.cDNA_pos+1; end
%      [AX,H1,H2] = plotyy(range,cov,range,cons);
%      set(get(AX(1),'Ylabel'),'String','Coverage (# samples)'); 
%      set(get(AX(2),'Ylabel'),'String','Conservation');
%      set(AX,'xlim',[1 len],'tickdir','out');
      scatter(range,cons,10,[0 0.5 0]);
      yl = ylim;
      xlim([1 len]); set(gca,'tickdir','out');
      title(G.name{g},'fontsize',20);
      xlabel('position (bp along concatenated exons)');
      ylabel('conservation','color',[0 0.5 0]);
      text(len*1.02,yl(2),num2str(max(cov)),'color',[0.2 0.2 1],'clipping','off');
      if plusstrand
        for i=2:2:nt
          rectangle('position',[T.start_pos(i-1) yl(1) T.start_pos(i)-T.start_pos(i-1) yl(2)-yl(1)],...
                    'facecolor',ones(1,3)*0.85);
        end
      else 
        for i=2:2:nt
          rectangle('position',[len-T.start_pos(i)+1 yl(1) T.start_pos(i)-T.start_pos(i-1) yl(2)-yl(1)],...
                    'facecolor',ones(1,3)*0.85);
        end
      end
%      legend({'cons','cov','muts'},'location','eastoutside');
      scatter(range,cons,10,[0 0.5 0]);
      scatter(range,cov/max(cov)*yl(2),10,[0.2 0.2 1]);
      scatter(mutpos,M.cons,30,[1 0 0],'filled');
      hold off
    subplot(3,1,2);
      bar(ucons,ch,'facecolor',[0.2 0.2 1]);
      xlim([mincons-0.5 maxcons+0.5]);
      ylabel('coverage (# bp)','color',[0.2 0.2 1]);
      xlabel('conservation','color',[0 0.5 0]);
      set(gca,'tickdir','out');
    subplot(3,1,3);
      bar(ucons,mh,'facecolor',[1 0 0]);
      xlim([mincons-0.5 maxcons+0.5]);
      ylabel('mutations (#)','color',[1 0 0]);
      xlabel('conservation','color',[0 0.5 0]);
      set(gca,'tickdir','out');
      ylim([0 max(2,max(mh))]);
      yl = ylim;
      if G.pval_conservation(g) < 0.001, fmt = '%d'; else fmt = '%0.3f'; end
      text(maxcons+0.5,yl(2)+diff(yl)*0.15,sprintf(['p = ' fmt],G.pval_conservation(g)),'color',[1 0 0],...
           'horizontalalignment','right','verticalalignment','middle','fontsize',15,'clipping','off');
      outname = [P.outdir sprintf('/gene%05d',g) '_conservation_' G.name{g} '.png'];
    print_to_file(outname);
  else
    fprintf('\n');
  end


end % next gene

reportname = [P.outdir '/conservation_report.txt'];
save_struct(G,reportname);

if P.generate_figures
  close all
end

if nargout>=1
  pval = G.pval_conservation;
end
if nargout>=2
  pval_cum = G.pval_conservation_cumulative;
end




