function n = write_hp_permstats(file,R,matcher,AC,G,cyto,treated_text,control_text)
if isfield(R,'pvalueL')
    filter = ~isnan(R.pvalueL);
else
    filter = ~isnan(R.pvalue);
end    
n = sum(filter);
[aci,cli] = find(matcher);
aci = aci(filter);
cli = cli(filter);
if isfield(R,'pvalueL')
    n = write_filtered_tabcols(file,[],...
                   {'hairpin', AC.hairpinID(aci)},...
                   {'gene', AC.geneID(aci)},...
                   {'band', {cyto(G.cyton(cli)).name}},...
                   {'chr', G.chrn(cli),'%d'},...
                   {'start base', G.gene_start(cli),'%d'},...
                   {'end base', G.gene_end(cli),'%d'},...
                   {[treated_text ' N'], R.ntreat(filter)},...
                   {[control_text ' N'], R.ncontrol(filter)},...
                   {[treated_text ' mean'], R.mean_t(filter)},...
                   {[control_text ' mean'], R.mean_c(filter)},...
                   {'left p-value',R.pvalueL(filter)},...
                   {'right p-value',R.pvalueR(filter)},...
                   {'right q-value',R.qvalueR(filter)},...
                   {'left q-value',R.qvalueL(filter)} );
else
    n = write_filtered_tabcols(file,[],...
                   {'hairpin', AC.hairpinID(aci)},...
                   {'gene', AC.geneID(aci)},...
                   {'band', {cyto(G.cyton(cli)).name}},...
                   {'chr', G.chrn(cli),'%d'},...
                   {'start base', G.gene_start(cli),'%d'},...
                   {'end base', G.gene_end(cli),'%d'},...
                   {[treated_text ' N'], R.ntreat(filter)},...
                   {[control_text ' N'], R.ncontrol(filter)},...
                   {[treated_text ' mean'], R.mean_t(filter)},...
                   {[control_text ' mean'], R.mean_c(filter)},...
                   {'p-value',R.pvalue(filter)},...
                   {'q-value',R.qvalue(filter)} );
end