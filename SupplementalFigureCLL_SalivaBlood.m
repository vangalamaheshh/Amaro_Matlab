Saliva = load_struct('~/amaro/deTiN_Figure_Code_and_PDFs/Supplement/CLL_194_Saliva_CD19_comparison/CLL_194_Saliva_Normal.call_stats.TiN.txt');
Blood = load_struct('~/amaro/deTiN_Figure_Code_and_PDFs/Supplement/CLL_194_Saliva_CD19_comparison/CLL_194_Blood_Normal.call_stats.TiN.txt');

true_positive_mutations = Saliva.xstart(ismember(Saliva.judgement,'KEEP'));
blood_mutations = Blood.xstart(ismember(Blood.judgement,'KEEP'));
sum(ismember(true_positive_mutations,blood_mutations))