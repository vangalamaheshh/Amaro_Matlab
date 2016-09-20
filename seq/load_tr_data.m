function x = load_tr_data(indiv,windowsize)

if ~exist('windowsize','var'), windowsize=100;end

dr = '/xchip/cga1/firehose_output/Individual';
files = regexprep({[dr '/xx/xx-Tumor/wgs/tr/xx-Tumor.for.wig'],[dr '/xx/xx-Tumor/wgs/tr/xx-Tumor.rev.wig'],...
   [dr '/xx/xx-Normal/wgs/tr/xx-Normal.for.wig'],[dr '/xx/xx-Normal/wgs/tr/xx-Normal.rev.wig']},'xx',indiv);
demand_file(files)

g=org.broadinstitute.cga.tools.seq.Genome();
g.loadWiggle(files{1}); tf = extract_wins(g,windowsize);
g.clear; g.loadWiggle(files{2}); tr = extract_wins(g,windowsize);
g.clear; g.loadWiggle(files{3}); nf = extract_wins(g,windowsize);
g.clear; g.loadWiggle(files{4}); nr = extract_wins(g,windowsize);

x = [tf';tr';nf';nr'];
