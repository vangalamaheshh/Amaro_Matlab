


p_values = cell(slength(M{1}.gene), 4);
report = fopen('testing/lung_testing_4_scratch', 'wt');
import('org.broadinstitute.cga.tools.seq.FixedWidthBinary');
forbidden_1 = [-1,220,230];
forbidden_2 = [-1,220,240];
nsets = length(M);

setnames = cell(nsets,1);
for si=1:nsets
      setnames{si} = M{si}.name;
end

ghist = zeros(M{1}.ng,1);
for si=1:nsets
    ghist = ghist + histc(M{si}.mut.gene(M{si}.use_nonsilent),1:M{1}.ng);
end
[tmpzz, gene_order] = sort(ghist,'descend');
count_1 = 0;
count_2 = 0;
for i = gene_order' %grep('^TP53$',M{1}.gene.name,1) %1:M{1}.ng
    
%    try


    % GENE characteristics    
    name = M{1}.gene.name{i};
    chr = M{1}.gene.chr(i);
    starts = M{1}.cov.targ.start(M{1}.cov.targ.gidx == i);
    ends = M{1}.cov.targ.end(M{1}.cov.targ.gidx == i);
    exons_matrix = [starts'; ends']';
    genelength = M{1}.gene.len(i);
    genomic_array = get_GenomicArray(exons_matrix, genelength);


       
    
    % POLYPHEN TRACKS
    fileA = FixedWidthBinary('/xchip/cga1/reference/reference-maf/hg19/flow_runs/step_3/A.total.all.fwb');
    fileC = FixedWidthBinary('/xchip/cga1/reference/reference-maf/hg19/flow_runs/step_3/C.total.all.fwb');
    fileG = FixedWidthBinary('/xchip/cga1/reference/reference-maf/hg19/flow_runs/step_3/G.total.all.fwb');
    fileT = FixedWidthBinary('/xchip/cga1/reference/reference-maf/hg19/flow_runs/step_3/T.total.all.fwb');
    
    
    polyphen_track = populate_polyphen(genomic_array, fileA, fileC, fileG, fileT, chr);
    tmp = polyphen_track(:);
    tmp(tmp==220) = [];
    tmp(tmp==250) = [];
    msg = sprintf('Gene: %s',name);
    ok = true;
    if mean(tmp==230)>=0.8
        msg = [msg sprintf('\tscript uniprot mapping problem')];
        count_1 = count_1+1;
        ok = false;
    end
    if mean(tmp==240)>=0.5
        msg = [msg sprintf('\tmissing polyphen data')];
        count_2 = count_2+1;
        ok = false;
    end
    if ok && mean(tmp==230 | tmp==240)>=0.8
        msg = [msg sprintf('\tmixed difficulties')];
        count_1 = count_1+1;
        ok = false;
    end
    if ok
        msg = [msg char(9) 'OK'];
    end
    fprintf('%s\n',msg);    

end
fprintf('script uniprot mapping problem cases: %d ; polyphen uknown cases: %d', count_1, count_2);


