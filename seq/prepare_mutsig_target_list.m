function prepare_mutsig_target_list(outfile,build,splicesiteflank,threshold,exome_interval_list)
% prepare_mutsig_target_list
%
% combines refseq with baitset target interval list
%

if ~exist('outfile','var'), error('outfile required'); end
if ~exist('build','var'), build = 'hg18'; end
if ~exist('threshold','var'), threshold = 0.20; end
if ~exist('splicesiteflank','var'), splicesiteflank = 2; end
if ~exist('exome_interval_list','var')
  exome_interval_list = ['/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/' ...
                      'whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list'];
end

disp(build);

G = get_refseq_introns_and_exons(build,splicesiteflank);
E = convert_to_exon_list(G);
E = sort_struct(E,{'chr','start'});
E.len = E.end-E.start+1;

dr = '/seq/references/HybSelOligos';
fprintf(['Using "' exome_interval_list '" as whole-exome for calculation of "membership" column.\n']);
T = load_target_interval_list(exome_interval_list);
T = sort_struct(T,{'chr','start'});
T.len = T.end-T.start+1;

% compare lists

nt = slength(T); ne = slength(E);
T.eidx = nan(nt,1); E.tidx = nan(ne,1);
ei = 1; ti = 1;
while(ti<=nt && ei<=ne)
  % see if they overlap
  if T.chr(ti)==E.chr(ei) && T.start(ti)<=E.end(ei) && T.end(ti)>=E.start(ei)
    % compute fraction overlap = len(E)+len(T) / 2*len(overlap)
    if T.start(ti)>E.start(ei), overlap = T.end(ti)-E.start(ei)+1;
    else overlap = E.end(ei)-T.start(ti)+1; end
    frac = (2*overlap) / (T.len(ti)+E.len(ei));
    if frac>=threshold % if sufficient overlap
      T.eidx(ti) = ei;  E.tidx(ei) = ti;  % mark these as a match
  end,end
  % choose which pointer to increment
  if T.chr(ti)<E.chr(ei) || (T.chr(ti)==E.chr(ei) && T.start(ti)<E.start(ei)), ti=ti+1; else ei=ei+1; end
end
fprintf('Out of %d targets in exome baitset, %d lack counterparts in Refseq exon list.\n',slength(T),sum(isnan(T.eidx)));
fprintf('Out of %d exons in Refseq, %d lack counterparts in exome baitset.\n',slength(E),sum(isnan(E.tidx)));

% compare to what's in C2K and C6K

fc2k = '/xchip/tcga_scratch/lawrence/capture/cancer_2000gene_shift170.targets.interval_list_NUMSONLY.txt';
fc6k = '/xchip/tcga_scratch/lawrence/capture/tcga_6k_genes.targets.interval_list_NUMSONLY.txt';
c2k = rename_fields(load_struct(fc2k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = rename_fields(load_struct(fc6k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = reorder_struct(c6k,c6k.chr>=1 & c6k.chr<=24);
E.mem = understand_nested_target_sets({E,c6k,c2k},0.90);

hist2d_fast(E.mem,double(~isnan(E.tidx)),1,3,0,1)
%          not in target set   in target set
% in we only        8209      113506
% in c6k             855       45872
% in c2k             343       26364

% (if reduce cutoff from 0.90 to 0.50, numbers change only slightly:)
%        8208      113442
%         856       45873
%         343       26427


% add GC content
fprintf('Annotating GC content\n');
E.gc = get_gc_content(E.chr,E.start,E.end,build);

% save final list

fprintf('Saving\n');
L = keep_fields(E,{'gene','chr','start','end','gc','len','mem'});
L.mem(isnan(E.tidx)) = 0;

% list is list of all Refseq exons.
%   membership field (mem):
%      0 = not in WE baitset
%      1 = in WE baitset; not in C6K baitset
%      2 = in WE baitset; in C6K baitset; not in C2K baitset
%      3 = in WE baitset; in C6K baitset; in C2K baitset

save_struct(L,outfile,'no_headers');



