function draw_mutation_spectrum_3d_barplot(Nn)
% draw_mutation_spectrum_3d_barplot(Nn)
%
% input is a 32x5 table
% where 32 rows = trinucleotide contexts around A or C
%                 as in the first 32 rows of /xchip/cga1/lawrence/db/context65/categs.txt
% columns:  column1 = coverage per context
%           column2 = number of mutations to A
%           column3 = # to C
%           column4 = # to G
%           column5 = # to T


% convert to table for display as 3d bargraph
% 2x3  rows: @C, @A;   cols = transition; flip transversion; skew transversion
%      4x4:   rows: 5'-base   cols: 3'-base

categs = load_struct('/xchip/cga1/lawrence/db/hg19/context65/categs.txt');
class = {'C->T', 'C->A', 'C->G', 'A->G', 'A->T', 'A->C'};
colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0.5 0.3 0.7;0 0.2 0.8];
bases = 'ACGT';base_order='TCAG';fridx = [find(base_order=='A') find(base_order=='C')];
transition='GTAC'; flip='TGCA'; skew='CATG'; change = {flip,skew,transition};
ttspace = 5; tcols = 3*4; trows = 2*4;
bign = nan(trows,tcols); bigN = bign; bigC = ones(trows,tcols,3); sx=1; sy=1;
ni = zeros(2*4,3*4);
for frombase=1:2, for base5=1:4, for base3=1:4
      row = find(strcmp(categs.name,[bases(frombase) ' in ' base_order(base5) '_' base_order(base3)]));
      for changetype=1:3
        tobase = change{changetype}(frombase); col = find(bases==tobase);
        xx = (sx-1)+(changetype-1)*4+base3; yy = (sy-1)+(2-frombase)*4+base5;
        bign(yy,xx) = Nn(row,col+1); bigN(yy,xx) = Nn(row,1);
        bigC(yy,xx,1:3) = colors(find(strcmp(class,[bases(frombase) '->' tobase])),:);
end,end,end,end
clf;bar3_with_colors(bign./bigN,bigC);
set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);


