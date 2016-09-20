function Ntot = calc_Ntot_from_wig_and_categ(W,Z,C)

maxsamps = double(max(cellfun(@max,W)));
Ntot = zeros(slength(Z),1);
for chr=1:24, fprintf('chr%d ',chr);
  minlen = min(length(C{chr}),length(W{chr}));
  Ntot = Ntot + hist2d_fast(double(C{chr}(1:minlen)),double(W{chr}(1:minlen)),1,slength(Z),0,maxsamps) * [0:maxsamps]';
end,fprintf('\n');
