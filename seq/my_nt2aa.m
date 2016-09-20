function aa=my_nt2aa(nt)

COD =  {    'TTT', 'TTC', 'TTA', 'TTG', 'TCT',...
            'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT',...
            'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC',...
            'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA',...
            'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG',...
            'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT',...
            'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC',...
            'GGA', 'GGG', 'TAA', 'TAG', 'TGA' };

AA = {      'F', 'F', 'L', 'L', 'S', 'S',...
            'S', 'S', 'Y', 'Y', 'C', 'C', 'W', 'L', 'L', 'L', 'L', 'P', 'P',...
            'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'I',...
            'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',...
            'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G',...
            'G', 'G', 'G', '*', '*', '*' };

aa = [];

nt=nt(1:3*floor(length(nt)/3));
nt=upper(nt);
for i = 1:+3:length(nt)
   pos = find(strcmp(COD,nt(i:i+2)));
   if ~isempty(pos)
      aa = [aa AA{pos}];
   else
      aa = [aa '?'];
   end
end

