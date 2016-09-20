function [N,G]=add_sampleinfo_genotypes(allelefile,sampleinfo,genotypefile,ABcn)
%ABcn=1 is flag to output original allele copy numbers with genotype variable, A is .adat(:,:1), B is .adat(:,:,2)
% addpath ~/CancerGenomeAnalysis/trunk/matlab/
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
% addpath ~/CancerGenomeAnalysis/trunk/matlab/POS
% addpath ~/CancerGenomeAnalysis/trunk/matlab/Startup


M1=read_modelled_data_file(allelefile,[],-1,1,0,0,0,0,10000);
M1.adat(:,:,1)=M1.dat(:,1:2:size(M1.dat,2));
M1.adat(:,:,2)=M1.dat(:,2:2:size(M1.dat,2));
M1.max=max(M1.adat,[],3);
M1.min=min(M1.adat,[],3);
M1.dat=M1.min+M1.max

if isempty(ABcn)
ABcn=0;
end

%read in sample info
SI=read_sample_info_file(sampleinfo,1);

%match sample info
N=M1;
N.sdesc=M1.sdesc(1:2:size(M1.sdesc,2));
M1.sdesc=M1.sdesc(1:2:size(M1.sdesc,2));
M1=rmfield(M1,'min');
M1=rmfield(M1,'max');

for k=1:size(N.adat,2);
    pat(k)=regexprep(N.sdesc(k),'_a(\w*)A', '');
end
idx=1:length(SI);
k=1:size(N.adat,2);
[Mt,m1,m2]=match_string_sets(pat(k),{SI(idx).array});
SIsub=SI(m2);
N.adat(:,:,1)=N.min;
N.adat(:,:,2)=N.max;
N2=reorder_D_cols(N, m1); 
N1=N2;
N1.sdesc=SIsub;
clear N2 
N1=rmfield(N1,'min');
N1=rmfield(N1,'max');

M1=reorder_D_cols(M1,m1);
M1.sdesc=SIsub;

%read in genotypes 
G=read_modelled_data_file(genotypefile,[],-1,0,0,0,0,0,10000)
[GN,idg,idn]=intersect(G.marker,N1.marker);
Gn=reorder_D_rows(G,idg);
Nn=reorder_D_rows(N1,idn);
M1=reorder_D_rows(M1,idn);
[aGN,n,g]=match_string_sets({N1.sdesc.array},G.sdesc);
Gn1=reorder_D_cols(Gn, g);

%need to match samples also, 
%if genotype file is only normals - tumors automatically are given value of 0 
Nn.calls(1:length(Nn.marker),1:size(Nn.adat,2))=0;
Nn.calls(:,n)=Gn1.dat;
Nn.cmarker=Gn1.marker;
clear G1 Gn G N1
M=Nn;
M.sis=M.sdesc;
M.sdesc={M.sis.array};

M1.sis=M1.sdesc;
M1.sdesc={M1.sis.array};
M1=order_by_pos(M1);

Mord=order_by_pos(M);
N.adat=Mord.adat;
N.marker=Mord.marker;
N.pos=Mord.pos;
N.chr=Mord.chr;
N.chrn=Mord.chrn;
N.sdesc=Mord.sdesc;
N.sis=Mord.sis;
N.dat=Mord.dat

G.marker=M.marker;
G.dat=M.calls;
G.sis=M.sis;
G.sdesc=M.sdesc;
G.pos=M.pos;
G.chr=M.chr;
G.chrn=M.chrn;
G=order_by_pos(G);
G=rmfield(G,'orig');
G=rmfield(G,'origidx');
G=rmfield(G,'gorigidx');
if ABcn==1
G.adat=M1.adat;
end