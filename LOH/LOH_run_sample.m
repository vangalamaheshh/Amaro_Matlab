function LOH_run_sample(normal_genotype_path, tumor_allele_path, outdir)
%addpath
%/xchip/cle/analysis/callahan/CancerGenomeAnalysis/sandbox/callahan/LOH/
addpath /xchip/cle/analysis/callahan/CancerGenomeAnalysis/trunk/matlab/LOH

disp('Hello world')

scratchdir = outdir;
cd ~/CancerGenomeAnalysis/trunk/matlab/
startup

% Load byallele file

%Hybridization REF                       SOURS_p_TCGAb22_SNP_N_GenomeWideSNP_6_H11_529942        SOURS_p_TCGAb22_SNP_N_GenomeWideSNP_6_H11_529942
%CompositeElement REF    Chromosome      PhysicalPosition        Signal_A        Signal_B
%SNP_A-8575125   1       554484  1.955   0.031
%SNP_A-8575115   1       554636  0.01    1.847

infid = fopen(tumor_allele_path,'r');
Sample = fscanf(infid, '%*s%*s%s', 1);
C = textscan(infid, '%s%d%d%f%f','Delimiter','\t','EndOfLine','\n','HeaderLines',2,'TreatAsEmpty',{'NA'});
[marker,chrn,pos,adat1,adat2]=deal(C{:});
fclose(infid);

N.marker = marker;
N.chr = strtrim(cellstr(num2str(chrn)));
N.chrn = chrn;
N.pos = pos;
N.adat(:,1,1) = min([adat1,adat2],[],2);
N.adat(:,1,2) = max([adat1,adat2],[],2);
N.sdesc = {Sample};

% Load Birdseed file

% 
% Hybridization REF       THING_p_TCGA_B19_SNP_N_GenomeWideSNP_6_A01_495116       THING_p_TCGA_B19_SNP_N_GenomeWideSNP_6_A01_495116
% CompositeElement REF    Call    Confidence
% SNP_A-2131660   0       0.0036
% SNP_A-1967418   2       0.0252


infid = fopen(normal_genotype_path,'r');
C = textscan(infid, '%s%d%f','Delimiter','\t','EndOfLine','\n','HeaderLines',2,'TreatAsEmpty',{'NA'});
[marker,call,conf]= deal(C{:});
fclose(infid);

G.marker = marker;
G.dat = call;
G.sdesc = {Sample};


% Align markers
[GN,idg,idn]=intersect(G.marker,N.marker);
G=reorder_D_rows(G,idg);
N=reorder_D_rows(N,idn);

G.pos = N.pos;
G.chr = N.chr;
G.chrn = N.chrn;

G = order_by_pos(G);
N = order_by_pos(N);


% NaN out the non hets
col = 1; % this assumes a single sample...
[non_het_idx, cols] = find (G.dat(:,col) ~= 1);
N.adat(non_het_idx,col,:) = NaN;

% Bundle up for CBS

Ntot=N;
Nmin=N;
Nmax=N;
Ntot.dat=N.adat(:,:,1)+Ntot.adat(:,:,2);
Nmin.dat = N.adat(:,:,1);
Nmax.dat = N.adat(:,:,2);
N.dat = Ntot.dat;
Ntot = rmfield(Ntot, 'adat');
Nmin = rmfield(Nmin, 'adat');
Nmax = rmfield(Nmax, 'adat');


%Run CBS
disp('running CBS...')
RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100*clock)));
sample_number=num2str(int32((rand(1, 1)+1)*100000));
fnTot = [scratchdir '/' Sample '_tot_' sample_number];
fnMin = [scratchdir '/' Sample '_min_' sample_number];
fnMax = [scratchdir '/' Sample '_max_' sample_number];

keyboard 

Nmin_path = run_cbs(fnMin,Nmin,[],1,0,100,1, 1,'local');
Nmax_path = run_cbs(fnMax,Nmax,[],1,0,100,1,1,'local');
Ntot_path = run_cbs(fnTot, Ntot, [], 1, 0, 100, 1, 1, 'local');

%read_segallele_data modified and renamed to smooth_segallele_data
Sc = smooth_segallele_data2(N, Nmin_path, Nmax_path, Ntot_path, 10, 'combine')

Sample =[Sample '_' sample_number];
[CN, L] = LOH_byallele2(Sc);
write_allele_data(Sc, [scratchdir '/' Sample], CN, [scratchdir '/' Sample '.cn']);
write_allele_data(Sc, [scratchdir '/' Sample], L, [scratchdir '/' Sample '.all']);




