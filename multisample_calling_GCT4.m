%multisample calling for GCT4

% data loading and processing
C1=load_struct('~/Downloads/GCT_4_MultiSample_Call/TGC-Testes_DFCI_4-TP-NT-SM-4PDD8-SM-4PDDA.call_stats.txt');
C2=load_struct('~/Downloads/GCT_4_MultiSample_Call/TGC-Testes_DFCI_4-TP-NT-SM-4PDD9-SM-4PDDA.call_stats.txt');

C1.tumor_f=str2double(C1.tumor_f);
C1.t_alt_count=str2double(C1.t_alt_count);
C1.t_ref_count=str2double(C1.t_ref_count);
C2.t_alt_count=str2double(C2.t_alt_count);
C2.t_ref_count=str2double(C2.t_ref_count);
C2.tumor_f=str2double(C2.tumor_f);
C1.n_alt_count=str2double(C1.n_alt_count);
C1.n_ref_count=str2double(C1.n_ref_count);
C2.n_alt_count=str2double(C2.n_alt_count);
C2.n_ref_count=str2double(C2.n_ref_count);
C1.normal_f=str2double(C1.normal_f);
C2.normal_f=str2double(C2.normal_f);


seg1=load_struct('~/Downloads/TCG-Testes_DFCI_4-Tumor-SM-4PDD8.tsv');
seg2=load_struct('~/Downloads/TCG-Testes_DFCI_4-Tumor-SM-4PDD9.tsv');

% values derived from ABSOLUTE and deTiN
Purity=[.28,.97];
TiN_uo=[.5,.2];
af=0:.01:1;

% priors based of TiN and Purity calculations
prior_germline=betapdf(af,25,25);      % prior on germline events 
prior_germline_lost=betapdf(af,40,60); % prior on LOHd alleles
prior_germline_gain=betapdf(af,60,40); 

prior_C1_somatic_clonal=betapdf(af,14+1,86+1); % prior on somatic events in tumor1
prior_C2_somatic_clonal=betapdf(af,50,50); % prior on somatic events in tumor2
prior_Normal_Somatic_clonal=betapdf(af,7,93); % prior on somatic events in normal

% represent events as beta pdfs
for i=1:slength(C1)
    C1.pdf(i,:)=betapdf(af,C1.t_alt_count(i)+1,C1.t_ref_count(i)+1);
    C1.normpdf(i,:)=betapdf(af,C1.n_alt_count(i)+1,C1.n_ref_count(i)+1);
end
for i=1:slength(C2)
    C2.pdf(i,:)=betapdf(af,C2.t_alt_count(i)+1,C2.t_ref_count(i)+1);
    C2.normpdf(i,:)=betapdf(af,C2.n_alt_count(i)+1,C2.n_ref_count(i)+1);
end


% likelihood for classification
for i=1:slength(C1)
    C1.p_C1_somatic(i,1)=dot(prior_C1_somatic_clonal,C1.pdf(i,:));
    C1.p_germ(i,1)=max([dot(prior_germline,C1.pdf(i,:)),dot(prior_germline_gain,C1.pdf(i,:)),dot(prior_germline_lost,C1.pdf(i,:))]);
    C1.p_n_somatic(i,1)=dot(prior_Normal_Somatic_clonal,C1.normpdf(i,:));
    C1.p_n_normal(i,1)=max([dot(prior_germline,C1.normpdf(i,:)),dot(prior_germline_gain,C1.normpdf(i,:)),dot(prior_germline_lost,C1.normpdf(i,:))]);
end


for i=1:slength(C2)
    C2.p_C2_somatic(i,1)=dot(prior_C2_somatic_clonal,C2.pdf(i,:));
    C2.p_germ(i,1)=max([dot(prior_germline,C2.pdf(i,:)),dot(prior_germline_gain,C2.pdf(i,:)),dot(prior_germline_lost,C2.pdf(i,:))]);
    C2.p_n_somatic(i,1)=dot(prior_Normal_Somatic_clonal,C2.normpdf(i,:));
    C2.p_n_normal(i,1)=max([dot(prior_germline,C2.normpdf(i,:)),dot(prior_germline_gain,C2.normpdf(i,:)),dot(prior_germline_lost,C2.normpdf(i,:))]);
end

% normalzing 
C1.p_C1_somatic=C1.p_C1_somatic./(C1.p_C1_somatic+C1.p_germ);
C1.p_germ=C1.p_germ./(C1.p_C1_somatic+C1.p_germ);
C1.p_n_somatic=C1.p_n_somatic./(C1.p_n_somatic+C1.p_n_normal);
C1.p_n_normal=C1.p_n_somatic./(C1.p_n_somatic+C1.p_n_normal);
C2.p_C2_somatic=C2.p_C2_somatic./(C2.p_C2_somatic+C2.p_germ);
C2.p_germ=C2.p_germ./(C2.p_C2_somatic+C2.p_germ);
C2.p_n_somatic=C2.p_n_somatic./(C2.p_n_somatic+C2.p_n_normal);
C2.p_n_normal=C2.p_n_somatic./(C2.p_n_somatic+C2.p_n_normal);


somaticC1=(C1.p_C1_somatic>C1.p_germ);
somaticN1=(C1.p_n_somatic>C1.p_n_normal);
C1.key=strcat(C1.position,C1.contig);
C2.key=strcat(C2.position,C2.contig);
positions=C1.key(somaticC1&somaticN1);
positions=positions(ismember(positions,C2.key));
score_s=zeros(length(positions),1);
score_g=zeros(length(positions),1);

for i=1:length(positions)
    if ~isempty(ismember(positions{i},C2.key))
        score_s(i,1)=C2.p_n_somatic(ismember(C2.key,positions{i}))+C2.p_C2_somatic(ismember(C2.key,positions{i}));
        score_g(i,1)=C2.p_n_normal(ismember(C2.key,positions{i}))+C2.p_germ(ismember(C2.key,positions{i}));
    end
end
positions=positions(score_s>score_g);

% some plots
subplot(1,2,1)
hold on
plot(C1.tumor_f,C1.normal_f,'o','Color',[.5 .5 .5])
plot(C1.tumor_f(somaticC1),C1.normal_f(somaticC1),'g.','MarkerSize',10)
plot(C1.tumor_f(somaticC1&somaticN1),C1.normal_f(somaticC1&somaticN1),'b.','MarkerSize',10)
plot(C1.tumor_f(ismember(C1.key,positions)),C1.normal_f(ismember(C1.key,positions)),'r.','MarkerSize',10)
subplot(1,2,2)
hold on
plot(C2.tumor_f,C2.normal_f,'o','Color',[.5 .5 .5])
plot(C2.tumor_f(ismember(C2.key,positions)),C2.normal_f(ismember(C2.key,positions)),'r.','MarkerSize',10)

figure()
hist(C2.p_C2_somatic,1000)




