

af=[0:.01:1];
 p53T=load_table('~/Documents/TumorCallStats.p53.cervical.txt');
 p53B=load_table('~/Documents/BenignCallStats.p53.cervical.txt');
 Bsamples=unique(p53B.tumor_name);

 
 % for the analysis we ended up using the Fisher's test not the likelihood
 % model
 
for Bs=1:length(Bsamples)
    % each of the normals; all the normals together; the PoN
    for i=1:slength(p53T)
        
        
        pMutation=binopdf(p53T.t_alt_count(i),p53T.total_reads(i),af);
        k=(p53B.position==p53T.position(i)&ismember(p53B.tumor_name,Bsamples{Bs}));
        pArtifact=betapdf(af,p53B.t_alt_count(k)+1,p53B.t_ref_count(k)+1);
        pArtifact=pArtifact./sum(pArtifact);
        p53T.pMut(i,Bs)=dot(pMutation,pArtifact);
        X=[p53T.t_alt_count(i),p53T.t_ref_count(i);p53B.t_alt_count(k),p53B.t_ref_count(k)];
        p53T.pMutFish(i,Bs)=FisherExtestx(X,'gt');
       
        p53T.Xmatrix{i,Bs}=X;
        p53T.expectedA(i,Bs)=((sum(X(1,:))/sum(sum(X)))*(sum(X(:,1))/sum(sum(X)))*sum(sum(X)));
        
    end
    
end
 for i=1:slength(P53TK)
      pMutation=binopdf(P53TK.t_alt_count(i),P53TK.total_reads(i),af);
       pArtifact=betapdf(af,sum(p53B.t_alt_count(k))+1,sum(p53B.t_ref_count(k))+1);
       pArtifact=pArtifact./sum(pArtifact);
        k=(p53B.position==P53TK.position(i)&ismember(p53B.tumor_name,Bsamples{Bs}));
        P53TK.pMut(i,8)=dot(pMutation,pArtifact);
        X=[P53TK.t_alt_count(i),P53TK.t_ref_count(i);sum(p53B.t_alt_count(k)),sum(p53B.t_ref_count(k))];
        P53TK.pMutFish(i,8)=FisherExtestx(X,'gt');
        P53TK.Xmatrix{i,8}=X;
        P53TK.expectedA(i,8)=((sum(X(1,:))/sum(sum(X)))*(sum(X(:,1))/sum(sum(X)))*sum(sum(X)));
 end

for i=1:slength(P53TK)
[v indx]=max(P53TK.pMutFish(i,:));
P53TK.max_af(i,1)=P53TK.Xmatrix{i,indx}(2,1)/(P53TK.Xmatrix{i,indx}(2,1)+P53TK.Xmatrix{i,indx}(2,2));
[phat pci]=binofit(P53TK.Xmatrix{i,indx}(2,1),(P53TK.Xmatrix{i,indx}(2,1)+P53TK.Xmatrix{i,indx}(2,2)),.025);
P53TK.low_max_af(i,1)=P53TK.max_af(i)-pci(1);
P53TK.high_max_af(i,1)=pci(2)-P53TK.max_af(i);
end
 % Plotting for power point:
 
 P53TK.x=xhg19(P53TK.contig,P53TK.position);
 
p53B=reorder_struct(p53B,ismember(p53B.position,P53TK.position));
p53B.x=xhg19(p53B.contig,p53B.position);
errorbar(P53TK.x,P53TK.max_af,P53TK.low_max_af,P53TK.high_max_af,'b.')
set(gca,'yscale','log')
hold on
plot(P53TK.x,P53TK.tumor_f,'r.')

xlabel('Genomic Position','FontSize',25)
ylabel('Allele Fraction','FontSize',25)

line([min(P53TK.x) max(P53TK.x)],[1e-6 1e-6])
legend({'Benign AF','Tumor AF','Accuracy limit Safe-SeqS'})

