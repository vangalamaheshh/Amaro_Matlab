



function SCNA_CCF_chi2(A,sample_name,ex_seg,ABS_seg)

purity=A.purity(ismember(A.sample,sample_name));
ploidy=A.ploidy(ismember(A.sample,sample_name));

CCF_dist=[0:0.01:1];

fudge_factor=20;
ex_seg.sigma_major=ex_seg.sigma_major*fudge_factor;
ex_seg.sigma_minor=ex_seg.sigma_minor*fudge_factor;
ex_seg.sigma_tau=ex_seg.sigma_tau*fudge_factor;

for i=1:slength(ex_seg)
    chi=repmat(inf,[5 5 101]);
    chimin=inf;
    min_preNA=NaN;
    min_preNB=NaN;
    chi_premin=NaN;
    minNA=NaN;
    minNB=NaN;
    for NA=0:4
        for NB=NA:4
            if ~isnan(ex_seg.mu_major(i))
                mua=((purity*(CCF_dist*(NA-1)+1))+1-purity)/(purity*(ploidy/2)+(1-purity));
                mub=((purity*(CCF_dist*(NB-1)+1))+1-purity)/(purity*(ploidy/2)+(1-purity));
                chi(NA+1,NB+1,:)=((ex_seg.mu_minor(i)-mua)/ex_seg.sigma_minor(i)).^2+((ex_seg.mu_major(i)-mub)/ex_seg.sigma_major(i)).^2;
                
                
            else
                if NA==1 || NB==1
                    tau=2*((purity*(CCF_dist*(NA+NB-2)+2)+2*(1-purity))/(purity*ploidy+2*(1-purity)));
                    chi(NA+1,NB+1,:)=((ex_seg.tau(i)-tau)/ex_seg.sigma_tau(i)).^2;
                    
                end
            end

            if min(chi(NA+1,NB+1,:),[],3)<chimin
                min_preNA=minNA;
                min_preNB=minNB;
                chi_premin=chimin;
                chimin=min(chi(NA+1,NB+1,:),[],3);
                minNA=NA;
                minNB=NB;
            end
            
        end
    end
    L=exp(-chi(minNA+1,minNB+1,:)/2);
    CL=cumsum(squeeze(L)./sum(L(:)));
    
    ex_seg.minNA(i,1)=minNA;
    ex_seg.minNB(i,1)=minNB;
    ex_seg.min_preNA(i,1)=min_preNA;
    ex_seg.min_preNB(i,1)=min_preNB;
    ex_seg.chiPreMin(i,1)=chi_premin;
    ex_seg.chiMin(i,1)=chimin;
    ex_seg.minNA(i,1)=minNA;
    ex_seg.minNB(i,1)=minNB;
    ex_seg.CCF_dist(i,:)=squeeze(chi(minNA+1,minNB+1,:));
    ex_seg.CCF(i,1)=CCF_dist(find(squeeze(chi(minNA+1,minNB+1,:))==min(squeeze(chi(minNA+1,minNB+1,:))),1,'first'));
    ex_seg.CCF_Low(i,1)=CCF_dist(max([find(CL<=.025,1,'Last') 1]));
    ex_seg.CCF_High(i,1)=CCF_dist(min([find(CL<=.975,1,'First') length(CL)]));
    ex_seg.IS_SCNA(i,1)=ex_seg.CCF_Low(i,1)>0;
    
end

end
% 
 function test
 
 A=load_table('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/FullSetABSOLUTE.Table.CLL8.txt');
ex_seg=load_table('~/Downloads/CLL-GCLL-0027-Tumor-SM-41JOW.tsv');
ABS_seg=load_table('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/ABSOLUTE_results/FixFor2D_Full_CLL8_dataset_3.7/reviewed/SEG_MAF/CLL-GCLL-0027-TR-NT-SM-41JOW-SM-41QC8.segtab.txt');



 
 function SCNA_CCF_chi2(A,ex_seg,ABS_seg)
 
ex_seg.abs_ccf_a1=repmat(NaN,slength(ex_seg),1);
ex_seg.abs_ccf_a2=repmat(NaN,slength(ex_seg),1);

figure()
hold on
ex_seg.xpose=xhg19(ex_seg.Chromosome,ex_seg.End_bp);
ex_seg.xposs=xhg19(ex_seg.Chromosome,ex_seg.Start_bp);
ABS_seg.xposs=xhg19(ABS_seg.Chromosome,ABS_seg.Start_bp);
[i m]=ismember(ABS_seg.xposs,ex_seg.xposs);
ex_seg.abs_ccf_a1(m(m>0))=ABS_seg.cancer_cell_frac_a1(i);
ex_seg.abs_ccf_a2(m(m>0))=ABS_seg.cancer_cell_frac_a2(i);

for i=1:slength(ex_seg)
plot([ex_seg.xposs(i) ex_seg.xpose(i)],repmat(ex_seg.CCF(i),2,1))
plot([ex_seg.xposs(i) ex_seg.xpose(i)],repmat(ex_seg.f(i),2,1),'r')
plot([ex_seg.xposs(i) ex_seg.xpose(i)],repmat(ex_seg.tau(i)/10,2,1),'k')
if ex_seg.IS_SCNA(i)
plot([ex_seg.xposs(i) ex_seg.xpose(i)],repmat((ex_seg.minNA(i)+ex_seg.minNB(i))/10,2,1),'c')
end
plot([ex_seg.xposs(i) ex_seg.xpose(i)],repmat(ex_seg.abs_ccf_a2(i),2,1),'m:')
plot([ex_seg.xposs(i) ex_seg.xpose(i)],repmat(ex_seg.abs_ccf_a1(i),2,1),'g:')

end
 end
% 
% end


%
%
% for i=1:slength(ex_seg)
%
%     if ex_seg.tau(i)>2.1 || ex_seg.tau(i)<1.9 || ex_seg.f(i)<.45
%         ex_seg.pc_mu_minor(i,1)=ex_seg.mu_minor(i)*(purity*ploidy/2+(1-purity)*1);
%         ex_seg.pc_mu_major(i,1)=ex_seg.mu_major(i)*(purity*ploidy/2+(1-purity)*1);
%
%         for NA=0:4
%             for NB=NA:4
%                 CiMinor(NA+1,NB+1)=(ex_seg.pc_mu_minor(i)-1)/(purity*(NA-1));
%                 CiMajor(NA+1,NB+1)=(ex_seg.pc_mu_major(i)-1)/(purity*(NB-1));
%             end
%         end
%
%
%     end
% end
