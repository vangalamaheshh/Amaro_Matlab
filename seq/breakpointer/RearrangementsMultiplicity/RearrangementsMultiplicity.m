function RearrangementsMultiplicity(sample,tumor_cr_fn,normal_cr_dir,P,dRangerFile,cn) 
% Yotam Drier, yotamd@gmail.com

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'quiet',1);
P=impose_default_value(P,'purity','*required*');
P=impose_default_value(P,'ploidy','*required*');
flds={'drnum','first_by_pair','second_by_pair','fusion_by_pair','first_by_span','second_by_span','fusion_by_span','contradict'};
tumorcr=load_struct(tumor_cr_fn);
d=dir(normal_cr_dir);
j=[d.isdir];
j([1 2])=false;
if (isfield(tumorcr,'NO_DATA') || isempty(j))
    rm.NO_DATA={};
else
    tumorcr=make_numeric(tumorcr,flds);
    d={d(j).name};
    l=length(d);
    normaltot=zeros(1,l);
    for i=1:l
        normalcr(i)=make_numeric(load_struct([normal_cr_dir d{i} '/support.onenormal.txt']),flds);
        normaltot(i)=dlmread([normal_cr_dir d{i} '/readstotnum.onenormal.txt']);
    end
    mytot=normaltot(strcmp(d,sample));
    normaltot=normaltot./mytot;
    dr=make_numeric(load_struct(dRangerFile),{'num','chr1','chr2','pos1','pos2'});    
    %t=[cn.Startbp cn.Endbp];    
    dr.cn1=nan(length(dr.num),2);
    dr.cn2=dr.cn1;
    for c=1:23
        cncl=cn.Chromosome==c;
        if any(cncl)
            %sm=repmat(cn.modal_A1(cncl)+cn.modal_A2(cncl),1,2); % the sum of both alleles
            %sm=reshape(sm',numel(sm),1);
            %loc=reshape(t(cncl,:)',numel(t(cncl,:)),1);
            drcl1=dr.chr1==c;
            drcl2=dr.chr2==c;
            %dr.cn1(drcl1)=interp1(loc,sm,dr.pos1(drcl1),'nearest');
            %dr.cn2(drcl2)=interp1(loc,sm,dr.pos2(drcl2),'nearest');
            sm=cn.modal_A1(cncl)+cn.modal_A2(cncl); % the sum of both alleles
            dr.cn1(drcl1,:)=get_values_from_segments(cn.Startbp(cncl),cn.Endbp(cncl),sm,dr.pos1(drcl1));
            dr.cn2(drcl2,:)=get_values_from_segments(cn.Startbp(cncl),cn.Endbp(cncl),sm,dr.pos2(drcl2));
        end
    end     
    [c,ic,id]=intersect(tumorcr.drnum,dr.num);    
    tumorcr.cn1=nan(length(tumorcr.drnum),2);
    tumorcr.cn2=tumorcr.cn1;
    tumorcr.cn1(ic,:)=dr.cn1(id,:);
    tumorcr.cn2(ic,:)=dr.cn2(id,:);
    rm.drnum=tumorcr.drnum;
    
    c1=repmat(tumorcr.first_by_pair+tumorcr.first_by_span,1,2)./(P.purity*tumorcr.cn1+(1-P.purity)*2);
    c2=repmat(tumorcr.second_by_pair+tumorcr.second_by_span,1,2)./(P.purity*tumorcr.cn2+(1-P.purity)*2);
    x=(cell2mat({normalcr.first_by_pair})+cell2mat({normalcr.first_by_span}))/2;
    x=x./repmat(normaltot,size(x,1),1);
    meanhalfnorm=repmat(mean(x,2),1,2);
    p1=2*normcdf(-abs(c1-meanhalfnorm),0,repmat(std(x,[],2),1,2));
    j=[p1(:,1)>p1(:,2) p1(:,1)<=p1(:,2)]';
    p1=p1';
    p1=p1(j);
    c1=c1';
    c1=c1(j);
    m1=meanhalfnorm';
    m1=m1(j);
%     [p1,j]=max(p1,[],2);
%     c1=c1(:,j);
    x=(cell2mat({normalcr.second_by_pair})+cell2mat({normalcr.second_by_span}))/2;
    x=x./repmat(normaltot,size(x,1),1);
    meanhalfnorm=repmat(mean(x,2),1,2);
    p2=2*normcdf(-abs(c2-meanhalfnorm),0,repmat(std(x,[],2),1,2));
    j=[p2(:,1)>p2(:,2) p2(:,1)<=p2(:,2)]';
    p2=p2';
    p2=p2(j);
    c2=c2';
    c2=c2(j);
    m2=meanhalfnorm';
    m2=m2(j);
%     [p2,j]=max(p2,[],2);
%     c2=c2(:,j);
    rm.p1=p1;
    rm.p2=p2;
    rm.coeff1=c1;
    rm.coeff2=c2;
    rm.first=((tumorcr.first_by_pair+tumorcr.first_by_span)./c1-2*(1-P.purity))/P.purity;
    rm.second=((tumorcr.second_by_pair+tumorcr.second_by_span)./c2-2*(1-P.purity))/P.purity;
    rm.first_by_normals=((tumorcr.first_by_pair+tumorcr.first_by_span)./m1-2*(1-P.purity))/P.purity;
    rm.second_by_normals=((tumorcr.second_by_pair+tumorcr.second_by_span)./m2-2*(1-P.purity))/P.purity;
    t=2/P.purity-2;    
    rm.fusion_tot=2*(P.ploidy+t)*(tumorcr.fusion_by_pair+tumorcr.fusion_by_span)./(normalcr.first_by_pair+normalcr.second_by_pair+normalcr.first_by_span+normalcr.second_by_span);    
    rm.first_tot=(P.ploidy+t)*((tumorcr.first_by_pair+tumorcr.first_by_span)./(normalcr.first_by_pair+normalcr.first_by_span))-t;
    rm.second_tot=(P.ploidy+t)*((tumorcr.second_by_pair+tumorcr.second_by_span)./(normalcr.second_by_pair+normalcr.second_by_span))-t;
    rm.somatic_strength_tot=(tumorcr.fusion_by_pair+tumorcr.fusion_by_span)./(normalcr.fusion_by_pair+normalcr.fusion_by_span);
    rm.first_tot(tumorcr.first_by_pair+tumorcr.first_by_span==0)=0;
    rm.second_tot(tumorcr.second_by_pair+tumorcr.second_by_span==0)=0;
    rm.somatic_strength_tot(tumorcr.fusion_by_pair+tumorcr.fusion_by_span==0)=0;
    rm.fusion_tot(tumorcr.fusion_by_pair+tumorcr.fusion_by_span==0)=0;
    rm.fusion_by_pair=2*(P.ploidy+t)*tumorcr.fusion_by_pair./(normalcr.first_by_pair+normalcr.second_by_pair);    
    rm.first_by_pair=(P.ploidy+t)*(tumorcr.first_by_pair./normalcr.first_by_pair)-t;
    rm.second_by_pair=(P.ploidy+t)*(tumorcr.second_by_pair./normalcr.second_by_pair)-t;
    rm.somatic_strength_by_pair=tumorcr.fusion_by_pair./normalcr.fusion_by_pair;
 %   rm.fusion_fraction_by_pair=2*tumorcr.fusion_by_pair./(tumorcr.first_by_pair+tumorcr.second_by_pair-2*t/(P.ploidy+t)); % tumorcr.first coefficient should be phi_fusion/phi_first, assuming =1 (for second as well)
    rm.first_by_pair(tumorcr.first_by_pair==0)=0;
    rm.second_by_pair(tumorcr.second_by_pair==0)=0;
    rm.somatic_strength_by_pair(tumorcr.fusion_by_pair==0)=0;
    rm.fusion_by_pair(tumorcr.fusion_by_pair==0)=0;
   % rm.fusion_fraction_by_pair(tumorcr.fusion_by_pair==0)=0;
    rm.fusion_by_span=2*(P.ploidy+t)*tumorcr.fusion_by_span./(normalcr.first_by_span+normalcr.second_by_span);    
    rm.first_by_span=(P.ploidy+t)*(tumorcr.first_by_span./normalcr.first_by_span)-t;
    rm.second_by_span=(P.ploidy+t)*(tumorcr.second_by_span./normalcr.second_by_span)-t;
    rm.somatic_strength_by_span=tumorcr.fusion_by_span./normalcr.fusion_by_span;
   % rm.fusion_fraction_by_span=2*tumorcr.fusion_by_span./(tumorcr.first_by_span+tumorcr.second_by_span-2*t/(P.ploidy+t));     
    rm.first_by_span(tumorcr.first_by_span==0)=0;
    rm.second_by_span(tumorcr.second_by_span==0)=0;
    rm.somatic_strength_by_span(tumorcr.fusion_by_span==0)=0;
    rm.fusion_by_span(tumorcr.fusion_by_span==0)=0;
   % rm.fusion_fraction_by_span(tumorcr.fusion_by_span==0)=0;
end
if exist('dRangerFile','var')
    if ~strcmp(dRangerFile,'none')
        dr=make_numeric(load_struct(dRangerFile),{'num'});
        [c,id,ir]=intersect(dr.num,rm.drnum);
        dr.RMfirst=nan(size(dr.num));
        dr.RMsecond=nan(size(dr.num));
        dr.RMfusion=nan(size(dr.num));
        dr.RMsomatic_strength=nan(size(dr.num));
        dr.RMfirst(id)=rm.first_tot(ir);
        dr.RMsecond(id)=rm.second_tot(ir);
        dr.RMfusion(id)=rm.fusion_tot(ir);
        dr.RMsomatic_strength(id)=rm.somatic_strength_tot(ir);
        %if exist('AbsoluteRes','var')
            dr=make_numeric(dr,{'chr1','chr2','pos1','pos2'});
            %cn=load_struct(cnFile);
            %j=strcmp(birdseedname,AbsoluteRes.Sample);
            %cn=structfun(@(x)x(j),AbsoluteRes,'UniformOutput',false);
            %cn=make_numeric(cn,{'chrom','locstart','locend','segmean'});
            %cn=make_numeric(cn,{'Chromosome','PhysicalPosition','Signal'});
            %cn=make_numeric(cn,{'Chromosome','Startbp','Endbp','modal_A1','modal_A2'});
            dr.cn1=zeros(size(dr.pos1));
            dr.cn2=dr.cn1;
            %t=[cn.locstart cn.locend];
            t=[cn.Startbp cn.Endbp];    
            %t=[AbsoluteRes.Startbp AbsoluteRes.Endbp];    
            %for c=1:24
            for c=1:23
                %c
                %cncl=cn.chrom==c;
                cncl=cn.Chromosome==c;
                %sm=repmat(AbsoluteRes.segmean(cncl),1,2);
                sm=repmat(cn.modal_A1(cncl)+cn.modal_A2(cncl),1,2);
                sm=reshape(sm',numel(sm),1);
                loc=reshape(t(cncl,:)',numel(t(cncl,:)),1);
                %sm=AbsoluteRes.Signal(cncl);
                %[loc,i]=unique(AbsoluteRes.PhysicalPosition(cncl));
                %sm=sm(i);
                drcl1=dr.chr1==c;
                drcl2=dr.chr2==c;                
                dr.cn1(drcl1)=interp1(loc,sm,dr.pos1(drcl1),'nearest');
                dr.cn2(drcl2)=interp1(loc,sm,dr.pos2(drcl2),'nearest');
            end
        %end
        save_struct(dr,'dRanger_results.detail.somatic.withRM.txt')
    end
end
save_struct(rm,'RearrangementsMultiplicity.support.txt')
