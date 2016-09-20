function [P_all,S_all]=differential_analysis(D,cls0,cls1,test_type,no_p,maxsnps)


chunkdims = getmemchunkdims(D,'dat',2);
maxsnps = chunkdims(1);



P_all = zeros(getsize(D,'dat',1),1);
S_all = zeros(getsize(D,'dat',1),1);

k = 0;

while k < getsize(D,'dat',1)
    fprintf(1,'.');
    thisloopendidx = min(k + maxsnps,getsize(D,'dat',1));
    
    dat = D.dat((k+1):thisloopendidx,:);

    n=size(dat,1);

    if ischar(test_type)
        tmp.method=test_type;
        test_type=tmp;
    end


    switch test_type.method
        case 'ttest'
            [P,S]=ttest2_many_nan(dat,cls0,cls1,0,[],no_p); % S positive if
            % cls0>cls1
        case 'ttest_minvar'
            [P,S]=ttest2_many_nan(dat,cls0,cls1,0,test_type.minvar,no_p); % S positive if
            % cls0>cls1
        case 'ttest1side_minvar'
            [P,S]=ttest2_many_nan(dat,cls0,cls1,1,test_type.minvar,no_p); %
            % S positive if cls0>cls1

        case 'ttest1side'
            [P,S]=ttest2_many_nan(dat,cls0,cls1,1,[],no_p);
        case 'ttest1side_rev'
            [P,S]=ttest2_many_nan(dat,cls1,cls0,1,[],no_p);
        case 'ranksum'
            [P,S]=ranksum_many(dat,cls0,cls1,0);

        case 'tnom'
            [P,S]=tnom_many_nan(dat,cls0,cls1);

        case 'ftest'  %F-test of variance (not Fisher's Exact!)
            v1=zeros(1,size(dat,2));
            v1(cls1)=1;
            v2=zeros(1,size(dat,2));
            if isfield(test_type,'confounding')
                v2=test_type.confounding;
            end
            [P,Br,msr,dfr,Bf,msf,dff]=ftest_many(dat,v1,v2);
            S=(msr./msf)';

        case 'snr'
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            m0=nanmean(dat(:,cls0),2);
            s0=nanstd(dat(:,cls0)')';
            m1=nanmean(dat(:,cls1),2);
            s1=nanstd(dat(:,cls1)')';
            S=(m1-m0)./(s0+s1);
        case 'gcsnr'
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            S=gcs2n(dat,cls0,cls1);

        case 'medmad_zprctile'
            m0=nanmedian(dat(:,cls0),2);
            s0=mad(dat(:,cls0),1,2)/norminv(0.75);
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            v=prctile(dat(:,cls1),test_type.prctile,2);
            v=(v-m0)./s0;
            S=prctile(v,test_type.prctile,2);

        case 'comp_oa'
            m0=nanmedian(dat(:,cls0),2);
            s0=mad(dat(:,cls0),1,2);
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            v=prctile(dat(:,cls1),test_type.prctile,2);
            v=(v-m0)./s0;
            S=v;

        case 'ort'
            p=prctile(dat(:,cls0),[25 75],2);
            cutoff=p(:,2)+diff(p,1,2);
            m1=median(dat(:,cls0),2);
            m2=median(dat(:,cls1),2);
            den=median([ abs(dat(:,cls0)-repmat(m1,1,length(cls0))) abs(dat(:,cls1)-repmat(m2,1,length(cls1)))],2);
            above=dat(:,cls1)>repmat(cutoff,1,length(cls1));
            num=sum(dat(:,cls1).*above,2);
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            den(den==0)=eps;
            S=num./den;

        case 'medmad_zprctile_all'
            m0=nanmedian(dat(:,[cls0 cls1]),2);
            s0=mad(dat(:,[cls0 cls1]),1,2)/norminv(0.75);
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            v=prctile(dat(:,cls1),test_type.prctile,2);
            v=(v-m0)./s0;
            S=prctile(v,test_type.prctile,2);

        case 'mediandiff'
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            v1=zeros(1,size(dat,2));
            v1(cls1)=1;
            if isfield(test_type,'confounding')
                v2=test_type.confounding;
            else
                v2=zeros(1,size(dat,2));
            end
            if ~isfield(test_type,'prctile0')
                test_type.prctile0=50;
            end
            if ~isfield(test_type,'prctile1')
                test_type.prctile1=50;
            end
            [u,ui,uj]=unique(v2);
            Sf=nan(size(dat,1),length(u));
            for i=1:length(u)
                ci=find(uj==i);
                c0=ci(find(v1(ci)==0));
                c1=ci(find(v1(ci)==1));
                if ~isempty(c0) && ~isempty(c1)
                    if test_type.prctile1==50
                        Sf(:,i)=nanmedian(dat(:,c1),2);
                    else
                        Sf(:,i)=prctile(dat(:,c1),test_type.prctile1,2);
                    end
                    if test_type.prctile0==50
                        Sf(:,i)=Sf(:,i)-nanmedian(prctile(dat(:,c0),2));
                    else
                        Sf(:,i)=Sf(:,i)-prctile(dat(:,c0),test_type.prctile0,2);
                    end
                end
            end
            S=nanmedian(Sf,2);

        case 'sumdiff'
            P=NaN*ones(size(dat,1),1); % no assymptotic P-values
            m0=sum(dat(:,cls0),2);
            m1=sum(dat(:,cls1),2);
            S=m1-m0;

        otherwise
            disp('No such method...using default');
            P=NaN*ones(size(dat,1),1);
            S=NaN*ones(size(dat,1),1);
            for i=1:size(dat,1)i
                dat0=dat(i,cls0);
                dat0(isnan(dat0))=[];
                dat1=dat(i,cls1);
                dat1(isnan(dat1))=[];
                if ~isempty(dat0) & ~isempty(dat1) & std([dat0 dat1])>0.001
                    [h,p,ci,stats]= ttest2(dat0,dat1,0.05,0);
                    P(i)=p;
                    S(i)=stats.tstat;
                else
                    P(i)=NaN;
                    S(i)=NaN;
                end
                if mod(i,1000)==0
                    verbose(num2str(i));
                end
            end
    end

    P_all(k+1:thisloopendidx) = P;
    S_all(k+1:thisloopendidx) = S;

    k = k+maxsnps;
end

