function [n_cluster,level, mode, member] = HMAC(data,sigmas)
% HMAC clusters 2-dimensional data based on HMAC algorithm described in [1]
% Jia Li, Surajit Ray & Bruce Lindsay. 
%
% %Inputs
% DATA which is a matrix whose rows contains the values for each
% observation.
% SIGMAS which contains the sequence of sigma used for
% kernel density estimate. 
%
% %Outputs
% N_CLUSTER is a vector of length(SIGMA) whose elements are number of
% cluster by corresponding sigma. 
% LEVEL is a vector of length(SIGMA) containing the level of clustering due
% to corresponding sigma.
% MODE is an structure that MODE.c1 includes modes of clusters at 1st level
% and MODE.c2 includes modes of clusters at 2nd level, and so on.
% MEMBER is a matrix (# of observations) by (# of levels) whose i-th
% column contains membership of each observation in i-th level. 
%                                                              
% Reference:
% [1] A Nonparametric Statistical Approach to Clustering via Mode   
% Identification, Jia Li, Surajit Ray and Bruce Lindsay, 2007
%
% Programmed by Yeojin Chung (ychung@psu.edu)                                                                                            %

[n,dim]=size(data);
s_hat=std(data);
n_cluster=[];G=data;member=1:n;member_n=[];
for j=1:length(sigmas)
    sigma=sigmas(j); %select j-th bandwidth
    M=ones(1,dim);g=1;clust=[];
    for i=1:size(G,1)
        x0=G(i,:);%set i-th row of G(cluster representative) as an initial point
        x0_old=1;d=1;
        while d>10^(-5)
            f=mvnpdf(x0,data,diag(sigma^2*ones(dim,1)));
            p=f/sum(f); %update p
            x0=p'*data; %update x
            if norm(x0)~=0
                d=norm(x0_old-x0)/norm(x0); %for stopping rule
            else d=10;
            end
            x0_old=x0;
        end
        %assign obs which lead to the same mode into a cluster 
        %clust:assigned cluster to each row of G. # of rows(clust)==# of
        %rows(G)
        temp=find(all([abs(M-ones(size(M,1),1)*x0)<ones(size(M,1),1)*s_hat*0.001],2));
        if length(temp)==0
            clust(i)=g;g=g+1;
        else
            if isempty(clust)
                clust=1;
            end
            clust(i)=min(clust(temp));
        end
        %M:all collection of modes with j-th bandwidth
        M(i,:)=x0;
        %member_n: assigned cluster to each observation of data by
        %          current bandwidth
        %member:   "   by previous bandwidth
        member_n(find(member==i))=clust(i);
    end
    %save # of clusters by each bandwidth
    n_cluster(j)=max(clust);
    %update the cluster representative, G
    [b,m,m1]=unique(clust);
    G=M(m,:);
    %save the variables
    eval(['center',int2str(j),'=G;'])
    eval(['member',int2str(j),'=member_n;'])
    %update membership 
    member=member_n;
end
[a1,a2,a3]=unique(n_cluster);
level=length(a1)-a3+1;
a1=a1(length(a1):-1:1);a2=a2(length(a2):-1:1);
mode=[];member=[];
for i=1:length(a1)
    eval(['mode.c',int2str(i),'=center',int2str(a2(i)),';'])
    member=[member; eval(['member',int2str(a2(i))])];
end
member=member';
        
