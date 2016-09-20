function [ridges,S,sep,density]=REM(data, mode,sigma, member, alpha)
% REM calculates ridgelines and separabilities between each pair of
% clusters. The outputs of HMAC can be used for input of this function.
%
% %Inputs
% DATA is a data matrix used for clustering. 
% MODE is a matrix whose rows are mode of corresponding cluster. 
% # of columns(DATA)==# of columns (MODE)
% SIGMA is a scalar used as a bandwidth in clustering.
% MEMBER is a vector whose length is the same as # of rows in DATA. Each
% element of MEMBER shows assigned cluster of corresponding observation.
% The maximum elements of MEMBER should be the same as # of rows of MODE.
% ALPHA is a sequence between 0 and 1. It is used as a grid of ridgelines.
%
% %Outputs
% RIDGES is a structure containing ridgelines. RIDGE.c1c2 is a matrix of
% ridgeline between cluster 1 and 2, and RIDGE.c1c3 is a ridgeline between
% cluster1 and 3, and so on.
% S is a (max(MEMBER))by(max(MEMBER)) matrix whose (i,j) element is a
% separability between i-th and j-th cluster. 
% SEP is a matrix whose first column contains the size of each cluster and
% the second column contains separabilities.
% DENSITY is a structure containing the densities of each ridgelines.
% 
% Reference:
% [1] A Nonparametric Statistical Approach to Clustering via Mode   
% Identification, Jia Li, Surajit Ray and Bruce Lindsay, 2007
%
% Programmed by Yeojin Chung (ychung@psu.edu)                                                                                            %

if size(mode,1)~=max(member) 
    'The number of rows of MODE should be the same with the maximum element in MEMBER.'
else
    n_a=length(alpha);[n,m]=size(data);
    S=[];ridges=[];density=[];
    for i=1:(size(mode,1)-1)
        for j=(i+1):size(mode,1)
            I1=find(member==i);I2=find(member==j);
            size1=size(I1,1);size2=size(I2,1);%get indexes of observations in i-th and j-th cluster 
            x=mode(i,:);%starting point of REM of first ridgeline
            ridge=[];g1=[];g2=[];
            for a=alpha
                d=1;fv_o=0;
                while d>10^(-5)
                    pdf1=mvnpdf(x,data(I1,:),diag(sigma^2*ones(m,1)));
                    pdf2=mvnpdf(x,data(I2,:),diag(sigma^2*ones(m,1)));
                    %update p's
                    p1=pdf1/sum(pdf1);
                    p2=pdf2/sum(pdf2);
                    %update x
                    x=(1-a)*p1'*data(I1,:)+a*p2'*data(I2,:);
                    fv=(sum(pdf1)/size1)^(1-a)*(sum(pdf2)/size2)^a;
                    if norm(fv)~=0
                        d=norm(fv_o-fv)/norm(fv);
                    else d=10;
                    end;
                    fv_o=fv;
                end
                ridge=[ridge;x]; %save the ridgeline 
                %save g_i and g_j for calculating pairwise separablility
                g1=[g1;sum(pdf1)/size1];
                g2=[g2;sum(pdf2)/size2];
            end
            eval(['ridges.c',int2str(i),'c',int2str(j),'=ridge;']) %save ridgelines
            eval(['density.c',int2str(i),'c',int2str(j),'=(size1*g1+size2*g2)/(size1+size2);'])
            %calculate pairwise separiability between cluster i and j
            S(i,j)=1-min(size1*g1+size2*g2)/(size1*g1(1)+size2*g2(1));
            S(j,i)=1-min(size1*g1+size2*g2)/(size1*g1(n_a)+size2*g2(n_a));
        end
        S(i,i)=1;S(j,j)=1; %put 1 into diagnals for calculating overall separability
    end

    % calculate overall separabilities
    sep=[];
    for i=1:size(mode,1)
        sep=[sep;length(find(member==i)),min(S(i,:))];
    end
end
