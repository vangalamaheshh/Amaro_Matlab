function [M,repidx,repnames,repuniqid] = find_duplicate_samples(M,thresh,show_figs)
%FIND_DUPLICATE_SAMPLES use genomic data (affy_calls) to identify
%duplicated samples.
% 
%   [M,REPIDX,REPNAMES] = FIND_DUPLICATE_SAMPLES(M,THRESH,SHOW_FIGS) returns
%   data structure M with updated supdat and cell arrays, REPIDX and REPNAMES, in
%   which each cell lists the index or names of samples that are duplicates. M is
%   the data structure, and THRESH is the euclidian distance between
%   genotypes below which samples are considered identical.  (The distance
%   from AA to AB is 1.  The distance from AA to BB is 2.)  The supdat of M
%   is updated with integer identifiers such that unique samples receive
%   unique identifiers and replicated samples receive the same identifier
%   in the REP field.
%    (i.e. if M.affy_calls(:,10) and M.affy_calls(:,50) are genomically
%    indistinguishable according to the set thresh, the 10th
%   and 50th elements of LAB will be identical.)  SHOW_FIGS is optional
%   flag to select for showing figures.  FIND_DUPLICATE_SAMPLES prints the
%   duplicate sample list to 'identical_sets.txt'.  Genotype distance
%   histogram and genotype distance matrix are printed to
%   'Genotype_distance_dist' and 'Genotype_distance_matrix'.  
%   
%   Revisions:
%       12 Oct 07 -- Function added.  Jen Dobson (jdobson@broad.mit.edu)
%
%
%---
% $Id$
% $Date: 2008-08-25 02:18:42 -0400 (Mon, 25 Aug 2008) $
% $LastChangedBy: rameen $
% $Rev$

if ~iscell(M)
    tmp{1} = M;
    M = tmp;
    clear tmp;
    wascell = 0;
else
    wascell = 1;
end


if ~exist('thresh','var')
    thresh = 0.4; % See supp. material of Beroukhim, Getz et al (2007) PNAS
end

if ~exist('show_figs','var')
    show_figs = 0;
end

for k = 1:length(M)
    x = M{k}.affy_calls';
%    x(find(isnan(x)))=4;
   % x = single(x);
    Z = make_D(x,[],[],M{k}.sdesc);

        Z.dat = double(Z.dat);
    
    
    
    %exchange 2 and 3 so that AA=1;AB=2;BB=3
    
    Z.dat(Z.dat==2) = 5; %5 is dummy
    Z.dat(Z.dat==3) = 2;
    Z.dat(Z.dat==5) = 3;
    Z.dat(Z.dat==4)=NaN;
    Z.dat(Z.dat==0)=NaN;
   

    dd{k} = dist(Z.dat,[],'euclidean')/sqrt(size(Z.dat,2));  %make NXN distance matrix dd (N is number of samples)
    
    yy = tril(dd{k},-1)+triu(nan(size(dd{k})),0);  %use only lower triangular part
    
    if show_figs
        figure(1); clf;
        hist(yy(~isnan(yy)),500);   
        print_D('Genotype_distance_dist',{{'fig'},{'pdf'}});

        figure(2); clf;

        [Zord,Zord.dend]=one_way_clustering(Z,'rows',struct('cluster','average','distmat',dd{k},'is_sim',0));
        imagesc(dd{k} (Zord.dend.idx,Zord.dend.idx))
        print_D('Genotype_distance_matrix',{{'fig'},{'pdf'}});

    end
    
    y = tril(dd{k}<thresh,-1);
    
    [edg(:,1),edg(:,2),tmp]=find(y);  %edg = ordered pairs of replicated samples
    lab=unionfind(edg,size(y,1));  %labels of samples (replicates have same label and lowest numbers)
    [u,ui,uj]=unique(lab);

    h = histc(uj,1:length(uj));
    large=find(h>1);

    if isempty(large)
        islarge = zeros(size(lab));  %no replicated samples

    else
        islarge = lab<=max(large);  %these are the groups with greater than one sample

    end

    repids = islarge.*lab;



    M{k} = add_D_sup(M{k},'REP','Replicate',repids,'col');  %add 1 to supdat if duplicated


    if max(islarge) == 0
        repnames{k} = {};
        repidx{k} = {};
    else
        for ll = large
            idxs = find(repids == large(ll));
            repnames{k}{ll} = M{k}.sdesc(idxs); %the names of duplicated samples
            repidx{k}{ll} = idxs; 
        end
    end

% 
%     
%     f = fopen('identical_sets.txt','w');
%         
%     for i=1:length(large)
%         fprintf(f,'------------------------------\n');
%         for j=find(uj==large(i))
%             fprintf(f,'%s\n',Z.sdesc{j});
%         end
%     end
%     fclose(f);

end

    

if ~wascell
    tmp = repnames{1};
    repnames = tmp;
    tmp = M{1};
    M = tmp;
    tmp = repidx{1};
    repidx = tmp;
    
end

clear tmp

if isfield(M.sis,'uniqid')
  uniqids = cellfun(@str2num,M.sis.uniqid);
  for i = 1:length(repidx)
    repuniqid{i}=unique(uniqids(repidx{i}));
  end
else
  repuniqid = [];
end













    
