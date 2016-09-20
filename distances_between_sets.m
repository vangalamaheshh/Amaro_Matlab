function [av_dists p mean_germ_dis mean_spor_dis] = distances_between_sets(set1,set2,euc_dis1,nperm)
% 
% Given distance matrix euc_dis1, calculates average distances within
% set1, within set2, and between the 2 sets, and the one- and two-sided
% p-values for whether the distances differ significantly
%  
% Also returns the results of the permutations, to allow combined
% permutations across several groups. 
%
%---
% $Id$
% $Date: 2010-02-15 12:50:37 -0500 (Mon, 15 Feb 2010) $
% $LastChangedBy: rameen $
% $Rev$
  
sum_VHL_euc(1, [1 2]) = 0;
for i = 1:length(set1)
  for j = 1:i-1
    sum_VHL_euc(1,1) = sum_VHL_euc(1,1)+euc_dis1(set1(i),set1(j));
    sum_VHL_euc(1,2) = sum_VHL_euc(1,2)+1;
    within_VHL(set1(i),set1(j)) = 1;
  end
end
av_dists(1) = sum_VHL_euc(1,1)/sum_VHL_euc(1,2);

sum_sporadic_euc(1, [1 2]) = 0;
for i = 1:length(set2)
  for j = 1:i-1
    sum_sporadic_euc(1,1) = sum_sporadic_euc(1,1)+euc_dis1(set2(i),set2(j));
    sum_sporadic_euc(1,2) = sum_sporadic_euc(1,2)+1;
    within_sporadics(set2(i),set2(j)) = 1;
  end
end
av_dists(2) = sum_sporadic_euc(1,1)/sum_sporadic_euc(1,2);


sum_VHLsporadic_euc(1, [1 2]) = 0;
for i =1:length(set1)
  for j = 1:length(set2)
    sum_VHLsporadic_euc(1,1) = sum_VHLsporadic_euc(1,1)+euc_dis1(set1(i),set2(j));
    sum_VHLsporadic_euc(1,2) = sum_VHLsporadic_euc(1,2)+1;
    VHL_sporadics(set1(i),set2(j))=1;
  end
end
av_dists(3) = sum_VHLsporadic_euc(1,1)/sum_VHLsporadic_euc(1,2);

%%%%%% permutations to calculate p-values

if ~exist('nperm','var') || isempty(nperm)
  nperm = 10000;
end

total_set = [set1 set2];

mean_germ_dis=zeros(nperm,1);
mean_spor_dis=mean_germ_dis;
mean_germspor_dis=mean_germ_dis;

for k = 1:nperm
  temp_set = randperm(length(total_set));
  temp_germs = temp_set(1:length(set1));
  temp_spors = temp_set(length(set1)+1:length(set1)+length(set2));
  
  within_VHL = zeros(size(euc_dis1));
  within_sporadics = within_VHL;
  VHL_sporadics = within_VHL;
  
  sum_VHL_euc(1, [1 2]) = 0;
  for i = 1:length(temp_germs)
    for j = 1:i-1
      sum_VHL_euc(1,1) = sum_VHL_euc(1,1)+euc_dis1(temp_germs(i),temp_germs(j));
      sum_VHL_euc(1,2) = sum_VHL_euc(1,2)+1;
      within_VHL(temp_germs(i),temp_germs(j)) = 1;
    end
  end
  mean_germ_dis(k)=sum_VHL_euc(1,1)/sum_VHL_euc(1,2);
  
  sum_sporadic_euc(1, [1 2]) = 0;
  for i = 1:length(temp_spors)
    for j = 1:i-1
      sum_sporadic_euc(1,1) = sum_sporadic_euc(1,1)+euc_dis1(temp_spors(i),temp_spors(j));
      sum_sporadic_euc(1,2) = sum_sporadic_euc(1,2)+1;
      within_sporadics(temp_spors(i),temp_spors(j)) = 1;
    end
  end
  mean_spor_dis(k)=sum_sporadic_euc(1,1)/sum_sporadic_euc(1,2);
  
  diff_mean(k) = mean_spor_dis(k)-mean_germ_dis(k);  
end

p(1) = length(find(diff_mean>(av_dists(2)-av_dists(1))))/nperm;
p(2) = length(find(diff_mean>(av_dists(1)-av_dists(2))))/nperm;
p(3) = length(find(abs(diff_mean)>(abs(av_dists(2)-av_dists(1)))))/nperm;
