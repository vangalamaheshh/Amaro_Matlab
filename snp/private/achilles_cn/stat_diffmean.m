function stats = stat_diffmean(x,y)
%STAT_DIFFMEAN Difference of the mean statistic
stats = (mean(x,1)-mean(y,1));