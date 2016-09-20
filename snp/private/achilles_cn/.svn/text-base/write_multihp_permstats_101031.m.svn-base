function n = write_multihp_permstats_101031(file,Ps,Qs,matcher,AC)

filter = ~all(all(isnan(Ps),2),3);

n = sum(filter);
[aci,cli] = find(matcher);
aci = aci(filter);
%cli = cli(filter);

npq = size(Ps,2);
pqcols = cell(2,npq,2);

for k = 1:npq
    kstr = num2str(k);
    pqcols{1,k,1} = {['lP' kstr],Ps(:,k,1)};
    pqcols{2,k,1} = {['lQ' kstr],Qs(:,k,1)};
    pqcols{1,k,2} = {['rP' kstr],Ps(:,k,1)};
    pqcols{2,k,2} = {['rQ' kstr],Qs(:,k,1)};
end

n = write_filtered_tabcols(file,[],...
               {'hairpin', AC.hairpinID(aci)},...
               {'gene', AC.geneID(aci)},...
               pqcols{:});
