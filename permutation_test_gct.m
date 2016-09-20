function permutation_test_gct
amps=load_table('GCT_Armlvl_data/amp_list_march.txt');
amps=rmfield(amps,{'header','headline'});
 CN.arms=amps.Arm;
%generate array
sample_ids=fieldnames(amps);
sample_ids={sample_ids{2:end}};
for i=1:length(sample_ids)
    for j=1:slength(amps)
        array(j,i)=amps.(sample_ids{i})(j);
    end
end
sums_per_patient=sum(array,1);
sums_per_arm_o=sum(array,2);

%permute
iter=100000;
sums_per_arm=zeros(iter,length(sums_per_arm_o));
for i=1:iter
    [a_s,ix(i,:)]=shuffle_gct_arms(array);
    sums_per_arm(i,:)=sum(a_s,2);
    if mod(10000,i)==0
        disp(sprintf('iteration %d',i))
    end
end
s=(reshape(sums_per_arm,[iter*39,1]));
test_a=[sums_per_arm_o;s];
test_a=test_a-mean(test_a(:));
test_a=test_a/std(test_a(:));
p=1-normcdf(test_a,0,1);
CN.p_amps=p(1:39);
CN.n_amps=sums_per_arm_o;
hist(s,26)
ylabel('Observations','FontSize',20)
xlabel('Counts','FontSize',20)
title('Permutation over Amps','FontSize',20)
dels=load_table('GCT_Armlvl_data/del_list_march.txt');
dels=rmfield(dels,{'header','headline'});
for i=1:length(sample_ids)
    for j=1:slength(dels)
        array(j,i)=dels.(sample_ids{i})(j);
    end
end
sums_per_patient=sum(array,1);
sums_per_arm_o=sum(array,2);
iter=100000;
sums_per_arm=zeros(iter,length(sums_per_arm_o));

for i=1:iter
    [a_s,ix(i,:)]=shuffle_gct_arms(array);
    sums_per_arm(i,:)=sum(a_s,2);
    if mod(10000,i)==0
        disp(sprintf('iteration %d',i))
    end
end
s=(reshape(sums_per_arm,[iter*39,1]));
test_a=[sums_per_arm_o;s];
test_a=test_a-mean(test_a(:));
test_a=test_a/std(test_a(:));
p=1-normcdf(test_a,0,1);
CN.p_dels=p(1:39);
CN.n_dels=sums_per_arm_o;
hist(s,22)

ylabel('Observations','FontSize',20)
xlabel('Counts','FontSize',20)
title('Permutation over Dels','FontSize',20)

cnloh=load_table('GCT_Armlvl_data/cnloh_list_march.txt');
cnloh=rmfield(cnloh,{'header','headline'});
for i=1:length(sample_ids)
    for j=1:slength(cnloh)
        array(j,i)=cnloh.(sample_ids{i})(j);
    end
end
sums_per_patient=sum(array,1);
sums_per_arm_o=sum(array,2);
iter=100000;
sums_per_arm=zeros(iter,length(sums_per_arm_o));

for i=1:iter
    [a_s,ix(i,:)]=shuffle_gct_arms(array);
    sums_per_arm(i,:)=sum(a_s,2);
    if mod(10000,i)==0
        disp(sprintf('iteration %d',i))
    end
end
s=(reshape(sums_per_arm,[iter*39,1]));
test_a=[sums_per_arm_o;s];
test_a=test_a-mean(test_a(:));
test_a=test_a/std(test_a(:));
p=1-normcdf(test_a,0,1);
CN.p_cnloh=p(1:39);
CN.n_cnloh=sums_per_arm_o;
hist(s,15)

ylabel('Observations','FontSize',20)
xlabel('Counts','FontSize',20)
title('Permutation over cnloh','FontSize',20)

end


