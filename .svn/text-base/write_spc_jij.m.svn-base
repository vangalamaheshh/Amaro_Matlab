function write_spc_jij(prefix,J,TrainLabels)

if exist('TrainLabels','var')
  n_train=size(TrainLabels,1);
  q=max(TrainLabels(:,2));
  J(:,1:2)=J(:,1:2)+q;
  J=sortrows([ TrainLabels(:,2) TrainLabels(:,1)+q repmat(1000,n_train,1); J]);
end
    

f=fopen([ prefix '.Jij'],'w');
for j=1:size(J,1)
  fprintf(f,'%d %d %f\n',J(j,1),J(j,2),J(j,3));
end
fclose(f);

f=fopen([ prefix '.ij'],'w');
for j=1:size(J,1)
  fprintf(f,' %d %d\n',J(j,1),J(j,2));
end
fclose(f);
