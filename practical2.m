class_ids2=load_struct_noheader('~/Github/cs181-practicals/practical_2/ClassIDs.txt.txt');
class_ids=rename_field(class_ids,'col1','ClassNumericalID');
Key={'Test','Agent', 'AutoRun', 'FraudLoad', 'FraudPack', 'Hupigon', 'Krap',...
           'Lipler', 'Magania', 'None', 'Poison', 'Swizzor', 'Tdss',...
           'VB', 'Virut', 'Zbot'};

class_feature_matrix=load_struct_noheader('~/Github/cs181-practicals/practical_2/MatrixTrainCatTest.txt');
cf=struct2cell(class_feature_matrix);
for i=1:length(cf)
    data_array(i,:)=str2double(cf{i});
end
data_array=data_array';
[ix in]=sort(class_ids.ClassNumericalID);
imagesc(data_array(in,:))
hold on
current=ix(1);
for i=1:length(ix)
    if ~isequal(ix(i),current)
        current
        plot([0,116],[i,i],'r-','MarkerSize',.5)
        current=ix(i);
        
    end
end
current
for i=1:size(data_array,2)
n_pos(i)=sum(data_array(:,i)>1);
end
data_array_t=data_array(:,n_pos<4500);
test_a=bsxfun(@minus,data_array_t,mean(data_array_t,1));
test_a=(bsxfun(@rdivide,test_a,std(test_a)));
test_a=test_a(:,sum(~isnan(test_a))==6810);

prj=pca(test_a,30);
imagesc(prj(in,:))
filt_mat=tf_mat(:,sum(tf_mat,1)>10&sum(tf_mat,1)<50000);
matrix_for_pca=[test_a(1:3086,:),filt_mat];
matrix_for_pca=bsxfun(@minus,matrix_for_pca,mean(matrix_for_pca,1));
%matrix_for_pca=(bsxfun(@rdivide,matrix_for_pca,std(matrix_for_pca)));
matrix_for_pca=matrix_for_pca(:,sum(~isnan(matrix_for_pca))==3086);
prj_exp=pca(matrix_for_pca,30);
imagesc(prj_exp(in,:))