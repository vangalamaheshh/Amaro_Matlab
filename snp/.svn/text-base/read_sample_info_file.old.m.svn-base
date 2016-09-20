function s=read_sample_info_file(fname);

cols={'array','name','type','cellxeno','primarymet','hormone','wga','gender',...
      'ploidy','lohctrl','dup','core','good','past_qc','rep','in100','loh','platform','batch',...
      'autopsy','age','survival','histology'};
% 'exp1','exp2',...
s=[];
d=read_dlm_file(fname);
for j=1:length(cols)
  disp([ d{1}{j} '-->' cols{j}]);
end
for i=2:length(d)
  for j=1:length(cols)
    if j<=length(d{i})
      s=setfield(s,{i-1},cols{j},fill_empty(d{i}{j},'EMPTY'));
    else
      s=setfield(s,{i-1},cols{j},'EMPTY');
    end
  end
  if mod(i,100)==0
    disp(i);
  end
end

%   s(i-1).name=fill_empty(d{i}{2},'EMPTY');
%   s(i-1).type=fill_empty(d{i}{3},'EMPTY');
%   s(i-1).cellxeno=fill_empty(d{i}{4},'EMPTY');
%   if length(d{i})>=7
%     s(i-1).wga=fill_empty(d{i}{7},'EMPTY');
%   else
%     s(i-1).wga='EMPTY';
%   end
%   if length(d{i})>=9
%     s(i-1).ploidy=fill_empty(d{i}{9},'EMPTY');
%   else
%     s(i-1).ploidy='EMPTY';
%   end
% end

