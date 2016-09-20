function T = print_mutation_list(M,is_germline)

if ~exist('is_germline','var') || isempty(is_germline)
  is_germline = (size(M,2)==9);
end

bases = 'ACGTN';
models={'AA','AC','AG','AT','CC','GG','TT','CG','CT','GT'};
T = [];
for i=1:size(M,1)
  if is_germline
    T = [T sprintf('%d\t%d\t%s\t%s\t%.2f\t%.2f\t%.2f\t%d\t%.2f\n',...
                   M(i,1:2),bases(M(i,3)),models{M(i,4)},M(i,5:9))];
  else
    T = [T sprintf('%d\t%d\t%s\t%s\t%.2f\t%.2f\t%.2f/%.2f\t%.2f/%.2f\t%d/%d\t%.2f\n',...
                   M(i,1:2),bases(M(i,3)),bases(M(i,4)),M(i,5:13))];
  end
end
