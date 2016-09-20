function X = load_and_process_isz_data(iszdir,lanetable,params)
% load_and_process_isz_data(iszdir,lanetable,params)
%
% Process InsertSizeByLane data
%
% Mike Lawrence 2009-06-25

if ~exist('params','var'), params=[]; end

fprintf('Chromosomes: ');for c=1:24, fprintf('%d/24 ',c);
  tmp = load_matrix([iszdir '/chr' num2str(c) '.isz']);
  if c==1, X.raw = tmp(:,2:end); else X.raw(:,:,c) = tmp(:,2:end); end
end, fprintf('\n');
X.dat = sum(X.raw,3);

X.lanes = load_struct(lanetable);
X.lanes = make_numeric(X.lanes,'lane');

X = process_isz_data(X,params);
