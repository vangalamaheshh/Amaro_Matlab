function X = add_newlines_to_long_strings(X,cutoff)

if ~exist('cutoff','var'), cutoff=20; end

if ~iscell(X), X={X}; end
nx = length(X);

for i=1:nx
  lines = split(X{i},char(10));
  len = cellfun('length',lines);
  for j=1:length(len)
    if len(j)>=cutoff
      z = find(lines{j}==' ' | lines{j}=='-');
      if isempty(z), z = 1:len(j); end  % have to split a word
      edgeness = abs(z-(len(j)/2));
      [tmp k] = min(edgeness);
      lines{j} = [lines{j}(1:z(k)) char(10) lines{j}(z(k)+1:end)];
    end
  end
  newname = cat(1,lines{:});
  newname = regexprep(newname,char(10),'\n');
  newname = regexprep(newname,['([a-zA-Z])' char(10)],'$1\-\n');
  newname = regexprep(newname,' \n','\n');
  X{i} = newname;
end
