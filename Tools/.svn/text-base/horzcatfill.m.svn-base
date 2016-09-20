function out = horzcatfill(varargin)

sizes = cellfun(@size,varargin,repmat({1},1,length(varargin)));
endsizes = max(sizes).*ones(1,length(sizes));
if ~isequal(sizes,endsizes)

    varargin = cellfun(@emptyfill,varargin,mat2cell(endsizes,1,ones(1,length(endsizes))),'UniformOutput',0);

end

out = horzcat(varargin);
    

function fill = emptyfill(x,l)

if isnumeric(x(1)) || ischar(x(1))
    fill = vertcat(x,cast(zeros(l-size(x,1),1),class(x(1))));
elseif iscell(x(1)) 
    fill = vertcat(x,repmat({[]},l-size(x,1),1));
elseif isstruct(x(1))
    emptystruct = x(1);
    for fl = fieldnames(emptystruct)'
        emptystruct.(char(fl)) = '';
    end
    fill = vertcat(x,repmat(emptystruct,l-size(x,1),1));
end

    