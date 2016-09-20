function out = collapse_1040_to_128(in)

if ~isnumeric(in), error('input should be numeric'); end

out = in(1:64,:)+ in(65:128,:);
