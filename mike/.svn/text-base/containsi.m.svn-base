function result = containsi(string,substring)
% case-insensitive version of contains

if ~ischar(string) | ~ischar(substring)
  error('containsi takes two strings as arguments');
end

result = contains(lower(string),lower(substring));
