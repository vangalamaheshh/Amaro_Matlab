function [categ_list categ_table] = get_effect13_categories_list

categ_list.num = (1:13)';
categ_list.name = {...
  'noncoding';                % 1
  'any change is silent';     % 2
  'change to A/C is silent';  % 3
  'change to A/G is silent';  % 4
  'change to A/T is silent';  % 5
  'change to C/G is silent';  % 6
  'change to C/T is silent';  % 7
  'change to G/T is silent';  % 8
  'change to A is silent';    % 9
  'change to C is silent';    % 10
  'change to G is silent';    % 11
  'change to T is silent';    % 12
  'any change is nonsilent';  % 13
};

if nargout>=2
  % nan = noncoding
  % 0 = silent
  % 1 = nonsilent
  %
  %  ->A  ->C  ->G  ->T
  categ_table = [...
     nan  nan  nan  nan ;...
      0    0    0    0  ;...
      0    0    1    1  ;...
      0    1    0    1  ;...
      0    1    1    0  ;...
      1    0    0    1  ;...
      1    0    1    0  ;...
      1    1    0    0  ;...
      0    1    1    1  ;...
      1    0    1    1  ;...
      1    1    0    1  ;...
      1    1    1    0  ;...
      1    1    1    1  ...
   ];     
end
