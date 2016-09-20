function k = map_categories_to_effect13(categs_txt)
% k = map_categories_to_effect13(categs_txt)
%
% categs_txt is path to file db/*/categs.txt listing the categories to be mapped
% collapses given categories to the standard list of 13 "effect" (silent/nonsilent) categories

demand_file(categs_txt);
C = load_struct(categs_txt);
demand_fields(C,{'num','name'});
if any(strcmp(C.num,'0')), error('%s has a category #0',categs_txt); end
if ~strcmp(C.num{1},'1'), error('%s first category is not #1',categs_txt); end

E = get_effect13_categories_list;
demand_field(E,'name');

k = nan(slength(C),1);
cidx = 1;
for ei=1:slength(E)
  idx = grep(E.name{ei},C.name,1);
  if isempty(idx), error('%s lacks a "%s" category!', categs_txt,E.name{ei}); end
  if any(~isnan(k(idx)))
    disp(C.name(idx(~isnan(k(idx)))))
    error('The above categories would be double counted!');
  end
  k(idx) = cidx;
  cidx=cidx+1;
end

if any(isnan(k(idx)))
  disp(C.name(idx(isnan(k(idx)))));
  fprintf('WARNING: the above categories were not counted.\n');
end



