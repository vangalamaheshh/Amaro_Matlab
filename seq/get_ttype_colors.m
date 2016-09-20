function ttcolors = get_ttype_colors(ttnames)
% assign colors to tumor types

synonyms = {...
    'BLAD'   'Bladder'
    'BLCA'   'Bladder'
    'BR'     'Breast'
    'BRCA'   'Breast'
    'ESO'    'Esoph'
    'HN'     'HeadNeck'
    'HNSC'   'HeadNeck'
    'KIRC'   'KidneyRC'
    'KIRP'   'KidneyRP'
    'MED'    'Medullo'
    'MEL'    'Melanoma'
    'SKCM'   'Melanoma'
    'STAD'   'Stomach'
    'UCEC'   'Endometrial'
    'COAD'   'Colon'
    'READ'   'Rectal'
    'LAML'   'AML'
    'THCA'   'Thyroid'
    'CESC'   'Cervical'
    'CARC'   'Carcinoid'
    'PAAD'   'Pancreas'
    'PRAD'   'Prostate'
    'RHAB'   'Rhabdoid'
};

aliases=[]; aliases.old = synonyms(:,1); aliases.new = synonyms(:,2);
ttnames = apply_aliases(ttnames,aliases);

pos_names = {...
      'Combined', 'NB',      'Ewing',   'Rhabdoid', 'LUSC' ,   'LUAD',...
      'Thyroid',  'CLL',     'KidneyRP','KidneyRC', 'HeadNeck','Melanoma',...
      'Carcinoid','AML',     'Breast',  'OV',       'Cervical','Bladder',...
      'LAM',      'MM',      'Pancreas','Prostate', 'Esoph',   'DLBCL',...
      'Medullo',  'LGG',     'GBM',     'CRC',      'Stomach', 'Endometrial',...
      'Colon',    'Rectal',
  };

pos_colors = [...
      -1  -1  -1;    022 010 042;   040 025 050;   060 060 065;   100 000 000;   070 010 000;   ...
      063 065 082;   080 000 080;   019 015 021;   033 023 011;   007 021 080;   000 000 000;   ...
      044 011 000;   000 065 085;   100 050 050;   080 050 070;   100 050 000;   090 090 000;   ...
      033 023 011;   050 051 020;   033 023 011;   070 070 055;   030 080 018;   013 000 055;   ...
      009 019 029;   055 035 040;   054 044 036;   000 060 035;   000 100 000;   060 000 060;   ...
      008 055 035;   000 040 000; ...
  ];

ttnames_long = get_long_ttype_names(ttnames);
posnames_long = get_long_ttype_names(pos_names);
ord = listmap(ttnames_long,posnames_long);
ttcolors = nansub(pos_colors,ord)/100;

missing = find(isnan(ord));
if length(missing)>=0.8*length(ord)
  fprintf('Ignoring conventional colors\n');
  ttcolors = distinct_colors(length(ttnames));
elseif ~isempty(missing)
  fprintf('Using randomly generated colors for the following tumor types:\n');
  disp(ttnames(missing));
  random_colors = distinct_colors(length(ttnames)+10);
  ttcolors(missing,:) = random_colors(10:9+length(missing),:);
end


