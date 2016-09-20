function gp_clustering_figure(varargin)
% gp_clustering_figure -o xxx.png -ot 1 -d xxx.gct -dp 1 -p Tissue Type -c xxx.cls -cc
% color_table.txt -gd xxx.gdend.odf -gdr 1 -sd xxx.sdend.odf -sdr 1
% -sf 5 -gf 10 -sx 0.1 0.3 0.4 -sy 0.2 0.1 0.1 0.4 0.1 -l 1 -lnc 1 -ltf 14 -lf 14 

a=handle_args({'o','ot','d','dp','p','c','cc','gd','gdr','sd','sdr','sf','gf','sx','sy','l','lnc','ltf','lf'},varargin);

if ~isempty(a.d)
  D=read_mit_gct_file(a.d{1});
  if ~isempty(a.dp)
    data_preprocessing=enum_param(a.dp{1},...
                                  {0,'none'; 1, ...
                        'row_centerd_and_normalize';2,'col_centered_and_normalize'});
  end
end


if ~isempty(a.c)
  D=read_mit_cls_file(D,a.c{1},1);
  if ~isempty(a.p)
    phen=sprintf('%s ',a.p{:});
  else
    phen='';
  end
  pos=find(D.supacc(1,:)==':');
  if ~isempty(pos)
    D.supacc=[ phen ':' D.supacc(1,(pos+1):end)];
  end
  if ~isempty(a.cc)
    c=read_colorscheme(a.cc{1});
  else
    c=hsv2rgb([ repmat((0:(1/9):1)',3,1) [ ones(20,1); 0.5*ones(10,1) ...
                   ] [ones(10,1); 0.5*ones(10,1); ones(10,1);]]);
  end
  D.supmark(1).colormap=c;
  D.supmark(1).height=1;
  D.supmark(1).patchwidth=1;      
  D.supmark(1).linewidth=0;
end

if ~isempty(a.gd)
  gdend=read_mit_dend_file(a.gd{1});
  if ~isempty(a.gdr)
    gdend_reduced=enum_param(a.gdr{1},{0,'all';1,'reduced_big'});
  else
    gdend_reduced='all';
  end
else
  gdend=[];
  gdend_reduced='';
end

if ~isempty(a.sd)
  sdend=read_mit_dend_file(a.sd{1});
  if isfield(D,'supdat')
    D.supdat=D.supdat(:,sdend.idx);
  end
  if ~isempty(a.sdr)
    sdend_reduced=enum_param(a.sdr{1},{0,'all';1,'reduced_big'});    
  else
    sdend_reduced='all';
  end
else
  sdend=[];
  sdend_reduced='';
end

if ~isempty(a.sf)
  sf=str2num(a.sf{1})';
else
  sf=5;
end

if ~isempty(a.gf)
  gf=str2num(a.gf{1})';
else
  gf=5;
end

if ~isempty(a.l)
  do_landscape=str2num(a.l{1})';
else
  do_landscape=1;
end

if ~isempty(a.lnc)
  leg_n_col=str2num(a.lnc{1})';
else
  leg_n_col=2;
end

if ~isempty(a.ltf)
  leg_title_f=str2num(a.ltf{1})';
else
  leg_title_f=8;
end

if ~isempty(a.lf)
  leg_f=str2num(a.lf{1})';
else
  leg_f=6;
end

if ~isempty(a.sx)
  sx=str2num(strvcat(a.sx))';
else
  sx=[1 1 1];
end

if ~isempty(a.sy)
  sy=str2num(strvcat(a.sy))';
else
  sy=[1 1 1 1 1];
end

sx2=sx;
sx2(sx==0)=eps;
disp_p.x.sizes=sx2;
disp_p.x.gaps=ones(1,length(sx)+1);
disp_p.x.gaps(sx==0)=0;
disp_p.x.border=0.05;

sy2=sy;
sy2(sy==0)=eps;
disp_p.y.sizes=sy2;
disp_p.y.gaps=ones(1,length(sy)+1);
disp_p.y.gaps(sy==0)=0;
disp_p.y.border=0.05;

% 3x5
full_items={{1:3,1:2,'ssuplegend',leg_n_col,leg_title_f,leg_f},...
            {1,3,'sdend',sdend_reduced,1},...
            {2,3,'ssupdat'},...
            {3,3,'sdesc','horizontal',sf,[],0,[],{'HorizontalAlignment','Right'}},...
            {4,3,'data','row_center_and_normalize'},...
            {4,1,'gdend',gdend_reduced,1},...
            {4,2,'gacc','vert',gf},...
            {5,3,'colorbar','horizontal',...
                   [1 11.5 22 32.5 43 53.5 64],...
                   {'-3\sigma','-2\sigma','-1\sigma','0',['+1\' ...
                    'sigma'],'+2\sigma','+3\sigma'},6}};

disp_p.items=full_items;
item_names=get_cell_elem(disp_p.items,3); 
item_names=cat(1,item_names{:});

% no gdend
if isempty(gdend)
  i=strmatch('gdend',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
end

% no sdend
if isempty(sdend)
  i=strmatch('sdend',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
end


if (sx(2)==0) || ~isfield(D,'gacc') % no gacc
  i=strmatch('gacc',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
end

if (sx(3)==0) || (~isfield(D,'dat') && ~isempty(gdend)) % no data, sdend, sdesc, colorbar, ssupdat, ssuplegend
    [i1,i2,i3]=find(match_string_sets({'data','sdend','sdesc',...
            'colorbar','ssupdat','ssuplegend'},item_names));
    disp_p.items(i2)=[];
    item_names(i2)=[];
end
    
if (sy(4)==0) || (~isfield(D,'dat') && ~isempty(sdend)) % no data, gdend, gacc, colorbar
    [i1,i2,i3]=find(match_string_sets({'data','gdend','gacc',...
            'colorbar'},item_names));
    disp_p.items(i2)=[];
    item_names(i2)=[];
end

if (sy(2)==0) || ~isfield(D,'supdat') % no ssupdat ssuplegend
  i=strmatch('ssupdat',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
  i=strmatch('ssuplegend',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
end

if (sy(3)==0) || ~isfield(D.sdesc) % no sdesc
  i=strmatch('sdesc',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
end

if sy(5)==0 % no colorbar 
  i=strmatch('colorbar',item_names);
  disp_p.items(i)=[];
  item_names(i)=[];
end

fname=a.o{1};
output_type={'png','-r180'};
if ~isempty(a.ot)
    if str2num(a.ot{1})==1
        output_type={'epsc'};
        if length(fname)<5 ||  ~strcmp(fname((end-2):end),'eps')
          fname=[fname '.eps'];
        end
    end
end
if strcmp(output_type{1},'png')
  if length(fname)<5 || ~strcmp(fname((end-2):end),'png')
    fname=[fname '.png'];
  end
end
  
close all
figure(1);
res=display_D(D,gdend,sdend,disp_p);
try
  cur_fig=get(1);
  print_D(fname,{output_type},do_landscape,[]);
  close all
catch
  close all
  figure(1);
  res=display_D(D,gdend,sdend,disp_p);
  try
    cur_fig=get(1);
    print_D(fname,{output_type},do_landscape,[]);
    close all
  catch
    error('Cannot generate figure.... try again');
  end
end

