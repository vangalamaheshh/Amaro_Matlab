function gp_plot(gctfile,rc,pd,dims,clsfile)

if ischar(rc) && ~isempty(str2num(rc))
  rc=str2num(rc);
end

if isnumeric(rc)  
    rc=enum_param(rc,{0,'row';1,'col'});
end

if ischar(pd) && ~isempty(str2num(pd))
  pd=str2num(pd);
end

if isnumeric(pd)  
    pd=enum_param(pd,{0,'dots';1,'profiles'});
end

if ischar(dims)
    dims=str2num(dims);
end

R=read_mit_gct_file(gctfile);

if exist('clsfile','var') && ~isempty(clsfile)
    R=read_mit_cls_file(R,clsfile,1);
    c=hsv2rgb([ repmat((0:(1/9):1)',3,1) [ ones(20,1); 0.5*ones(10,1)] ...
        [ones(10,1); 0.5*ones(10,1); ones(10,1);]]);
    c=c(mod((1:7:(size(c,1)*7))-1,size(c,1))+1,:);
    R=add_supmark(R,c);
    R.supmark(1).marker=cellstr(repmat('o',size(R.supmark(1).marker,1),1));
end

X=R.dat;
if is_col(rc)
    X=X';
end
X=X(:,dims);
if ~iscell(R.gacc)
  R.gacc=cellstr(R.gacc);
end
if ~iscell(R.gdesc)
  R.gdesc=cellstr(R.gdesc);
end
if ~iscell(R.sdesc)
  R.sdesc=cellstr(R.sdesc);
end


switch pd
 case 'dots'
  if length(dims)<2 || length(dims)>3
    error('Dots support only 2 or 3 dimension');
  end
  figure(1); clf;
  if isfield(R,'supdat')
    scatter_plot_D(X,R,1,1);
    [typeacc,typedesc,R1,range,non_empty]=decollapse_supdat(R,1);
    legend(deblank(cellstr(R1.supacc(range,:))));
  else
    scatter_plot_D(X,R);
  end
  fs=14;
  if is_col(rc)
    lh(1)=xlabel(R.gacc{1},'FontSize',fs);
    lh(2)=ylabel(R.gacc{2},'FontSize',fs);
    if length(dims)==3
      lh(3)=zlabel(R.gacc{3},'FontSize',fs);
    end
  else
    lh(1)=xlabel(R.sdesc{1},'FontSize',fs);
    lh(2)=ylabel(R.sdesc{2},'FontSize',fs);
    if length(dims)==3
      lh(3)=zlabel(R.sdesc{3},'FontSize',fs);
    end    
  end
  set(lh,'Interpreter','none');
  grid on
  
 case 'profiles'
  figure(1); clf;
  plot(1:length(dims),X,'o-'); %% HERE!!
  set(gca,'XTick',1:length(dims),'XTickLabel',cellstr(num2str((1:length(dims))')));
  if is_col(rc)
    lh=legend(R.sdesc);
  else
    lh=legend(R.gacc);
  end
  set(lh,'Interpreter','none');
end
set(gca,'Visible','On');
disp('end of plot');
fprintf(2,'stderr: end of plot');



% mcc -m gp_scatter_plot
