function [r,fet,s,p,fpr,fwer,rankp,fdr,q,pi0,T,idx]=gp_marker_selection(D,supid,test_type,nperm)

global my_gp

if ~isempty(my_gp)  
  write_mit_res_file('temp.res',D);
  supdat2cls(D,supid,'temp'); % write_mit_cls_file('temp.cls',D,supid);
  if ischar(test_type)
    test_type.method=test_type;
  end
  gp_params=cell2struct({'temp.res',['temp.' deblank(D.supacc(supid,:)) '.cls'],test_type.method,num2str(nperm)},...
                        {'input_filename','cls_filename','test_statistic','number_of_permutations'},2);
  if isfield(test_type,'params')
    gp_params=add_struct(gp_params,test_type.params);
  end
  res=runAnalysis(my_gp,'ComparativeMarkerSelection',gp_params);
  
  [r,fet,s,p,fpr,fwer,rankp,fdr,q,pi0,T]= ...
      read_mit_marker_selection_results('temp.comp.marker.odf');
  idx=cell2mat(findstrings_list(strvcat(D.gacc),fet));
  r(idx)=r;
  fet(idx,:)=fet;
  s(idx)=s;
  p(idx)=p;
  fpr(idx)=fpr;
  fwer(idx)=fwer;
  rankp(idx)=rankp;
  fdr(idx)=fdr;
  q(idx)=q;
else
  error('no my_gp');
end

