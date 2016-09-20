function gp_classifier_xvalidation(gctfile,clsfile,temp_pref,outfile,...
    test_type,min_std,fets,sig_type,sigs,prior_type,prior,dist_type,matrix_file)

% read gctfile
D=read_mit_gct_file(gctfile);

% read clsfile
D=read_mit_cls_file(D,clsfile,1);

if ischar(test_type)
    test_type=str2num(test_type);
end

if ischar(min_std)
    min_std=str2num(min_std);
end

if ischar(fets)
    fets=colon2vals(fets);
end

if ischar(sigs)
    sigs=colon2vals(sigs);
end

if ischar(sig_type)
    sig_type=str2num(sig_type);
end

if ischar(prior_type)
    prior_type=str2num(prior_type);
end  

if ischar(prior)
    prior=str2num(prior);
end

dist_type_num=str2num(dist_type);
if ~isempty(dist_type_num)
    switch(dist_type_num)
        case 0
            dist_type='euclidean';
        case 1
            dist_type='cosine';
        otherwise
            error('No such distance type.');
    end
end

N=size(D.dat,2);
n_cls=D_n_cls(D,1);

n_fets=length(fets);
n_sigs=length(sigs);

nerr=zeros(n_cls,N,n_fets,n_sigs);
fp=zeros(n_cls,N,n_fets,n_sigs);
fn=zeros(n_cls,N,n_fets,n_sigs);
logscore=zeros(n_cls,N,n_fets,n_sigs);
use_files=0;
% LOO
global FET_SEL_SI
FET_SEL_SI=cell(N,n_cls);
for looi=1:N
    Dtrain=reorder_D_cols(D,setdiff(1:N,looi));
    Dtest=reorder_D_cols(D,looi);
    Dtrain.fet_sel_id=looi;
    for fi=1:n_fets
        f=fets(fi);
        for si=1:n_sigs
            s=sigs(si);
            if (use_files)
                % temp_fname=[temp_pref '.' num2str(looi) '.' num2str(fi) '.' num2str(si) ];
                temp_fname=[temp_pref '.xvtemp'];
                write_mit_gct_file([temp_fname '.train.gct'],Dtrain);
                write_mit_cls_file([temp_fname '.train.cls'],Dtrain,1,1);
                write_mit_gct_file([temp_fname '.test.gct'],Dtest);
                write_mit_cls_file([temp_fname '.test.cls'],Dtest,1,1);
                
                gp_pnn([temp_fname '.train.gct'],[temp_fname '.train.cls'],...
                    [temp_fname '.test.gct'],[temp_fname '.test.cls'],...
                    [temp_fname '.full.odf'],[temp_fname '.pred.odf'],f,...
                    test_type,min_std,sig_type,s,prior_type,prior,dist_type,0,'');
                res=read_mit_odf_file([temp_fname '.full.odf']);                
            else
                [res,fullres]=gp_pnn(Dtrain,[],Dtest,[],[],[],f,...
                    test_type,min_std,sig_type,s,prior_type,prior,dist_type,0,'');
            end 
            nerr(find(fullres.data{8}==0),looi,fi,si)=1;
            fp(find(str2num(strvcat(fullres.data{4}))==0 & ...
                    str2num(strvcat(fullres.data{5}))==1),looi,fi,si)=1;
            fn(find(str2num(strvcat(fullres.data{4}))==1 & ...
                    str2num(strvcat(fullres.data{5}))==0),looi,fi,si)=1;            
            logscore(:,looi,fi,si)=fullres.data{7};
            
        end
    end
    fprintf(1,'%d ',looi);
end
fprintf(1,'\r\n');
FET_SEL_SI={};

% write chosen parameters
ne=squeeze(sum(nerr,2)); % sum over LOO
tfp=squeeze(sum(fp,2));
tfn=squeeze(sum(fn,2));
mls=squeeze(mean(logscore,2));

if exist('matrix_file','var')
  save(matrix_file,'ne','tfp','tfn','mls');
end
% save last_point_x

for c=1:n_cls
  ne_c=squeeze(ne(c,:,:));
  tfp_c=squeeze(tfp(c,:,:));
  tfn_c=squeeze(tfn(c,:,:));
  mls_c=squeeze(mls(c,:,:));
  [min_err,chosen_i]=min(ne_c(:));
  with_min_err=find(ne_c(:)==min_err);
  if length(with_min_err)>1
    [max_logscore,max_i]=max(mls_c(with_min_err));
    chosen_i=with_min_err(max_i);
  else
    max_logscore=mls_c(chosen_i);
  end
  [chosen_fet_i,chosen_sig_i]=ind2sub(size(ne_c),chosen_i);
  chosen_params.data{1}(c)=c;
  chosen_params.data{2}(c)=fets(chosen_fet_i);
  chosen_params.data{3}(c)=test_type;
  chosen_params.data{4}(c)=min_std;
  chosen_params.data{5}(c)=sig_type;
  chosen_params.data{6}(c)=sigs(chosen_sig_i);
  chosen_params.data{7}(c)=prior_type;
  chosen_params.data{8}(c)=prior;
  chosen_params.data{9}{c}=dist_type;
  chosen_params.data{10}(c)=min_err;
  chosen_params.data{11}(c)=tfp_c(chosen_i);
  chosen_params.data{12}(c)=tfn_c(chosen_i);  
  chosen_params.data{13}(c)=max_logscore;
end  
chosen_params.col_names={'Class' 'N Features' 'Test type' ...
                    'Minimal std.' 'Sigma type' 'Sigma' ...
                    'Prior type' 'Prior' 'Distance type' ...
                    'CV num. errors' 'CV false positives' ...
                    'CV false negatives' 'CV mean logscore'};
chosen_params.col_types={'float','float','float','float',...
                    'float','float','float','float','string', ...
                    'float','float','float','float'};
chosen_params.Model='PNNXValidation Optimization';
write_mit_odf_file(outfile,chosen_params);                    

