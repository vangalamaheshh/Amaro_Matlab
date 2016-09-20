function res=classify_bel_cor(res,method_name,delta,q,Jij,labels)
% classify based on beliefs and cor
% 1..q classes
% 0 can't decide
% -1... -Inf new class
%
% labels: NaN for unknown or don't test

bel=getfield(res,method_name,'bel');
cor=getfield(res,method_name,'cor');

n=length(labels);
known=find(~isnan(labels));
has_labels=min(find(Jij(:,1)>q));
lbled=Jij(1:(has_labels-1),2)-q;
known_no_lbl=setdiff(known,lbled);

% change Noam
spcJij=Jij(has_labels:end,:);
spcJij(:,1:2)=spcJij(:,1:2)-q;

% add Noam
if length(unique(labels))>q
  Q=length(unique(labels));
else
  Q=q;
end

for i=1:length(cor) % bel can be empty
  %if i==41
  %  keyboard
  %end
  
  if ~isempty(bel)
        bb=bel{i};
        [tmp,cl]=max(bb);
        bb(cl+(0:(size(bb,2)-1))*size(bb,1))=-1;
        [tmp2,cl2]=max(bb);
        not_sure=(tmp-tmp2<delta);
    else
        cl=ones(1,n);
        not_sure=cl;
    end
    cannot_decide=find(not_sure);    

    % add part on finding new classes !!!!!
    conf=cell(n,1);
    cl2=cl;
    if ~isempty(cannot_decide)
        cc=cor{i};
        %keyboard
        edg=spcJij(find(cc>0.5*(1+1/q)),:);
        labs=unionfind(edg,n);
        new_class=q+1;
        for j=1:max(labs)
            undecide_in_lab=intersect(cannot_decide,find(labs==j));
            decide_in_lab=intersect(setdiff(1:n,cannot_decide),find(labs==j));
            if ~isempty(undecide_in_lab)  
                if ~isempty(decide_in_lab)
                    type_hist=hist(cl(decide_in_lab),1:q);
                    [mt,mti]=max(type_hist);
                    if mt==sum(type_hist)
                        cl2(undecide_in_lab)=mti;
                    else
                        cl2(undecide_in_lab)=-1; %confused
                        for k=1:length(undecide_in_lab)
                            conf{undecide_in_lab(k)}=type_hist;
                        end
                    end
                else
                    cl2(undecide_in_lab)=new_class;
                    new_class=new_class+1;
                end
            end
        end
    end
    
    sols{i}.cl2=cl2;
    sols{i}.conf=conf;
    sols{i}.not_sure=not_sure;
    sols{i}.cl=cl;
    sols{i}.delta=delta;
    
    if (0)
        known_no_lbl_decide=setdiff(known_no_lbl,cannot_decide);
        sols{i}.n_cannot_decide=length(cannot_decide);
        tmp=getfield(res,method_name,'Ts');
        sols{i}.T=tmp(i);
        
        % 1) performance on all
        % 2) performance on sure
        % 3) best performance on not confused
        
        % 1) 
        test_on=known_no_lbl;
        tp=cl;
        %keyboard
        [p1,a1]=purity_accuracy2(tp(test_on),labels(test_on),q);
        e=calc_configuration_energy([1:q tp],Jij,q,has_labels);
        f_half=(2*p1*a1)/(p1+a1);
        cm=zeros(Q,Q);
        cm1=crosstab(tp(test_on),labels(test_on));
        cm(1:size(cm1,1),1:size(cm1,2))=cm1;
        sols{i}.perf{1}.cm=cm;
        sols{i}.perf{1}.cl=tp;
        sols{i}.perf{1}.f=f_half;
        sols{i}.perf{1}.n_test_on=length(test_on);
        sols{i}.perf{1}.e=e;
        sols{i}.perf{1}.nerr=sum(sum(cm-diag(diag(cm))));
        
        % 2) 
        test_on=setdiff(known_no_lbl,find(not_sure==1));
        tp=cl;
        [p1,a1]=purity_accuracy(tp(test_on),labels(test_on),q);
        e=calc_configuration_energy([1:q tp],Jij,q,has_labels);
        f_half=(2*p1*a1)/(p1+a1);
        cm=zeros(Q,Q);
        cm1=crosstab(tp(test_on),labels(test_on));
        cm(1:size(cm1,1),1:size(cm1,2))=cm1;
        sols{i}.perf{2}.cm=cm;
        sols{i}.perf{2}.cl=tp;
        sols{i}.perf{2}.f=f_half;
        sols{i}.perf{2}.n_test_on=length(test_on);
        sols{i}.perf{2}.e=e;
        sols{i}.perf{2}.nerr=sum(sum(cm-diag(diag(cm))));
        
        % 3) 
        test_on=setdiff(known_no_lbl,find(cl2==-1));
        tp=cl2;
        LABS=labels;
        LABS(isnan(labels))=0;
        for j=(q+1):max(cl2)
            new_class_j=find(cl2==j);
            h1=hist(LABS(new_class_j),0:q);
            [tmp,mxi]=max(h1(2:end));
            tp(new_class_j)=mxi;
        end
        [p1,a1]=purity_accuracy(tp(test_on),labels(test_on),q);
        e=calc_configuration_energy([1:q tp],Jij,q,has_labels);
        f_half=(2*p1*a1)/(p1+a1);
        cm=zeros(Q,Q);
        cm1=crosstab(tp(test_on),labels(test_on));
        cm(1:size(cm1,1),1:size(cm1,2))=cm1;
        sols{i}.perf{3}.cm=cm;
        sols{i}.perf{3}.cl=tp;
        sols{i}.perf{3}.f=f_half;
        sols{i}.perf{3}.n_test_on=length(test_on);
        sols{i}.perf{3}.e=e;
        sols{i}.perf{3}.nerr=sum(sum(cm-diag(diag(cm))));
    end
end

res=setfield(res,method_name,'sol',sols);


