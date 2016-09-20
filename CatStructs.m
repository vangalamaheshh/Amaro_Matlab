function S=CatStructs(struct1, struct2)
fields1=fieldnames(struct1);
fields2=fieldnames(struct2);
if isequal(fields1,fields2)
S=concatStructs_sameFields(struct1,struct2,fields1,fields2);
else
    nOV1=fields1(~ismember(fields1,fields2));
    nOV2=fields2(~ismember(fields2,fields1));
    CommonFields=fields1(ismember(fields1,fields2));
    S=concatStructs_differentFields(struct1,struct2,CommonFields,nOV1,nOV2);
end
end



function S=concatStructs_sameFields(struct1,struct2,fields1,fields2)
for i=1:slength(struct1)
    for j=1:length(fields1)
        if iscell(struct1.(fields1{j})(i))
        S.(fields1{j}){i,1}=struct1.(fields1{j}){i};
        else
                S.(fields1{j})(i,1)=struct1.(fields1{j})(i);
        end

    S.source_struct(i,1)=1;
    end
end
counter=0;
    for z=slength(struct1)+1:((slength(struct1)+slength(struct2)))
        counter=counter+1;
        i=counter;
        for j=1:length(fields1)
             if iscell(struct2.(fields1{j})(i))
                S.(fields1{j}){z,1}=struct2.(fields1{j}){i};
            else
                S.(fields1{j})(z,1)=struct2.(fields1{j})(i);
            end
                S.source_struct(z,1)=2;

        end
    end


end


function S=concatStructs_differentFields(struct1,struct2,CommonFields,nOV1,nOV2)
for i=1:slength(struct1)
    for j=1:length(CommonFields)
       
            
        if iscell(struct1.(CommonFields{j})(i))
        S.(CommonFields{j}){i,1}=struct1.(CommonFields{j}){i};
        else
                S.(CommonFields{j})(i,1)=struct1.(CommonFields{j})(i);
        end

    end
            S.source_struct(i,1)=1;

end

counter=0;
    for z=slength(struct1)+1:((slength(struct1)+slength(struct2)))
        counter=counter+1;
        i=counter;
        for j=1:length(CommonFields)
             if iscell(struct2.(CommonFields{j})(i))
                S.(CommonFields{j}){z,1}=struct2.(CommonFields{j}){i};
            else
                S.(CommonFields{j})(z,1)=struct2.(CommonFields{j})(i);
            end

        end
                        S.source_struct(z,1)=2;

    end

c1=1;
c2=1;
for k=1:slength(S)
    if S.source_struct(k)==1;
        for j=1:length(nOV2)
            if iscell(struct2.(nOV2{j})(1))
                S.(nOV2{j}){k,1}='NA';
            else
                S.(nOV2{j})(k,1)=nan;
            end
        end
        for j=1:length(nOV1)
                if iscell(struct1.(nOV1{j})(1))        
                    S.(nOV1{j}){k,1}=struct1.(nOV1{j}){c1};
                    c1=c1+1;
                else
                    S.(nOV1{j})(k,1)=struct1.(nOV1{j})(c1);
                    c1=c1+1;
                end
    
        end
    end  
    if S.source_struct(k)==2;
        for j=1:length(nOV1)
            if iscell(struct1.(nOV1{j})(1))
                S.(nOV1{j}){k,1}='NA';
            else
                S.(nOV1{j})(k,1)=nan;
            end
        end
        for j=1:length(nOV2)
                if iscell(struct2.(nOV2{j})(1))        
                    S.(nOV2{j}){k,1}=struct2.(nOV2{j}){c2};
                    c2=c2+1;
                else
                    S.(nOV2{j})(k,1)=struct2.(nOV2{j})(c2);
                    c2=c2+1;
                end
    
        end
    end
    
    end
end