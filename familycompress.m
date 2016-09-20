for i=1:slength(mafc)                                                                                
if length(find(ismember(mafc.uniq,mafc.uniq{i}))==1)>1 && isempty(regexp(mafc.patient{i},'-TP-N','match'))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
 
    reject(i)=1;
    
end
end







    