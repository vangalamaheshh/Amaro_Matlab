%Quck and Dirty Figure for Cathy

patients=unique(MTp1.individual);

MTp1.key=strcat(MTp1.Chromosome,MTp1.Start_position);
MTp2.key=strcat(MTp2.Chromosome,MTp2.Start_position);
MTp1.i_tumor_f=str2double(MTp1.i_tumor_f);
MTp2.i_tumor_f=str2double(MTp2.i_tumor_f);

for i=1:length(patients)
    matches=unique(MTp2.Tumor_Sample_Barcode(ismember(MTp2.individual,patients{i})));
    xmaf=reorder_struct(MTp1,ismember(MTp1.individual,patients{i}));
    for j=1:length(matches)
    ymaf=reorder_struct(MTp2,ismember(MTp2.Tumor_Sample_Barcode,matches{j}));
    figure()
    ylabel(matches{j});
    xlabel(strcat(patients{i},'-TP1'));
    hold on
    for m=1:slength(xmaf)
    
        x=xmaf.i_tumor_f(m);%allele fraction TP1
        k=find(ismember(ymaf.key,xmaf.key{m}));
        if isempty(k)
            y=0;
        else
            y=ymaf.i_tumor_f(k);%allele fraction TP2;
        end
       if ismember(xmaf.Hugo_Symbol{m}, CLLsig)
       sigGene=find(ismember(CLLsig,xmaf.Hugo_Symbol{m}));
       text(x,y,CLLsig{sigGene});
       plot(x,y,'r.');
       else
       plot(x,y,'k.') ;
       end
    
    
    end
    saveas(gcf,sprintf('/Users/amaro/Documents/Quick2Dplotsforcathy/%s-%s.jpg',patients{i},matches{j}))
    close(gcf);
    
    end
    

end