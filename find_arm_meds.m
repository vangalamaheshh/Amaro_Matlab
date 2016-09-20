function [armsize, Xmeds, Ymeds, DAT, MEDS] =find_arm_meds(Y, cyto)
    pos = [Y.pos];
    DAT = [Y.dat];
    names = {cyto.name};
    chrn = [Y.chrn];
    cytostarts=[cyto.start];
    
    a=size(DAT);
    a=a(2);
    pmed=zeros(1, a);
    qmed=zeros(1, a);
    ind = regexprep(names, '.+p.+', '0');
    ind = regexprep(ind, '.+q.+', '1');
    ind = cell2mat(ind);
    ind2 = [0 ind(1:end-1)];
    ind3 = abs(ind2 - ind);
    ind3(1)=0;
    start_ind = find(ind3);
    starts = [0 cytostarts(start_ind)];
    MEDS = zeros(1,a);
    armsize=zeros(1, 2);
    starts = starts([2 24 32 34 36 38 40 42 44 4 6 8 10 12 14 16 18 20 22 26 28 30 46 48]);
    for i = 1:24
                disp(i)
                D_chr = DAT((chrn==i),:);
                D_pos=pos((chrn==i));
                p=D_chr(D_pos < starts(i),:);
                q=D_chr(D_pos >= starts(i),:);
                pmed = nanmedian(p);
                qmed = nanmedian(q);
                if(i==23)
                    Xmeds=[pmed qmed];
                end
                if (i == 24)
                    Ymeds=[pmed qmed];
                end
                
                [psize n] = size(p);
                [qsize n] = size(q);
                pmed = ones(psize, 1)*pmed;
                qmed = ones(qsize, 1)*qmed;
                MEDS = [MEDS;pmed ;qmed];
                armsize=[armsize;psize qsize];
    end
    armsize=armsize(2:end,:);
    MEDS=MEDS(2:end,:);
    DAT=DAT-MEDS;
     
end
                