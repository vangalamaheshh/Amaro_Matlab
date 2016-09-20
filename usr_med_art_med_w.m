d = dataset('file', '~/Downloads/train-2.csv', 'ReadObsNames', false, 'ReadVarNames', true, 'Delimiter', ',');
art = dataset('file','~/Downloads/artists.csv','ReadObsNames',false,'ReadVarNames',true,'Delimiter',',');
load('~/Downloads/m-2.mat')

art_m_order=unique(d.artist);
users=unique(d.user);


CV=10;
artist_med_cv=zeros(2000,10);
usr_med_cv=zeros(round(233286/CV),10);

for cvs=1:CV
    disp(sprintf('Generating Fold %d',cvs))
    cv_index=randi(length(m),[round(233286/CV),1]);
    cv_user(:,cvs)={users{cv_index}};
    artist_med_cv(:,cvs)=nanmedian(m(cv_index,:))';
    usr_med_cv(:,cvs)=nanmedian(m(cv_index,:),2);
    
end

[x y]=ismember(d.artist,art_m_order);
for cvs=1:CV
    disp(sprintf('Running Fold %d',cvs))
    [i mm_cv]=ismember(d.user,cv_user(:,cvs)); 
    [x y_cv]=ismember(d.artist(i),art_m_order);
    for a_p=1:length(alpha)
        predict_played=(usr_med_cv(mm_cv(mm_cv>0),cvs)*alpha(a_p))+(artist_med_cv(y_cv,cvs)*(1-alpha(a_p)));
        alpha_loss(cvs,a_p)=mean(abs(d.plays(i)-predict_played(:,a_p)));
        if mod(a_p,.2)==0
            a_p
        end
    end
end
   


artist_med=nanmedian(m);
artist_mean=nanmean(m);
usr_med=nanmedian(m,2);

alpha=[.01:.01:.99];
stat_pro=sort_struct(stat_pro,'user_id',1);
[i mm]=ismember(d.user,stat_pro.user_id);
[x y]=ismember(d.artist,art_m_order);
predict_played=zeros(length(d.user),length(alpha));
for i=1:length(alpha)
    predict_played(:,i)=(usr_med(mm)*alpha(i))+(artist_med(y)'*(1-alpha(i)));
    alpha_loss(i)=mean(abs(d.plays-predict_played(:,i)));
    i
end


test = dataset('file', '~/Downloads/test-2.csv', 'ReadObsNames', false, 'ReadVarNames', true, 'Delimiter', ',');

[i mm]=ismember(test.user,stat_pro.user_id);
[x y]=ismember(test.artist,art_m_order);
predict_played=zeros(length(test.user),1);
predict_played=(usr_med(mm)*.7)+(artist_med(y)'*(1-.3));









int_d=string_to_categ_var(d.artist);
