function t=s2n(dat,cl1,cl2)
    n_1=length(cl1);
    m_1=mean(dat(:,cl1),2);
    s_1=std(dat(:,cl1),0,2);

    n_2=length(cl2);
    m_2=mean(dat(:,cl2),2);
    s_2=std(dat(:,cl2),0,2);

    t=(m_1-m_2)./(s_1+s_2);
