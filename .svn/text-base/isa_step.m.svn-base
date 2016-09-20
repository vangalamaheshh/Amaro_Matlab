function [g1,c1]=isa_step(g0,t_g,t_c,EG,ECT,down)
verbose=0;

sr=[]; sc=[];

c_proj=EG*g0;
m_c_proj=mean(c_proj);
s_c_proj=std(c_proj);


if down
  c1=(c_proj-m_c_proj+eps)./(s_c_proj+eps);
  sc1=sign(c1);
  c1=abs(c1)-t_c;
  c1(c1<0)=0;
  c1=sparse(c1);
  c1=c1.*sc1;
else
  c1=(c_proj-m_c_proj+eps)./(s_c_proj+eps)-t_c;
  c1(c1<0)=0;
  c1=sparse(c1);
end

g_proj=ECT*c1;
m_g_proj=mean(g_proj);
s_g_proj=std(g_proj);

if down
  g1=(g_proj-m_g_proj+eps)./(s_g_proj+eps);
  sg1=sign(g1);
  g1=abs(g1)-t_g;
  g1(g1<0)=0;
  g1=sparse(g1);
  g1=g1.*sg1;
else
  g1=(g_proj-m_g_proj+eps)./(s_g_proj+eps)-t_g;
  g1(g1<0)=0;
  g1=sparse(g1);
end

if(verbose)
  figure(1)
  subplot(2,1,1);
  hist(c_proj,100);
  ax=axis;
  line([ m_c_proj m_c_proj ],[ax(3:4)],'Color','red');
  line([ m_c_proj+t_c*s_c_proj m_c_proj+t_c*s_c_proj ],[ax(3:4)],'Color','blue');
  
  subplot(2,1,2);
  hist(g_proj,100);
  ax=axis;
  line([ m_g_proj m_g_proj ],[ax(3:4)],'Color','red');
  line([ m_g_proj+t_g*s_g_proj m_g_proj+t_g*s_g_proj ],[ax(3:4)],'Color','blue');
end




