M=load_struct('/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/Stilgenbauer.ICGC.DFCI.538.Noncoding.AT.maf'); 
m=reorder_struct(M,ismember(M.Dan_clonal,'1'));
 save_struct(m,'/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/clonal.maf') 
 
 m=reorder_struct(M,ismember(M.Dan_clonal,'0'));
 save_struct(m,'/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/subclonal.maf')
 
 m=reorder_struct(M,ismember(M.IGHV,'0'));      
 save_struct(m,'/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/IGHV.negative.maf')
 
 m=reorder_struct(M,~ismember(M.IGHV,'0'));
 save_struct(m,'/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/IGHV.positive.maf')
 
 m=reorder_struct(M,~ismember(M.IGHV,'0')&ismember(M.Dan_clonal,'1'));
 save_struct(m,'/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/clonal_IGHVneg.maf')   

 m=reorder_struct(M,ismember(M.IGHV,'0')&ismember(M.Dan_clonal,'0'));
 save_struct(m,'/xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/subclonal_IGHVneg.maf')
 
 m=reorder_struct(M,~ismember(M.IGHV,'0')&ismember(M.Dan_clonal,'0'));
 save_struct(m,' /xchip/cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/subclonal_IGHVpos.maf')