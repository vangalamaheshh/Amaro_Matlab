function [beta,variance]=beta_pdf_on_call_stats(alt,ref)
af=[0:.01:1];
beta=betapdf(af,alt+1,ref+1); 
variance=((ref+1)*(alt+1))/((alt+ref+2)^2*(ref+alt+3));

end