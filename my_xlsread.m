function s=my_xlsread(fname)


javaaddpath /home/radon00/gadgetz/java/jexcel/jexcelapi/jxl.jar
import java.io.File;
import jxl.*;

s={};
f=File(fname);
if f.exists()
  workbook=Workbook.getWorkbook(f);
  sh=workbook.getSheet(0);
  r=sh.getRows;
  c=sh.getColumns;
  disp(['r=' num2str(r) '; c=' num2str(c)]);
  for ri=1:r
    for ci=1:c
      cl=sh.getCell(ci-1,ri-1);
      s{ri,ci}=cl.getContents.toCharArray;
    end
  end
end
