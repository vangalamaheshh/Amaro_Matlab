g = Genome16bit();

R = [];
R.chr = [1 2]';
R.start = [1000 5001]';
R.end =   [1007 5006]';
iname = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_index.txt';
save_struct_noheader(R,iname);

r1 = round(32767*rand(8,1));
r2 = round(32767*rand(6,1));

g.setContents(1,1000,1007,r1);
g.setContents(2,5001,5006,r2);

width=16
fwbname = ['/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_width' num2str(width) '.fwb'];
g.saveFixedWidthBinary(fwbname,width,0,iname);
a{width} = load_textfile(fwbname); binary(a{width})
aa=double(reshape(a{width},2,14))
[[r1;r2] (aa(1,:)*256+aa(2,:))']

width=24
fwbname = ['/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_width' num2str(width) '.fwb'];
g.saveFixedWidthBinary(fwbname,width,0,iname);
a{width} = load_textfile(fwbname); binary(a{width})
aa=double(reshape(a{width},3,14))
[[r1;r2] (aa(2,:)*256+aa(3,:))']



%%%%%%%%%%%%%%%%


g = Genome();

R = [];
R.chr = [1 2 2]';
R.start = [1000 5001 5011]';
R.end =   [1007 5006 5019]';
iname = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_index.txt';
save_struct_noheader(R,iname);

r1 = round(1.5*rand(8,1));
r2 = round(1.5*rand(6,1));
r3 = round(1.5*rand(9,1));

g.setContents(1,1000,1007,r1);
g.setContents(2,5001,5006,r2);
g.setContents(2,5011,5019,r3);

width=1;
fwbname = ['/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_width' num2str(width) '.fwb'];
iname = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_index.txt';
g.saveFixedWidthBinary(fwbname,width,0,iname);

width=1;
fwbname = ['/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_width' num2str(width) '.fwb'];
iname = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_index.txt';
fwb = FixedWidthBinary(fwbname,width,iname)
tic;fwb.get([1 2],[990 5001],[1007 5006]),toc







%%%%%%%%%%%%%%%%

width=1;
fwbname = '/xchip/cga1/lawrence/ov/analysis/20100818_3centers/ov316.fwb/OV-61-2102.somatic_coverage.fwb';
iname = '/xchip/cga1/lawrence/ov/analysis/20100818_3centers/ov316.fwb/ov316.fwi';
fwb = FixedWidthBinary(fwbname,width,iname,2,3,4)








width=24
fwbname = ['/xchip/cga2/lawrence/cga/trunk/matlab/seq/test_width' num2str(width) '.fwb'];
g.saveFixedWidthBinary(fwbname,width,0,iname);
a{width} = load_textfile(fwbname); binary(a{width})
aa=double(reshape(a{width},3,14))
[[r1;r2] (aa(2,:)*256+aa(3,:))']


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = '/xchip/cga1/lawrence/db/context65/all.fwb';
width = 8;
iname = '/xchip/cga1/lawrence/db/context65/all.fwi'
fwb = FixedWidthBinary(fname,width,iname);

fname = '/xchip/cga1/lawrence/ov/analysis/20100818_3centers/ov316.fwb/OV-61-2113.somatic_coverage.fwb';
width = 1;
iname = '/xchip/cga1/lawrence/ov/analysis/20100818_3centers/ov316.fwb/ov316.fwi';
fwb = FixedWidthBinary(fname,width,iname,2,3,4);

