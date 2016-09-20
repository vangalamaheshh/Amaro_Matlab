function a = add_memfield(a,fieldname,data)

a.fieldnames = [a.fieldnames fieldname];
a.fielddata = [a.fielddata data];
a.storagetype = [a.storagetype 'memory'];
a.colmapping = [a.colmapping {[]} ];
a.rowmapping = [a.rowmapping {[]} ];