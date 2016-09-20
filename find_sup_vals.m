function v=find_sup_vals(D,supid,sup_vals)

[typeacc,typedesc,D,range,non_empty]=decollapse_supdat(D,supid,1:max(D.supdat(supid,:)),0);
D=reorder_D_sup(D,'cols',range);
v=find_supid(D,sup_vals);

