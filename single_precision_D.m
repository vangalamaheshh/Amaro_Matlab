function D=single_precision_D(D)

flds={'dat','affy_call','affy_calls','adat','supdat','gsupdat','raw','sm1','sm2','sm2j','sm3','smooth','cbs','cbs_fixed'};

for i=1:length(flds)
  if isfield(D,flds{i})
    D=setfield(D,flds{i},single(getfield(D,flds{i})));
  end
end
