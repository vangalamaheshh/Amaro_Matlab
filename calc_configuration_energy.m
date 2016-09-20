function e=calc_configuration_energy(S,J,q,has_labels)

if has_labels
    e=sum(J( find(S(J(has_labels:end,1)) ~= S(J(has_labels:end,2)))+has_labels-1 ,3));    
else
    e=sum(J( S(J(:,1)) ~= S(J(:,2)) ,3));
end