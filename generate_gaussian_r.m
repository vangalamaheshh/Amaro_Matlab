function reward=generate_gaussian_r(tree_top,tree_bot,tree_dist,m_pos,vel)%monkey_top,monkey_bot,vel)
target_y=tree_top-(tree_top-tree_bot)/2;
%m_pos=monkey_top-(monkey_top-monkey_bot)/2;
t=tree_dist/vel;
reward=normpdf(m_pos,target_y,40*t+.001);

end

function example
tree_top=400;
tree_bot=200;
tree_dist=5;


end