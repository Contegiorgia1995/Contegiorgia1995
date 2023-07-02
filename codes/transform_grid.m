function [grid2,MM] = transform_grid(grid_orig)
% Transform par.gridh and par.grids into a matrix to use in parfor
% grid2(j,1:MM(j)) = grid{j}

NN = length(grid_orig);
MM = zeros(NN,1);       % length of each par.gridh{j}

for j=1:NN
    MM(j) = length(grid_orig{j});
end

MMmax = max(MM);
grid2 = nan(NN,MMmax);

for j=1:NN
   grid2(j,1:MM(j)) = grid_orig{j}; 
end


end
