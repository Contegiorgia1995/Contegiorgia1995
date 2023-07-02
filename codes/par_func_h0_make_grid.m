function eq = par_func_h0_make_grid(gridh0,desired_points,desired_probs)

gridh0  = sort(gridh0);
n_data  = length(desired_points);
eqaux   = nan(1,n_data);

for ii = 1:n_data
    aux_dif     = (gridh0 - desired_points(ii)).^2;
    eqaux(ii)	= min(aux_dif);
end

eq              = eqaux * desired_probs;

end

