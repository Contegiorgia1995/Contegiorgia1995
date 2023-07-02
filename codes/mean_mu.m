function out = mean_mu(mu,par)

[H0, par.PG]         = discretize_ar_rouwen(mu,par.beta_inc, par.sigmah0, par.N_fe);
par.H0               = H0;


aux = par.PG^1000;
out = aux(1,:)*exp(H0);
end
