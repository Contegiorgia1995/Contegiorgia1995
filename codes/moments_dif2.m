function [ out ] = moments_dif2( x,data )
    rho     = exp(x(1))/(1+exp(x(1)));
    var_z0  = exp(x(2));
    var_m   = exp(x(3));
    var_v   = exp(x(4));
    
    age_group  = data(:,1);
    lag        = data(:,2);
    cov        = data(:,3);
    
    out        = zeros(length(cov),1);
    
    for j = 1:length(cov)
        if age_group(j) == 0
            if lag(j) == 0
                out(j) = (var_z0 + var_m - cov(j));
            elseif lag(j) >= 1
                out(j) = (rho^lag(j).*var_z0 - cov(j));
            end
        elseif age_group(j) > 0
            if lag(j) == 0
                out(j) = (rho.^(2*age_group(j))*var_z0 + var_m + var_v * sum(rho.^(2.*(0:age_group(j)-1))) - cov(j));
            elseif lag(j) >= 1
                out(j) = (rho.^(2*age_group(j)+lag(j))*var_z0 + rho^lag(j).* var_v * sum(rho.^(2.*(0:age_group(j)-1))) - cov(j));
            end
        end    
    end
    out = sum(out.^2)./length(out);
end

