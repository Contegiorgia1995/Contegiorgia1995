function [Vval,Vpos,W] = interp_inv(VEC,val)
% Inverse interpolation: VEC is a vector where we find the nearest right
% and left values such that Vp(1)<= val <Vp(2) and W are the weight for
% Vp
% VEC must be an increasing vector

Vpos = zeros(1,2);
Vval = zeros(1,2);
W    = zeros(1,2);


if val< VEC(1)
%     fprintf('val< VEC(1)')
    Vpos = [1 2];
    Vval = VEC(Vpos);
    W    = [1 0];
elseif VEC(end)<= val
%     fprintf('VEC(end) < val')
    Vpos = [length(VEC)-1 length(VEC)];
    Vval = VEC(Vpos);
    W    = [0 1];
else
    Vpos(1) = find(val>= VEC,1,'last');
    Vpos(2) = find(val< VEC,1,'first');
    Vval(1) = VEC(Vpos(1));
    Vval(2) = VEC(Vpos(2));
    
    W(1) = (Vval(2) - val)   / (Vval(2)-Vval(1));
    W(2) = (val   - Vval(1)) / (Vval(2)-Vval(1));    
end

end