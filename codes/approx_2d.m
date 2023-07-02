function y = approx_2d(X, Y, x)


Nx = length(X);
% Mx = round(Nx/2);
j = 2;
y = nan(length(x),1);
for i = 1: length(x)
	while (j < Nx && x(i) > X(j))
		j = j + 1;
	end
	
	% First order approximation
	dy = (Y(j)-Y(j-1))/(X(j)-X(j-1));
	dx = x(i)-X(j-1);
	y(i) = Y(j-1) + dy*dx;

    % approximation outside bounds: linear extrapolation (OLS)
    Nextrap = 10;
    Nextrap = min(Nextrap,Nx);
if (x(i) < X(1))
    X_reg = [ones(Nextrap,1) reshape(X(1:Nextrap),Nextrap,1)];
    Y_reg = reshape(Y(1:Nextrap),Nextrap,1);
    beta  = (X_reg'*X_reg)\X_reg'*Y_reg;
    y(i)  = [1 x(i)] * beta;
%     j = 2;
% 	dy = (Y(j)-Y(j-1))/(X(j)-X(j-1));
% 	dx = x(i)-X(j-1);
% 	y(i) = Y(j-1) + dy*dx;
elseif (x(i) > X(Nx))
    X_reg = [ones(Nextrap,1) reshape(X(Nx-Nextrap+1:end),Nextrap,1)];
    Y_reg = reshape(Y(Nx-Nextrap+1:end),Nextrap,1);
    beta  = (X_reg'*X_reg)\X_reg'*Y_reg;
    y(i)  = [1 x(i)] * beta;

%     j = Nx;
% 	dy = (Y(j)-Y(j-1))/(X(j)-X(j-1));
% 	dx = x(i)-X(j);
%     y(i) = Y(j) + dy*dx;
end

end

% y_fn = griddedInterpolant(X,Y);
% 
% y = y_fn(x);
% y =y';

% % Second order approximation outside bounds
% if (x(i) < X(1))
%     d2y = (Y(3)-Y(2))/(X(3)-X(2)) - (Y(2)-Y(1))/(X(2)-X(1));
%     d2y = d2y*2/(X(3)-X(1));
%     y(i) = y(i) + 0.5*d2y*dx^2;
% elseif (x(i) > X(Nx))
%     d2y = (Y(Nx)-Y(Nx-1))/(X(Nx)-X(Nx-1)) - (Y(Nx-1)-Y(Nx-2))/(X(Nx-1)-X(Nx-2));
%     d2y = d2y*2/(X(Nx)-X(Nx-2));
%     y(i) = y(i) + 0.5*d2y*dx^2;
% end