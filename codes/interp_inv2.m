function Q = interp_inv2(VEC_or,val_or)
%#codegen
% VEC and val must be row vectors
% VEC must be sorted as increasing vector


%% Reshape val and VEC
if iscolumn(val_or) == 1
    val = val_or';
else
    val = val_or;
end

if iscolumn(VEC_or) == 1
    VEC = VEC_or';
else
    VEC = VEC_or;
end

N0 = length(val);
N1 = length(VEC);
val2 = repmat(val',1,N1);
VEC2 = repmat(VEC,N0,1);

%% Corners:
% Left corner: val<=VEC(1)
pos_LC = (val2<=VEC2);
pos_LC = sum(pos_LC,2);
pos_LC = (pos_LC == N1);
Q_LC   = zeros(N0,N1);
Q_LC(pos_LC,1)   = 1;

% Right corner: val>=VEC(end)
pos_RC = (val2>=VEC2);
pos_RC = sum(pos_RC,2);
pos_RC = (pos_RC == N1);
Q_RC   = zeros(N0,N1);
Q_RC(pos_RC,end)   = 1;

%% Interior values
% Right: VEC>= val , first
A = (VEC2>=val2);
[~,pos_R] = max(A,[],2);
% Left: R_pos - 1
pos_L = (pos_R>1) .* (pos_R-1) + (pos_R==1) .* 2;

% Weights
Vval_L = VEC(pos_L);
Vval_R = VEC(pos_R);


W_L    = (Vval_R-val)./(Vval_R-Vval_L);
W_L    = W_L .* (pos_LC' ~= 1) .* (pos_RC' ~= 1);

W_R    = (val-Vval_L)./(Vval_R-Vval_L);
W_R    = W_R .* (pos_LC' ~= 1) .* (pos_RC' ~= 1);


NC = (pos_LC' ~= 1) .* (pos_RC' ~= 1);

Q_L = zeros(N0,N1);
Q_R = zeros(N0,N1);

size_W   = size(Q_L);
ind_aux1 = cumprod([1 size_W(1:end-1)]);

pos      = [(1:N0)' pos_L];
index    =  ind_aux1 * (pos' - repmat([0; ones(numel(pos(1,:))-1, 1)],1,size(pos,1)));
Q_L(index) = W_L;

pos      = [(1:N0)' pos_R];
index    =  ind_aux1 * (pos' - repmat([0; ones(numel(pos(1,:))-1, 1)],1,size(pos,1)));
Q_R(index) = W_R;

%% Transition matrix:

Q = Q_L + Q_R + Q_LC + Q_RC;
end
