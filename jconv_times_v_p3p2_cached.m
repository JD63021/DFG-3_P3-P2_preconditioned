function y = jconv_times_v_p3p2_cached(cache, Re, rho, U, v, boundaryInfo, corner)
% Matrix-free action of the convective Jacobian for P3–P2 triangles.
% Robust to cache without Npr: it infers Npr from length(v).

Nxy = cache.Nxy;
Npr = length(v) - 2*Nxy;                     % infer pressure size
if Npr < 0
    error('jconv_times_v_p3p2_cached: length(v) inconsistent with cache.Nxy');
end

y_u = zeros(Nxy,1);
y_v = zeros(Nxy,1);

U1 = U(1:Nxy);          U2 = U(Nxy+1:2*Nxy);
v1 = v(1:Nxy);          v2 = v(Nxy+1:2*Nxy);

for e = 1:cache.numEl
    K   = cache.el(e).Kvel;
    Ni  = cache.el(e).Ni;          % (10 x G)
    dNx = cache.el(e).dNx;         % (10 x G)
    dNy = cache.el(e).dNy;         % (10 x G)
    W   = cache.el(e).W(:);        % (G x 1)

    U1e = U1(K);  U2e = U2(K);
    v1e = v1(K);  v2e = v2(K);

    a1   = Ni.'  * U1e;   a2   = Ni.'  * U2e;
    du1x = dNx.' * U1e;   du1y = dNy.' * U1e;
    du2x = dNx.' * U2e;   du2y = dNy.' * U2e;

    b1   = Ni.'  * v1e;   b2   = Ni.'  * v2e;
    dv1x = dNx.' * v1e;   dv1y = dNy.' * v1e;
    dv2x = dNx.' * v2e;   dv2y = dNy.' * v2e;

    cx = a1.*dv1x + a2.*dv1y + b1.*du1x + b2.*du1y;
    cy = a1.*dv2x + a2.*dv2y + b1.*du2x + b2.*du2y;

    y_u(K) = y_u(K) + rho*Re * (Ni * (cx .* W));
    y_v(K) = y_v(K) + rho*Re * (Ni * (cy .* W));
end

% append pressure zeros
y = [y_u; y_v; zeros(Npr,1)];

% Dirichlet velocity rows → identity in J
b = [];
if isfield(boundaryInfo,'allVelNodes') && ~isempty(boundaryInfo.allVelNodes)
    b = boundaryInfo.allVelNodes(:);
elseif isfield(boundaryInfo,'allNodes') && ~isempty(boundaryInfo.allNodes)
    b = boundaryInfo.allNodes(:);
end
if ~isempty(b)
    b = unique(b); b = b(b>=1 & b<=Nxy);
    y(b)       = v(b);
    y(Nxy + b) = v(Nxy + b);
end

% pinned pressure row
idxP = 2*Nxy + corner;
if idxP >= 1 && idxP <= length(v)
    y(idxP) = v(idxP);
end
end
