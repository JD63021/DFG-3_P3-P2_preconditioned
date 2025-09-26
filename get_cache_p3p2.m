function cache = get_cache_p3p2(nodeInfo, elemInfo)
%GET_CACHE_P3P2  Geometry/shape cache for P3–P2 (triangles, vel=P3 with 10 DOFs/elem).
% Returns a struct with:
%   NV   : (10 x nGp)   P3 velocity shape values at triangle Gauss points (reference → pulled to physical via detJ only)
%   dNx{e}, dNy{e} : (10 x nGp) physical-gradients of P3 shapes for element e
%   wdet{e}        : (1  x nGp) Gauss weights * det(J) per element
%
% Rebuilds only when mesh size (#vel nodes / #elements) changes.
%
% Requires on path:
%   precomputeShapeFunctionsP3_Tri() -> [NV, dNxiV, dNetaV, wV]
%   p3ShapeDerivativesAllNodes(xe, ye, dNxi_col, dNeta_col) -> [dNx_col, dNy_col, detJ]

persistent C Nxy_s nEl_s

Nxy  = length(nodeInfo.velocity.x);
nEl  = size(elemInfo.velElements, 1);

% Rebuild if cache is empty or mesh size changed
if isempty(C) || ~isequal(Nxy, Nxy_s) || ~isequal(nEl, nEl_s)
    % P3 reference shapes/derivs/weights on the reference triangle
    [NV, dNxiV, dNetaV, wV] = precomputeShapeFunctionsP3_Tri();   % NV: (10 x nGp)
    nGp = numel(wV);

    C          = struct;
    C.NV       = NV;                 % same for all elements (reference evals)
    C.dNx      = cell(nEl,1);
    C.dNy      = cell(nEl,1);
    C.wdet     = cell(nEl,1);

    % Build per-element physical gradients and detJ*weights
    for e = 1:nEl
        Kvel = elemInfo.velElements(e,:);                   % 10 velocity nodes (P3)
        xV = nodeInfo.velocity.x(Kvel);  yV = nodeInfo.velocity.y(Kvel);

        dNx = zeros(10, nGp);
        dNy = zeros(10, nGp);
        wdet = zeros(1,  nGp);

        for gp = 1:nGp
            % reference derivatives at this GP
            dNxi  = dNxiV(:, gp);
            dNeta = dNetaV(:, gp);
            % map to physical space
            [dNx(:,gp), dNy(:,gp), detJ] = p3ShapeDerivativesAllNodes(xV, yV, dNxi, dNeta);
            wdet(gp) = wV(gp) * detJ;
        end

        C.dNx{e}  = dNx;
        C.dNy{e}  = dNy;
        C.wdet{e} = wdet;
    end

    % remember sizes to detect changes
    Nxy_s = Nxy;
    nEl_s = nEl;
end

cache = C;
end
