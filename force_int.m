function [F_x, F_y] = force_int(u, v, p, ...
    nodeInfo, elemInfo, boundaryInfo, flag, mu, cRef)
% FORCE_INT_P3P2 - Boundary traction integration for a P3–P2 triangular approach.
%
% This function integrates the stress vector on a specified boundary (flag)
% to compute drag (F_x) and lift (F_y).
%
% Inputs:
%   u, v, p    -> Current velocity (u,v) and pressure (p) solutions
%                 in global DOF ordering:
%                   - u, v each size (#velNodes), for P3
%                   - p is size (#presNodes), for P2
%   nodeInfo   -> .velocity.x, .velocity.y for the P3 mesh
%                 .pressure.x, .pressure.y for the P2 mesh
%   elemInfo   -> .velElements (#elems×10) for P3 triangles
%                 .presElements (#elems×6)  for P2 triangles
%   boundaryInfo -> .velLine4Elements.flag_<flag> => Nx4 lines
%                   (each row has [cornerA, cornerB, mid1, mid2], possibly,
%                   or at least columns 1,2 are corner nodes).
%   flag       -> integer specifying which boundary to integrate
%   mu         -> dynamic viscosity
%   cRef       -> (optional) [cx, cy], a reference point to ensure outward normal
%                 flips correctly; default is [0,0] if not provided
%
% Outputs:
%   F_x, F_y   -> total force (x- and y-components) on that boundary.

if nargin < 9
    cRef = [0.2, 0.2];  % default for outward normal check
end
cx = cRef(1);
cy = cRef(2);

F_x = 0;
F_y = 0;

%% 1) Access the boundary line data for the chosen 'flag'
lineField = sprintf('flag_%d', flag);
if ~isfield(boundaryInfo, 'velLine4Elements') || ...
   ~isfield(boundaryInfo.velLine4Elements, lineField)
    fprintf('No P3 boundary lines found for flag %d.\n', flag);
    return;
end
lineEls = boundaryInfo.velLine4Elements.(lineField);
if isempty(lineEls)
    fprintf('No line elements stored for flag %d.\n', flag);
    return;
end

%% 2) Choose a 1D Gauss rule for [0..1] along each boundary segment
%    (Here a 4-point rule)
sGP = [0.0694318, 0.3300095, 0.6699905, 0.9305682];
wGP = [0.1739274, 0.3260726, 0.3260726, 0.1739274];

%% 3) Loop over each boundary segment
%    lineEls(i,:) might be [cornerA, cornerB, mid1, mid2],
%    but we only need cornerA= lineEls(i,1), cornerB= lineEls(i,2].
for iLn = 1:size(lineEls,1)
    rawNodes = lineEls(iLn,:);
    cornerA  = rawNodes(1);
    cornerB  = rawNodes(2);

    % 3.1) Find which P3 triangle + which local side has these 2 corner nodes
    [elID, locSide] = findP3TriangleSide(elemInfo.velElements, [cornerA, cornerB]);
    if elID < 1
        % Not found => skip
        continue;
    end

    % 3.2) Gather that triangle's 10 velocity DOFs
    Kvel = elemInfo.velElements(elID,:);  % 10 node IDs for P3
    xcoord = nodeInfo.velocity.x(Kvel);
    ycoord = nodeInfo.velocity.y(Kvel);

    % 3.3) Gather that triangle's 6 pressure DOFs
    Kprs  = elemInfo.presElements(elID,:); % 6 node IDs for P2

    % 3.4) Physical coords of corners
    xA = nodeInfo.velocity.x(cornerA);
    yA = nodeInfo.velocity.y(cornerA);
    xB = nodeInfo.velocity.x(cornerB);
    yB = nodeInfo.velocity.y(cornerB);

    Fx_el = 0;
    Fy_el = 0;

    % 4) Integrate in 1D param s \in [0..1]
    for ig = 1:length(sGP)
        s = sGP(ig);
        w = wGP(ig);

        % 4.1) Physical coordinates on the segment
        xLoc = xA + s*(xB - xA);
        yLoc = yA + s*(yB - yA);

        % tangent vector
        tx = (xB - xA);
        ty = (yB - yA);
        ds_1d = sqrt(tx^2 + ty^2);

        % 4.2) Convert (xLoc,yLoc) to local (xi,eta) for the reference triangle
        %      using a simple side param approach
        [xi, eta] = triLocalParamMapping(xcoord, ycoord, xLoc, yLoc, locSide, s);

        % 4.3) Evaluate P3 shape & derivatives at (xi,eta) => velocity gradients
        [Nv, dNxi_v, dNeta_v] = p3basisGmsh(xi, eta);
        % Build the Jacobian for velocity tri
        dX_dxi  = sum(xcoord .* dNxi_v);
        dX_deta = sum(xcoord .* dNeta_v);
        dY_dxi  = sum(ycoord .* dNxi_v);
        dY_deta = sum(ycoord .* dNeta_v);

        Jmat = [dX_dxi, dY_dxi; dX_deta, dY_deta];
        detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
        if abs(detJ) < 1e-14
            continue;  % degenerate
        end
        invJ = inv(Jmat);

        % physical derivatives of velocity shape
        dNdx_v = invJ(1,1)*dNxi_v + invJ(1,2)*dNeta_v;
        dNdy_v = invJ(2,1)*dNxi_v + invJ(2,2)*dNeta_v;

        % velocity at this point
        uLocal = u(Kvel);
        vLocal = v(Kvel);
        dudx = sum(uLocal .* dNdx_v);
        dudy = sum(uLocal .* dNdy_v);
        dvdx = sum(vLocal .* dNdx_v);
        dvdy = sum(vLocal .* dNdy_v);

        % 4.4) Evaluate pressure (P2) at (xi,eta)
        [Np2] = p2basisGmsh(xi, eta);
        pLocal = p(Kprs);
        p_g = sum(pLocal .* Np2);

        % 4.5) Stress tensor => sigma = -p I + mu (gradU + gradU^T)
        T_xx = -p_g + 2*mu*dudx;
        T_xy = mu*(dudy + dvdx);
        T_yy = -p_g + 2*mu*dvdy;

        % 4.6) Outward normal
        nx =  ty/ds_1d;
        ny = -tx/ds_1d;
        rx = xLoc - cx;
        ry = yLoc - cy;
        if (rx*nx + ry*ny) < 0
            nx = -nx;
            ny = -ny;
        end

        % traction = sigma·n
        tx_g = T_xx*nx + T_xy*ny;
        ty_g = T_xy*nx + T_yy*ny;

        Fx_el = Fx_el + tx_g*(ds_1d*w);
        Fy_el = Fy_el + ty_g*(ds_1d*w);
    end

    % 5) Accumulate total
    F_x = F_x + Fx_el;
    F_y = F_y + Fy_el;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [elID, locSide] = findP3TriangleSide(allTri10, cornerSet)
% findP3TriangleSide:
%   Finds which P3 triangle (10 nodes) has corners cornerSet (2 node IDs)
%   among the first 3 corners of that element. 
%   Output:
%     elID    -> index of that triangle in allTri10
%     locSide -> 'side12','side23','side31' depending on which pair it is.

elID   = 0;
locSide= '';

for e = 1:size(allTri10,1)
    row = allTri10(e,:);
    % corners are row(1), row(2), row(3)
    c1 = row(1);
    c2 = row(2);
    c3 = row(3);

    if all(ismember(cornerSet, [c1,c2]))
        elID = e; locSide='side12'; return;
    elseif all(ismember(cornerSet, [c2,c3]))
        elID = e; locSide='side23'; return;
    elseif all(ismember(cornerSet, [c3,c1]))
        elID = e; locSide='side31'; return;
    end
end
end

function [xi, eta] = triLocalParamMapping(xcoord, ycoord, xLoc, yLoc, locSide, s)
% triLocalParamMapping:
%   For side12 => corner1= (xi,eta)=(0,0), corner2=(1,0)
%   side23 => corner2= (1,0), corner3=(0,1)
%   side31 => corner3= (0,1), corner1=(0,0)
%   We param from corner1->corner2, etc. 
%   Here we do a simple approach: side12 => (xi,eta)=(s,0), etc.

switch locSide
    case 'side12'
        xi  = s; eta = 0;
    case 'side23'
        xi  = 1 - s; eta = s;
    case 'side31'
        xi  = 0; eta = 1 - s;
    otherwise
        % fallback
        xi= s; eta=0;
end
end

function [N] = p2basisGmsh(xi, eta)
% p2basisGmsh:
%   Returns the 6 shape functions for a quadratic (P2) triangle,
%   anticlockwise node ordering:
%    (1) corner1, (2) corner2, (3) corner3,
%    (4) edge(1-2), (5) edge(2-3), (6) edge(3-1).
%
%   We do not need derivatives here, only N for pressure.
%
zeta = 1 - xi - eta;
N1 = zeta*(2*zeta - 1);
N2 = xi*(2*xi - 1);
N3 = eta*(2*eta - 1);
N4 = 4*xi*zeta;
N5 = 4*xi*eta;
N6 = 4*eta*zeta;

N = [N1; N2; N3; N4; N5; N6];
end

function [N, dNxi, dNeta] = p3basisGmsh(xi, eta)
% p3basisGmsh:
%   Returns the 10 shape functions + derivatives for a cubic (P3) triangle,
%   anticlockwise node ordering:
%     node1-> corner1, node2->corner2, node3->corner3,
%     node4,node5-> mid-edge(1->2),
%     node6,node7-> mid-edge(2->3),
%     node8,node9-> mid-edge(3->1),
%     node10-> bubble.
%
%   If you have a precompute routine, you can call it. Here we do inline:

L1 = 1 - xi - eta;
L2 = xi;
L3 = eta;

N = zeros(10,1);
% corners
N(1) = L1*(3*L1 - 1)*(3*L1 - 2)/2;
N(2) = L2*(3*L2 - 1)*(3*L2 - 2)/2;
N(3) = L3*(3*L3 - 1)*(3*L3 - 2)/2;

% mid-edges
c = 9/2;
N(4) = c * L1*L2*(3*L1 - 1);
N(5) = c * L1*L2*(3*L2 - 1);
N(6) = c * L2*L3*(3*L2 - 1);
N(7) = c * L2*L3*(3*L3 - 1);
N(8) = c * L3*L1*(3*L3 - 1);
N(9) = c * L3*L1*(3*L1 - 1);

% bubble
N(10) = 27*L1*L2*L3;

% derivatives
dNxi = zeros(10,1);
dNeta= zeros(10,1);

% partials of L1=1-xi-eta => dL1dxi=-1, dL1deta=-1
% L2=xi => dL2dxi=1, dL2deta=0
% L3=eta=> dL3dxi=0, dL3deta=1
dL1_dxi=-1; dL1_deta=-1;
dL2_dxi= 1; dL2_deta= 0;
dL3_dxi= 0; dL3_deta= 1;

% For corners, we can do it similarly to your prior code.
% For mid-edges & bubble we can do product rule. Below is one example:

% corner1
f = (3*L1 -1)*(3*L1 -2)/2;
df_dL1 = (9*(2*L1 -1))/2;
dNxi(1)  = dL1_dxi*f + L1*(df_dL1*dL1_dxi);
dNeta(1) = dL1_deta*f + L1*(df_dL1*dL1_deta);

% corner2
f2 = (3*L2 -1)*(3*L2 -2)/2;
df2dL2 = (9*(2*L2 -1))/2;
dNxi(2)  = dL2_dxi*f2 + L2*(df2dL2*dL2_dxi);
dNeta(2) = dL2_deta*f2 + L2*(df2dL2*dL2_deta);

% corner3
f3 = (3*L3 -1)*(3*L3 -2)/2;
df3dL3 = (9*(2*L3 -1))/2;
dNxi(3)  = dL3_dxi*f3 + L3*(df3dL3*dL3_dxi);
dNeta(3) = dL3_deta*f3 + L3*(df3dL3*dL3_deta);

% mid-edge4: c * L1*L2*(3L1 -1)
[dNxi(4), dNeta(4)] = d_of_prod3(c, L1, L2, (3*L1 -1), ...
                                 dL1_dxi, dL1_deta, ...
                                 dL2_dxi, dL2_deta, ...
                                 3*dL1_dxi, 3*dL1_deta);

% mid-edge5: c * L1*L2*(3L2 -1)
[dNxi(5), dNeta(5)] = d_of_prod3(c, L1, L2, (3*L2 -1), ...
                                 dL1_dxi, dL1_deta, ...
                                 dL2_dxi, dL2_deta, ...
                                 3*dL2_dxi, 3*dL2_deta);

% mid-edge6: c * L2*L3*(3L2 -1)
[dNxi(6), dNeta(6)] = d_of_prod3(c, L2, L3, (3*L2 -1), ...
                                 dL2_dxi, dL2_deta, ...
                                 dL3_dxi, dL3_deta, ...
                                 3*dL2_dxi, 3*dL2_deta);

% mid-edge7: c * L2*L3*(3*L3 -1)
[dNxi(7), dNeta(7)] = d_of_prod3(c, L2, L3, (3*L3 -1), ...
                                 dL2_dxi, dL2_deta, ...
                                 dL3_dxi, dL3_deta, ...
                                 3*dL3_dxi, 3*dL3_deta);

% mid-edge8: c * L3*L1*(3*L3 -1)
[dNxi(8), dNeta(8)] = d_of_prod3(c, L3, L1, (3*L3 -1), ...
                                 dL3_dxi, dL3_deta, ...
                                 dL1_dxi, dL1_deta, ...
                                 3*dL3_dxi, 3*dL3_deta);

% mid-edge9: c * L3*L1*(3*L1 -1)
[dNxi(9), dNeta(9)] = d_of_prod3(c, L3, L1, (3*L1 -1), ...
                                 dL3_dxi, dL3_deta, ...
                                 dL1_dxi, dL1_deta, ...
                                 3*dL1_dxi, 3*dL1_deta);

% bubble10: 27*L1*L2*L3
[dNxi(10), dNeta(10)] = d_of_prod3(27, L1, L2, L3, ...
                                   dL1_dxi, dL1_deta, ...
                                   dL2_dxi, dL2_deta, ...
                                   dL3_dxi, dL3_deta);
end

function [dNdxi, dNeta] = d_of_prod3(scale, X, Y, Z, ...
                                     dXdxi, dXdeta, ...
                                     dYdxi, dYdeta, ...
                                     dZdxi, dZdeta)
% Product rule for N=scale*X*Y*Z
dNdxi = scale*( dXdxi*Y*Z + X*dYdxi*Z + X*Y*dZdxi );
dNeta = scale*( dXdeta*Y*Z + X*dYdeta*Z + X*Y*dZdeta );
end

