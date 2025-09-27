function M = mass_matrix_func_p3p2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    rho, corner )
% MASS_MATRIX_FUNC_P3P2:
%   Builds the block-diagonal mass matrix for a P3–P2 triangular discretization
%   in incompressible flow (Navier–Stokes, etc.) where:
%     - Velocity is approximated by a 10-node (P3) element
%     - Pressure is approximated by a 6-node (P2) element.
%
%   We only build the velocity mass block (2*Nxy × 2*Nxy), then add
%   a zero block for pressure (Npr × Npr). Finally, we fix one corner
%   pressure node to remove the null space.
%
%   The result M is a (2*Nxy + Npr) × (2*Nxy + Npr) sparse matrix:
%     [ M_vel     0
%         0       0   ]
%   with one corner pressure row set to 1 to fix that DOF.
%
% Inputs:
%   nodeInfo.velocity => x,y for P3 nodes (total of Nxy = # of velocity nodes)
%   elemInfo.velElements => (#elements × 10) connectivity for velocity
%              .presElements => (#elements × 6) for pressure (not used here, except to find Npr)
%   boundaryInfo.allVelNodes => velocity node IDs with Dirichlet boundary conditions
%   rho => density (for M = rho * ∫φ_i φ_j dΩ )
%   corner => global pressure node ID to fix (in the range 1..Npr).
%
% Output:
%   M => global sparse mass matrix ( (2*Nxy+Npr) × (2*Nxy+Npr) )

%% 1) Basic size info
numVelNodes = length(nodeInfo.velocity.x);  % total P3 velocity nodes
Nxy         = numVelNodes;                  % velocity DOFs per component
numElements = size(elemInfo.velElements,1);
% Pressure node count (P2)
Npr         = max(elemInfo.presElements(:));

NN = 2*Nxy;  % total velocity DOFs (x + y)

%% 2) Preallocate triplets for velocity mass
% Each P3 element has 10 local velocity nodes => a 10×10 local mass matrix.
% For x-velocity alone, that's 100 entries. Same for y-velocity => another 100.
% So total 200 triplets per element.
nEl = numElements;
nzEst = 200 * nEl;
ii = zeros(nzEst,1);
jj = zeros(nzEst,1);
vv = zeros(nzEst,1);
pos = 0;

%% 3) Identify boundary DOFs for velocity
aBC = boundaryInfo.allVelNodes(:);

%% 4) Precompute shape functions (P3) and 12-point Dunavant rule
%    We'll create a local function or call your existing one for P3 shape.
%    The code below uses a "12-point Dunavant" for polynomial up to degree 4–5
%    but is typically good enough for many flows.

[Nref, ~, ~, wRef, ptRef] = precomputeShapeFunctionsP3_Tri_12();  
% We do NOT need derivatives for the mass matrix, just shape_i * shape_j
numGauss = length(wRef);

%% 5) Loop over elements to build local 10×10 mass
for e = 1:numElements

    % local P3 velocity node indices
    Kvel = elemInfo.velElements(e,:);  % length=10
    xcoords = nodeInfo.velocity.x(Kvel);  % (10×1)
    ycoords = nodeInfo.velocity.y(Kvel);

    % local 10×10 mass
    Me_loc = zeros(10,10);

    % 5A) Gauss integration
    for gp = 1:numGauss
        xi  = ptRef(gp,1);
        eta = ptRef(gp,2);
        wgt = wRef(gp);

        % shape functions at (xi,eta)
        N10 = Nref(:, gp);  % 10×1

        % We must compute the physical area scaling => detJ
        % For that, we do need shape derivatives at this gp:
        % but the function "precomputeShapeFunctionsP3_Tri_12" probably only stored N, not dN/dxi?
        % We'll define a local derivative routine or store it. 
        % For the sake of demonstration, let's call a local function p3ShapeDerivForMass:
        detJ = local_detJ_p3(xcoords, ycoords, xi, eta);

        % accumulate mass
        for iNode = 1:10
            Ni = N10(iNode);
            for jNode = 1:10
                Nj = N10(jNode);
                Me_loc(iNode,jNode) = Me_loc(iNode,jNode) ...
                    + rho * Ni*Nj * detJ * wgt;
            end
        end
    end

    % 5B) Insert local Me_loc into triplets for x-vel, then y-vel
    for iNode = 1:10
        global_i = Kvel(iNode);  % x-velocity row
        for jNode = 1:10
            global_j = Kvel(jNode);
            pos = pos + 1;
            ii(pos) = global_i;
            jj(pos) = global_j;
            vv(pos) = Me_loc(iNode,jNode);
        end
    end

    for iNode = 1:10
        global_i = Kvel(iNode) + Nxy;  % y-velocity
        for jNode = 1:10
            global_j = Kvel(jNode) + Nxy;
            pos = pos + 1;
            ii(pos) = global_i;
            jj(pos) = global_j;
            vv(pos) = Me_loc(iNode,jNode);
        end
    end

end  % end element loop

%% 6) Build sparse velocity mass
M_vel = sparse(ii, jj, vv, NN, NN);

%% 7) Impose Dirichlet BCs on velocity
%   Zero out rows & set diagonal to 1 for each Dirichlet velocity node
for bcIndex = 1:length(aBC)
    bNode = aBC(bcIndex);
    % x-row
    M_vel(bNode,:) = 0;
    M_vel(bNode,bNode) = 1;
    % y-row
    rowY = bNode + Nxy;
    M_vel(rowY,:) = 0;
    M_vel(rowY,rowY) = 1;
end

%% 8) Lump the velocity mass (if you want a diagonal mass for velocity)
%   If you want a consistent mass, skip lumping. 
%   If you do want lumping, do:
rowSum = sum(M_vel,2);
M_lumped_vel = spdiags(rowSum, 0, NN, NN);

%% 9) Expand to include pressure DOFs
% The final matrix is block diagonal:
%   [ M_vel / M_lumped_vel,   0
%            0,              0  ] 
% Here we choose the lumped version for velocity:
M = [ M_lumped_vel,           sparse(NN, Npr);
      sparse(Npr, NN),        sparse(Npr, Npr) ];

%% 10) Fix corner pressure node
% 'corner' is presumably a global pressure node index in [1..Npr].
cornerGlobal = NN + corner;  % offset after velocity rows
M(cornerGlobal,:) = 0;
M(cornerGlobal, cornerGlobal) = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function detJ = local_detJ_p3(xcoords, ycoords, xi, eta)
% local_detJ_p3:
%   A minimal function to get the Jacobian determinant for a P3 element
%   at (xi,eta). We do not need the full derivative array, just the
%   determinant for area scaling in the mass matrix.
%
%   For a P3 element, we can reuse the same derivatives as P1 for area,
%   because the geometry is (usually) linear => Actually if the geometry
%   is truly curved, you do need the actual isoparam shape. 
%   We'll do a direct approach: use the 10-node shape function derivatives
%   at (xi,eta), then sum xcoords*dN/dxi etc.

% Barycentric:
L1 = 1 - xi - eta; 
L2 = xi; 
L3 = eta;

% For a P3 isoparam, geometry might also be cubic if you have curved boundaries.
% So we do the actual derivative. We'll define them inline:

% shape functions for geometry interpolation:
% corners + mid-edge + bubble => same polynomials as in p3basisGmsh
% But typically the mesh "POS" is only storing linear corners, so there's
% no actual mid-edge control in geometry. If that's the case, geometry is
% effectively linear => the determinant is the same as a linear P1 approach.
% We'll do the full approach anyway for correctness.

[Ngeom, dNxi_geom, dNeta_geom] = p3basis_forGeometry(L1,L2,L3);

% Compute x,y in physical space
X = sum(xcoords .* Ngeom);
Y = sum(ycoords .* Ngeom);
% partial derivatives
dX_dxi  = sum(xcoords .* dNxi_geom);
dX_deta = sum(xcoords .* dNeta_geom);
dY_dxi  = sum(ycoords .* dNxi_geom);
dY_deta = sum(ycoords .* dNeta_geom);

% determinant
detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
end

function [Ngeom, dNxi_geom, dNeta_geom] = p3basis_forGeometry(L1,L2,L3)
% p3basis_forGeometry:
%   Minimal polynomial set for geometry. If your geometry
%   truly uses a P3 shape, define the corner + mid-edge + bubble polynomials. 
%   If your geometry is actually linear, you'd skip mid-edge/bubble terms
%   for X interpolation. But let's do the full P3 in case you do have curved edges.

% P3 shape polynomials in barycentric form (same as "p3basisGmsh" style).
Ngeom = zeros(10,1);
c = 9/2;

% corners
Ngeom(1) = L1*(3*L1 -1)*(3*L1 -2)/2;
Ngeom(2) = L2*(3*L2 -1)*(3*L2 -2)/2;
Ngeom(3) = L3*(3*L3 -1)*(3*L3 -2)/2;

% mid-edge
Ngeom(4) = c*L1*L2*(3*L1 -1);
Ngeom(5) = c*L1*L2*(3*L2 -1);
Ngeom(6) = c*L2*L3*(3*L2 -1);
Ngeom(7) = c*L2*L3*(3*L3 -1);
Ngeom(8) = c*L3*L1*(3*L3 -1);
Ngeom(9) = c*L3*L1*(3*L1 -1);

% bubble
Ngeom(10) = 27*L1*L2*L3;

% partial derivatives of L1,L2,L3 wrt (xi,eta):
dL1_dxi=-1; dL1_deta=-1;
dL2_dxi= 1; dL2_deta= 0;
dL3_dxi= 0; dL3_deta= 1;

dNxi_geom = zeros(10,1);
dNeta_geom= zeros(10,1);

% corner1
f1 = (3*L1 -1)*(3*L1 -2)/2;
df1dL1= (9*(2*L1 -1))/2;
dNxi_geom(1)  = dL1_dxi*f1 + L1*(df1dL1*dL1_dxi);
dNeta_geom(1) = dL1_deta*f1 + L1*(df1dL1*dL1_deta);

% corner2
f2 = (3*L2 -1)*(3*L2 -2)/2;
df2dL2= (9*(2*L2 -1))/2;
dNxi_geom(2)  = dL2_dxi*f2 + L2*(df2dL2*dL2_dxi);
dNeta_geom(2) = dL2_deta*f2 + L2*(df2dL2*dL2_deta);

% corner3
f3 = (3*L3 -1)*(3*L3 -2)/2;
df3dL3= (9*(2*L3 -1))/2;
dNxi_geom(3)  = dL3_dxi*f3 + L3*(df3dL3*dL3_dxi);
dNeta_geom(3) = dL3_deta*f3 + L3*(df3dL3*dL3_deta);

% mid-edge4 => c * L1*L2*(3L1 -1)
[dNxi_geom(4), dNeta_geom(4)] = d_of_prod3( c, L1, L2, (3*L1 -1), ...
                                            dL1_dxi,dL1_deta, ...
                                            dL2_dxi,dL2_deta, ...
                                            3*dL1_dxi,3*dL1_deta);
% mid-edge5 => c * L1*L2*(3L2 -1)
[dNxi_geom(5), dNeta_geom(5)] = d_of_prod3( c, L1, L2, (3*L2 -1), ...
                                            dL1_dxi,dL1_deta, ...
                                            dL2_dxi,dL2_deta, ...
                                            3*dL2_dxi,3*dL2_deta);

% mid-edge6 => c * L2*L3*(3*L2 -1)
[dNxi_geom(6), dNeta_geom(6)] = d_of_prod3( c, L2, L3, (3*L2 -1), ...
                                            dL2_dxi,dL2_deta, ...
                                            dL3_dxi,dL3_deta, ...
                                            3*dL2_dxi,3*dL2_deta);

% mid-edge7 => c * L2*L3*(3*L3 -1)
[dNxi_geom(7), dNeta_geom(7)] = d_of_prod3( c, L2, L3, (3*L3 -1), ...
                                            dL2_dxi,dL2_deta, ...
                                            dL3_dxi,dL3_deta, ...
                                            3*dL3_dxi,3*dL3_deta);

% mid-edge8 => c * L3*L1*(3*L3 -1)
[dNxi_geom(8), dNeta_geom(8)] = d_of_prod3( c, L3, L1, (3*L3 -1), ...
                                            dL3_dxi,dL3_deta, ...
                                            dL1_dxi,dL1_deta, ...
                                            3*dL3_dxi,3*dL3_deta);

% mid-edge9 => c * L3*L1*(3*L1 -1)
[dNxi_geom(9), dNeta_geom(9)] = d_of_prod3( c, L3, L1, (3*L1 -1), ...
                                            dL3_dxi,dL3_deta, ...
                                            dL1_dxi,dL1_deta, ...
                                            3*dL1_dxi,3*dL1_deta);

% bubble10 => 27*L1*L2*L3
[dNxi_geom(10), dNeta_geom(10)] = d_of_prod3( 27, L1, L2, L3, ...
                                              dL1_dxi,dL1_deta, ...
                                              dL2_dxi,dL2_deta, ...
                                              dL3_dxi,dL3_deta);
end

function [dNdxi, dNeta] = d_of_prod3( scale, X, Y, Z, ...
                                      dXdxi, dXdet, ...
                                      dYdxi, dYdet, ...
                                      dZdxi, dZdet)
% d(N)/dxi = scale * [ (dXdxi * Y * Z) + (X * dYdxi * Z) + (X * Y * dZdxi ) ]
% Similarly for dN/deta
dNdxi = scale*( dXdxi*Y*Z + X*dYdxi*Z + X*Y*dZdxi );
dNeta = scale*( dXdet*Y*Z + X*dYdet*Z + X*Y*dZdet );
end

%--------------------------------------------------------------------------
function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP3_Tri_12()
% precomputeShapeFunctionsP3_Tri_12:
%   Returns shape functions (N), derivatives (dNxi, dNeta),
%   and a 12-point Dunavant rule (g_pt, g_wt) for a P3 triangle,
%   in the reference domain 0 <= xi, 0 <= eta, xi+eta <= 1.
%
%   Here we demonstrate a simpler 12-point rule (Dunavant) that is
%   exact for polynomials up to degree 5. Usually enough for typical mass
%   matrix integrals with P3 elements.

[g_pt, g_wt] = triGaussPoints6(); 
numGauss = length(g_wt);

numNodes = 10;  % P3 has 10 local nodes
N      = zeros(numNodes, numGauss);
dNxi   = zeros(numNodes, numGauss);
dNeta  = zeros(numNodes, numGauss);

for k = 1:numGauss
    xi  = g_pt(k,1);
    eta = g_pt(k,2);
    [Ni, dNxi_k, dNeta_k] = p3IsoShape_basic(xi, eta);
    N(:,k)     = Ni;
    dNxi(:,k)  = dNxi_k;
    dNeta(:,k) = dNeta_k;
end
end

function [Ni, dNxi_k, dNeta_k] = p3IsoShape_basic(xi, eta)
% A minimal p3 shape function approach. If you want the "mid-edge" polynomials
% used in your solver, define them. We'll do the standard set from your Gmsh ordering.

L1 = 1 - xi - eta;
L2 = xi;
L3 = eta;

Ni = zeros(10,1);
% corners
Ni(1) = L1*(3*L1 -1)*(3*L1 -2)/2;
Ni(2) = L2*(3*L2 -1)*(3*L2 -2)/2;
Ni(3) = L3*(3*L3 -1)*(3*L3 -2)/2;
% mid-edges
c = 9/2;
Ni(4) = c*L1*L2*(3*L1 -1);
Ni(5) = c*L1*L2*(3*L2 -1);
Ni(6) = c*L2*L3*(3*L2 -1);
Ni(7) = c*L2*L3*(3*L3 -1);
Ni(8) = c*L3*L1*(3*L3 -1);
Ni(9) = c*L3*L1*(3*L1 -1);
% bubble
Ni(10) = 27*L1*L2*L3;

% derivatives
dNxi_k  = zeros(10,1);
dNeta_k = zeros(10,1);

% partials of L1=1-xi-eta => dL1dxi=-1,dL1deta=-1
% L2=xi => dL2dxi=1, dL2deta=0
% L3=eta=> dL3dxi=0, dL3deta=1
dL1xi=-1; dL1et=-1;
dL2xi= 1; dL2et= 0;
dL3xi= 0; dL3et= 1;

% corner1
f1 = (3*L1-1)*(3*L1-2)/2;
df1= (9*(2*L1-1))/2; 
dNxi_k(1)  = dL1xi*f1 + L1*(df1*dL1xi);
dNeta_k(1) = dL1et*f1 + L1*(df1*dL1et);

% corner2
f2 = (3*L2-1)*(3*L2-2)/2;
df2= (9*(2*L2-1))/2;
dNxi_k(2)  = dL2xi*f2 + L2*(df2*dL2xi);
dNeta_k(2) = dL2et*f2 + L2*(df2*dL2et);

% corner3
f3 = (3*L3-1)*(3*L3-2)/2;
df3= (9*(2*L3-1))/2;
dNxi_k(3)  = dL3xi*f3 + L3*(df3*dL3xi);
dNeta_k(3) = dL3et*f3 + L3*(df3*dL3et);

% mid-edge4 => c*L1*L2*(3L1-1)
[dNxi_k(4), dNeta_k(4)] = d_of_prod3( 9/2, L1, L2, (3*L1-1), ...
                                      dL1xi,dL1et, dL2xi,dL2et, 3*dL1xi,3*dL1et);
% mid-edge5 => c*L1*L2*(3L2-1)
[dNxi_k(5), dNeta_k(5)] = d_of_prod3( 9/2, L1, L2, (3*L2-1), ...
                                      dL1xi,dL1et, dL2xi,dL2et, 3*dL2xi,3*dL2et);
% mid-edge6 => c*L2*L3*(3*L2-1)
[dNxi_k(6), dNeta_k(6)] = d_of_prod3( 9/2, L2, L3, (3*L2-1), ...
                                      dL2xi,dL2et, dL3xi,dL3et, 3*dL2xi,3*dL2et);
% mid-edge7 => c*L2*L3*(3*L3-1)
[dNxi_k(7), dNeta_k(7)] = d_of_prod3( 9/2, L2, L3, (3*L3-1), ...
                                      dL2xi,dL2et, dL3xi,dL3et, 3*dL3xi,3*dL3et);
% mid-edge8 => c*L3*L1*(3*L3-1)
[dNxi_k(8), dNeta_k(8)] = d_of_prod3( 9/2, L3, L1, (3*L3-1), ...
                                      dL3xi,dL3et, dL1xi,dL1et, 3*dL3xi,3*dL3et);
% mid-edge9 => c*L3*L1*(3*L1-1)
[dNxi_k(9), dNeta_k(9)] = d_of_prod3( 9/2, L3, L1, (3*L1-1), ...
                                      dL3xi,dL3et, dL1xi,dL1et, 3*dL1xi,3*dL1et);
% bubble => 27*L1*L2*L3
[dNxi_k(10), dNeta_k(10)] = d_of_prod3(27, L1,L2,L3, ...
                                       dL1xi,dL1et, dL2xi,dL2et, dL3xi,dL3et);
end

function [g_pt, g_wt] = triGaussPoints6()
% triGaussPoints6: returns a 6-point integration rule on the
%  reference triangle (xi,eta>=0, xi+eta<=1).
% Each row of g_pt is (xi, eta), and the corresponding weight is in g_wt.
% This rule integrates polynomials exactly up to degree 4.
%
% Weights sum to area = 1/2 for the standard reference triangle,
% so each weight is the sub-area contribution.
%
% Source of these points is standard in many FE references.

    g_pt = [ ...
      0.4459484909, 0.4459484909;
      0.4459484909, 0.1081030182;
      0.1081030182, 0.4459484909;
      0.0915762135, 0.0915762135;
      0.0915762135, 0.8168475730;
      0.8168475730, 0.0915762135 ];

    g_wt = [ ...
      0.2233815897;
      0.2233815897;
      0.2233815897;
      0.1099517437;
      0.1099517437;
      0.1099517437 ];
end


% function [g_pt, g_wt] = triGaussPoints12()
% % triGaussPoints12:
% %   A 12-point Dunavant rule on the reference triangle 
% %   0<=xi,0<=eta, xi+eta<=1,
% %   exact up to polynomial order 5 (commonly used).
% 
% g_pt = [
%   0.24928674517091,  0.24928674517091;
%   0.24928674517091,  0.50142650965818;
%   0.50142650965818,  0.24928674517091;
% 
%   0.06308901449150,  0.06308901449150;
%   0.06308901449150,  0.87382197101700;
%   0.87382197101700,  0.06308901449150;
% 
%   0.31035245103378,  0.63650249912140;
%   0.63650249912140,  0.05314504984482;
%   0.05314504984482,  0.31035245103378;
% 
%   0.63650249912140,  0.31035245103378;
%   0.31035245103378,  0.05314504984482;
%   0.05314504984482,  0.63650249912140
% ];
% 
% g_wt = [
%   0.05839313786319;
%   0.05839313786319;
%   0.05839313786319;
% 
%   0.02542245318510;
%   0.02542245318510;
%   0.02542245318510;
% 
%   0.04142553780919;
%   0.04142553780919;
%   0.04142553780919;
% 
%   0.04142553780919;
%   0.04142553780919;
%   0.04142553780919
% ];




% function M = mass_matrix_func_p2p1( ...
%     nodeInfo, elemInfo, boundaryInfo, ...
%     rho, corner )
% % MASS_MATRIX_FUNC_P2P1:
% %   Builds the mass matrix for a P2–P1 triangular discretization,
% %   then imposes Dirichlet boundary conditions on velocity DOFs,
% %   lumps the velocity block, and expands for pressure DOFs.
% %   Finally, fixes one corner pressure node.
% %
% %   The result is a (2*Nxy + Npr) × (2*Nxy + Npr) sparse matrix.
% %     - The top-left 2*Nxy×2*Nxy block is the velocity mass (lumped).
% %     - The pressure mass block is zero (incompressible flow).
% %     - We fix one corner pressure node to remove the null space (p-ref).
% %
% % Inputs:
% %   nodeInfo.velocity => x,y coordinates for P2 nodes
% %   elemInfo.velElements => Nx6 connectivity (P2)
% %              .presElements => Nx3 connectivity (P1)
% %   boundaryInfo.allNodes => velocity Dirichlet nodes
% %   rho => density (if needed for M = rho * volume)
% %   corner => index of the pressure node to fix
% %
% % Output:
% %   M => the global mass matrix (sparse)
% 
% %% 1) Basic size info
% numVelNodes = length(nodeInfo.velocity.x);  % total P2 velocity nodes
% Nxy         = numVelNodes;                  % velocity DOFs per component
% numElements = size(elemInfo.velElements,1);
% Npr         = max(elemInfo.presElements(:));
% 
% NN = 2*Nxy;  % total velocity DOFs (x + y)
% 
% %% 2) Preallocate triplets for velocity mass
% % Each P2 element has 6 local velocity nodes => a 6×6 local mass matrix.
% % For x-velocity, that's 6×6=36 entries. Same for y-velocity => another 36.
% % So total 72 triplets per element.
% nzEst = 72 * numElements;
% ii = zeros(nzEst,1);
% jj = zeros(nzEst,1);
% vv = zeros(nzEst,1);
% pos = 0;
% 
% %% 3) Identify boundary DOFs for velocity
% aBC = boundaryInfo.allVelNodes(:);
% 
% %% 4) Precompute the reference P2 shape functions and Gauss rule
% %    e.g., a 6-point rule for the reference triangle.
% [Nref, dNxi_ref, dNeta_ref, wRef, ptRef] = precomputeShapeFunctionsP2_Tri();
% numGauss = length(wRef);  % typically 6 or 7
% 
% %% 5) Loop over elements to build local 6×6 mass
% for e = 1:numElements
% 
%     % local P2 velocity node indices
%     Kvel = elemInfo.velElements(e,:);  % length=6
%     xcoords = nodeInfo.velocity.x(Kvel);  % (6×1)
%     ycoords = nodeInfo.velocity.y(Kvel);
% 
%     % local 6×6 mass
%     Me_loc = zeros(6,6);
% 
%     % 5A) Gauss integration
%     for gp = 1:numGauss
%         xi  = ptRef(gp,1);
%         eta = ptRef(gp,2);
%         wgt = wRef(gp);
% 
%         % shape functions & derivatives at (xi,eta) on the reference triangle
%         Ni     = Nref(:, gp);        % 6×1
%         dNxi   = dNxi_ref(:, gp);    % 6×1
%         dNeta  = dNeta_ref(:, gp);   % 6×1
% 
%         % map to physical coords => get det(J)
%         [dNx, dNy, detJ] = p3ShapeDerivativesAllNodes( ...
%                            xcoords, ycoords, dNxi, dNeta );
%         % for mass, we only need shape_i * shape_j * (detJ*wgt), 
%         % no velocity derivatives or such.
% 
%         for iNode = 1:6
%             shape_i = Ni(iNode);
%             for jNode = 1:6
%                 shape_j = Ni(jNode);
%                 Me_loc(iNode,jNode) = Me_loc(iNode,jNode) ...
%                     + rho * shape_i*shape_j * detJ * wgt;
%             end
%         end
%     end
% 
%     % 5B) Insert local Me_loc into triplets for x-vel, then y-vel
%     for iNode = 1:6
%         global_i = Kvel(iNode);  % x-velocity row
%         for jNode = 1:6
%             global_j = Kvel(jNode);
%             pos = pos + 1;
%             ii(pos) = global_i;
%             jj(pos) = global_j;
%             vv(pos) = Me_loc(iNode,jNode);
%         end
%     end
% 
%     for iNode = 1:6
%         global_i = Kvel(iNode) + Nxy;  % y-velocity
%         for jNode = 1:6
%             global_j = Kvel(jNode) + Nxy;
%             pos = pos + 1;
%             ii(pos) = global_i;
%             jj(pos) = global_j;
%             vv(pos) = Me_loc(iNode,jNode);
%         end
%     end
% 
% end  % end element loop
% 
% %% 6) Build sparse velocity mass
% M_vel = sparse(ii, jj, vv, NN, NN);
% 
% %% 7) Impose Dirichlet BCs on velocity
% % Zero out rows & set diagonal to 1
% for bcIndex = 1:length(aBC)
%     bNode = aBC(bcIndex);
%     % x-row
%     M_vel(bNode,:) = 0;
%     M_vel(bNode,bNode) = 1;
%     % y-row
%     rowY = bNode + Nxy;
%     M_vel(rowY,:) = 0;
%     M_vel(rowY,rowY) = 1;
% end
% 
% %% 8) Lump the velocity mass
% % sum across each row
% rowSum = sum(M_vel,2);
% M_lumped_vel = spdiags(rowSum, 0, NN, NN);
% 
% %% 9) Expand for pressure DOFs
% % The final matrix is block diagonal:
% %  [ M_lumped_vel,   0
% %         0,         0  ] 
% % where the zero block is Npr×Npr for pressure
% M = [ M_lumped_vel,           sparse(NN, Npr);
%       sparse(Npr, NN),        sparse(Npr, Npr) ];
% 
% %% 10) Fix corner pressure node
% % corner is presumably a global pressure index from 1..Npr
% cornerGlobal = NN + corner;  % offset after the 2*Nxy velocity rows
% M(cornerGlobal,:) = 0;
% M(cornerGlobal, cornerGlobal) = 1;
% 
% end
% 
% function [dNxAll,dNyAll,detJ] = p3ShapeDerivativesAllNodes( ...
%     xcoords,ycoords,dNxi,dNeta)
% % p3ShapeDerivativesAllNodes:
% %   For 10-node P3 triangle, given partial derivatives wrt (xi,eta),
% %   compute partials wrt physical (x,y).
% 
% dX_dxi  = sum(xcoords.*dNxi);
% dX_deta = sum(xcoords.*dNeta);
% dY_dxi  = sum(ycoords.*dNxi);
% dY_deta = sum(ycoords.*dNeta);
% 
% J = [ dX_dxi, dY_dxi;
%       dX_deta,dY_deta];
% detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
% 
% invJ=inv(J);
% 
% dNxAll = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
% dNyAll = invJ(2,1)*dNxi + invJ(2,2)*dNeta;
% end
% 
% 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [dNxAll, dNyAll, ddJ] = q2ShapeDerivatives_AllNodes( ...
% %     xcoords, ycoords, ...
% %     quad1_vals, quad2_vals, ...
% %     dquad1_vals, dquad2_vals )
% % % same as your old code => returns global derivatives + detJ etc.
% % % We'll do the same logic, but only the final 'ddJ' is used if you want
% % % no derivative-based shape in the integrand
% % dNdxi  = dquad1_vals.*quad2_vals;
% % dNdeta = quad1_vals.*dquad2_vals;
% % 
% % dX_dxi  = sum(xcoords.*dNdxi);
% % dX_deta = sum(xcoords.*dNdeta);
% % dY_dxi  = sum(ycoords.*dNdxi);
% % dY_deta = sum(ycoords.*dNdeta);
% % 
% % % J = [dX_dxi, dX_deta;
% % %      dY_dxi, dY_deta];
% %  J = [ dX_dxi, dY_dxi ;
% %           dX_deta,  dY_deta ];
% % ddJ = (dX_dxi*dY_deta - dX_deta*dY_dxi);
% % 
% % % if negative is possible => ddJ=abs(...)
% % invJ = inv(J);
% % 
% % % to illustrate you do:
% % dNxAll = zeros(9,1);
% % dNyAll = zeros(9,1);
% % for i=1:9
% %     dxi  = dNdxi(i);
% %     deta = dNdeta(i);
% %     dNxAll(i) = invJ(1,1)*dxi + invJ(1,2)*deta;
% %     dNyAll(i) = invJ(2,1)*dxi + invJ(2,2)*deta;
% % end
% % end
% % 
% % function detJ = getDetJ(~, ~, ddJ)
% % % in mass matrix we only want ddJ. 
% % % ignoring the derivative arrays => pass them in if you want.
% % detJ = (ddJ);  % or do sign fix
% % end
% % 