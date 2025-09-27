function K = build_stiffness_operator_stationary( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    quadratic1, quadratic2, del_quad1, del_quad2, ...
    lin1, lin2, ...
    g_wt, ...
    D, Re, o, U1, U2, U3, gamma, mu, rho,corner )
%BUILD_STIFFNESS_OPERATOR_S3
%  Assembles the Q2–Q1 Navier–Stokes stiffness matrix using element 
%  connectivity from 'elemInfo' and node coordinates from 'nodeInfo'.
%
%  This version uses triplet-based assembly, vectorized BC imposition,
%  and FLATTENS shapeAll/dNxAll/dNyAll from (9x5x5) to (9x25) to reduce
%  3D-slice overhead.

%% 1) Basic DOF counts
numVelNodes = length(nodeInfo.velocity.x);  % total Q2 velocity nodes
Nxy = numVelNodes;                          % velocity DOFs per component
numElements = size(elemInfo.quadElements,1);% how many Q2 elements
Npr = max(elemInfo.presElements(:));        % Q1 pressure node count

%% 2) Pre-allocate triplet arrays for each sub-block
n81 = 81*numElements;  % For 9x9 blocks
n36 = 36*numElements;  % For 9x4 or 4x9 blocks

% Jx1: (Nxy x Nxy)
I_Jx1 = zeros(n81,1);  J_Jx1 = zeros(n81,1);  S_Jx1 = zeros(n81,1);
ptr_Jx1 = 1;

% Jx2: (Nxy x Nxy)
I_Jx2 = zeros(n81,1);  J_Jx2 = zeros(n81,1);  S_Jx2 = zeros(n81,1);
ptr_Jx2 = 1;

% Jx3: (Nxy x Npr) -> 9x4 = 36 entries/element
I_Jx3 = zeros(n36,1);  J_Jx3 = zeros(n36,1);  S_Jx3 = zeros(n36,1);
ptr_Jx3 = 1;

% Jy1: (Nxy x Nxy)
I_Jy1 = zeros(n81,1);  J_Jy1 = zeros(n81,1);  S_Jy1 = zeros(n81,1);
ptr_Jy1 = 1;

% Jy2: (Nxy x Nxy)
I_Jy2 = zeros(n81,1);  J_Jy2 = zeros(n81,1);  S_Jy2 = zeros(n81,1);
ptr_Jy2 = 1;

% Jy3: (Nxy x Npr) -> 9x4 = 36 entries/element
I_Jy3 = zeros(n36,1);  J_Jy3 = zeros(n36,1);  S_Jy3 = zeros(n36,1);
ptr_Jy3 = 1;

% Jc1: (Npr x Nxy) -> 4x9=36 entries/element
I_Jc1 = zeros(n36,1);  J_Jc1 = zeros(n36,1);  S_Jc1 = zeros(n36,1);
ptr_Jc1 = 1;

% Jc2: (Npr x Nxy)
I_Jc2 = zeros(n36,1);  J_Jc2 = zeros(n36,1);  S_Jc2 = zeros(n36,1);
ptr_Jc2 = 1;

% Jc3 is always zero for incompressibility
Jc3 = spalloc(Npr, Npr, 0);

%% 3) Identify boundary velocity DOFs
aBC = boundaryInfo.allNodes;  % Dirichlet BC for these global velocity nodes

%% 4) Element loop
for e = 1:numElements
    
    % local velocity node indices (9 Q2 nodes)
    Kvel = elemInfo.quadElements(e,:);
    % local pressure node indices (4 Q1 nodes)
    Kpr  = elemInfo.presElements(e,:);

    xcoords = nodeInfo.velocity.x(Kvel);
    ycoords = nodeInfo.velocity.y(Kvel);

    % Local sub-blocks
    Kx1e = zeros(9,9);
    Kx2e = zeros(9,9);
    Kx3e = zeros(9,4);

    Ky1e = zeros(9,9);
    Ky2e = zeros(9,9);
    Ky3e = zeros(9,4);

    Kc1e = zeros(4,9);
    Kc2e = zeros(4,9);
    % Kc3e = zeros(4,4);  % typically 0 for incompressible

    %-----------------------------------------------------------------
    % 4A) Precompute shape derivatives + shapeAll in 3D, then flatten
    %     shapeAll(9,5,5), dNxAll(9,5,5), etc.
    %-----------------------------------------------------------------
    shapeAll    = zeros(9,5,5);
    dNxAll      = zeros(9,5,5);
    dNyAll      = zeros(9,5,5);
    detJgauss   = zeros(5,5);

    for gd = 1:5
        for ga = 1:5
            for kNode=1:9
                shapeAll(kNode, gd, ga) = ...
                    quadratic1(kNode, gd) * quadratic2(kNode, ga);
            end

            [dNx, dNy, ddJ] = q2ShapeDerivatives_AllNodes( ...
                xcoords, ycoords, ...
                quadratic1(:,gd), quadratic2(:,ga), ...
                del_quad1(:,gd), del_quad2(:,ga) );

            dNxAll(:,gd,ga) = dNx;
            dNyAll(:,gd,ga) = dNy;
            detJgauss(gd,ga) = ddJ;
        end
    end

    %==== FLATTEN to 2D arrays => (9 x 25), (1 x 25) =====
    shapeAll2 = reshape(shapeAll, [9, 25]);
    dNxAll2   = reshape(dNxAll,   [9, 25]);
    dNyAll2   = reshape(dNyAll,   [9, 25]);
    detJall   = reshape(detJgauss,[1, 25]);

    % Local copies of U1, U2 at the 9 velocity nodes
    U1el = U1(Kvel);
    U2el = U2(Kvel);

    %% 4B) Integrate for the velocity-velocity blocks
    for iNode = 1:9
        for kNode = 1:9

            jxa_val = 0; jxb_val = 0;
            jya_val = 0; jyb_val = 0;
            jca_val = 0; jcb_val = 0;

            %---- Loop 1: Diffusion + (any other terms not requiring a1,a2)
            for gd = 1:5
                for ga = 1:5
                    idx = 5*(gd-1) + ga;  % flatten (gd,ga) => 1..25
                    wt   = g_wt(gd)*g_wt(ga);
                    detJ = detJall(idx);

                    dNx_i = dNxAll2(iNode, idx);
                    dNy_i = dNyAll2(iNode, idx);
                    dNx_k = dNxAll2(kNode, idx);
                    dNy_k = dNyAll2(kNode, idx);

                    % Viscous terms for x-momentum wrt x-velocity
                    visc_xx = mu * ( 2*(dNx_i*dNx_k) + (dNy_i*dNy_k) ) * detJ * wt;
                    jxa_val = jxa_val + visc_xx;
                    graddiv_xx = gamma*(dNx_i*dNx_k)*detJ*wt;
                    jxa_val    = jxa_val + graddiv_xx;

                    % Viscous terms x-momentum wrt y-velocity
                    visc_xy = mu*( dNy_i*dNx_k )*detJ*wt;
                    jxb_val = jxb_val + visc_xy;
                    graddiv_xy = gamma*(dNx_i*dNy_k)*detJ*wt;
                    jxb_val    = jxb_val + graddiv_xy;

                    % Viscous terms y-momentum wrt x-velocity
                    visc_yx = mu*( dNx_i*dNy_k )*detJ*wt;
                    jya_val = jya_val + visc_yx;
                    graddiv_yx = gamma*(dNy_i*dNx_k)*detJ*wt;
                    jya_val    = jya_val + graddiv_yx;

                    % Viscous terms y-momentum wrt y-velocity
                    visc_yy = mu*( (dNx_i*dNx_k) + 2*(dNy_i*dNy_k) )*detJ*wt;
                    jyb_val = jyb_val + visc_yy;
                    graddiv_yy = gamma*(dNy_i*dNy_k)*detJ*wt;
                    jyb_val    = jyb_val + graddiv_yy;

                    % Continuity wrt U?
                    if iNode <= 4  % first 4 Q2 nodes => Q1 corners
                        cont_cx = dNx_k * lin1(iNode,gd) * lin2(iNode,ga) * detJ * wt;
                        jca_val = jca_val - cont_cx;

                        cont_cy = dNy_k * lin1(iNode,gd) * lin2(iNode,ga) * detJ * wt;
                        jcb_val = jcb_val - cont_cy;
                    end
                end
            end

            %---- Loop 2: Convection terms (a1,a2)
            for gd = 1:5
                for ga = 1:5
                    idx = 5*(gd-1) + ga;
                    wt   = g_wt(gd)*g_wt(ga);
                    detJ = detJall(idx);

                    shape_i = shapeAll2(iNode, idx);
                    dNx_k   = dNxAll2(kNode, idx);
                    dNy_k   = dNyAll2(kNode, idx);

                    % a1, a2 from local velocity
                    shVec = shapeAll2(:, idx);  % 9 x 1
                    a1 = shVec' * U1el;  
                    a2 = shVec' * U2el;

                    conv_x = (a1*dNx_k + a2*dNy_k);
                    jxa_val = jxa_val + rho*Re*( shape_i * conv_x ) * detJ * wt;

                    conv_y = (a1*dNx_k + a2*dNy_k);
                    jyb_val = jyb_val + rho*Re*( shape_i * conv_y ) * detJ * wt;
                end
            end

            % Store in local blocks
            Kx1e(iNode,kNode) = jxa_val;
            Kx2e(iNode,kNode) = jxb_val;
            Ky1e(iNode,kNode) = jya_val;
            Ky2e(iNode,kNode) = jyb_val;

            if iNode <= 4
                Kc1e(iNode,kNode) = jca_val;
                Kc2e(iNode,kNode) = jcb_val;
            end
        end
    end

    %% 4C) velocity--pressure blocks (9x4 => Kx3e, Ky3e) 
    for iNode = 1:9
        for pNode = 1:4
            jxc_val=0; jyc_val=0;
            for d=1:5
                for a=1:5
                    idx = 5*(d-1) + a;  % flatten
                    wt = g_wt(d)*g_wt(a);
                    detJ = detJall(idx);

                    dNx_i = dNxAll2(iNode, idx);
                    dNy_i = dNyAll2(iNode, idx);

                    jxc_val = jxc_val - dNx_i * lin1(pNode,d) * lin2(pNode,a) * detJ * wt;
                    jyc_val = jyc_val - dNy_i * lin1(pNode,d) * lin2(pNode,a) * detJ * wt;
                end
            end
            Kx3e(iNode,pNode) = jxc_val;
            Ky3e(iNode,pNode) = jyc_val;
        end
    end

    %% 4D) Accumulate into triplet arrays
    % Jx1: 9x9 block
    [iMat, jMat] = ndgrid(Kvel, Kvel);   % 9x9 => 81 entries
    rng = ptr_Jx1:(ptr_Jx1+81-1);
    I_Jx1(rng) = iMat(:); 
    J_Jx1(rng) = jMat(:);
    S_Jx1(rng) = Kx1e(:);
    ptr_Jx1 = ptr_Jx1 + 81;

    % Jx2
    rng = ptr_Jx2:(ptr_Jx2+81-1);
    I_Jx2(rng) = iMat(:);
    J_Jx2(rng) = jMat(:);
    S_Jx2(rng) = Kx2e(:);
    ptr_Jx2 = ptr_Jx2 + 81;

    % Jx3: 9x4 => 36
    [iMat, jMat] = ndgrid(Kvel, Kpr);
    rng = ptr_Jx3:(ptr_Jx3+36-1);
    I_Jx3(rng) = iMat(:);
    J_Jx3(rng) = jMat(:);
    S_Jx3(rng) = Kx3e(:);
    ptr_Jx3 = ptr_Jx3 + 36;

    % Jy1
    rng = ptr_Jy1:(ptr_Jy1+81-1);
    [iMat2, jMat2] = ndgrid(Kvel, Kvel);
    I_Jy1(rng) = iMat2(:);
    J_Jy1(rng) = jMat2(:);
    S_Jy1(rng) = Ky1e(:);
    ptr_Jy1 = ptr_Jy1 + 81;

    % Jy2
    rng = ptr_Jy2:(ptr_Jy2+81-1);
    I_Jy2(rng) = iMat2(:);
    J_Jy2(rng) = jMat2(:);
    S_Jy2(rng) = Ky2e(:);
    ptr_Jy2 = ptr_Jy2 + 81;

    % Jy3: 9x4 => 36
    [iMat3, jMat3] = ndgrid(Kvel, Kpr);
    rng = ptr_Jy3:(ptr_Jy3+36-1);
    I_Jy3(rng) = iMat3(:);
    J_Jy3(rng) = jMat3(:);
    S_Jy3(rng) = Ky3e(:);
    ptr_Jy3 = ptr_Jy3 + 36;

    % Jc1: (4x9 => 36). row = Kpr, col = Kvel
    [iMat4, jMat4] = ndgrid(Kpr, Kvel);
    rng = ptr_Jc1:(ptr_Jc1+36-1);
    I_Jc1(rng) = iMat4(:);
    J_Jc1(rng) = jMat4(:);
    S_Jc1(rng) = Kc1e(:);
    ptr_Jc1 = ptr_Jc1 + 36;

    % Jc2: (4x9 => 36)
    rng = ptr_Jc2:(ptr_Jc2+36-1);
    I_Jc2(rng) = iMat4(:);
    J_Jc2(rng) = jMat4(:);
    S_Jc2(rng) = Kc2e(:);
    ptr_Jc2 = ptr_Jc2 + 36;
end

%% 5) Build the sparse sub-blocks just once
Jx1 = sparse(I_Jx1, J_Jx1, S_Jx1, Nxy, Nxy);
Jx2 = sparse(I_Jx2, J_Jx2, S_Jx2, Nxy, Nxy);
Jx3 = sparse(I_Jx3, J_Jx3, S_Jx3, Nxy, Npr);

Jy1 = sparse(I_Jy1, J_Jy1, S_Jy1, Nxy, Nxy);
Jy2 = sparse(I_Jy2, J_Jy2, S_Jy2, Nxy, Nxy);
Jy3 = sparse(I_Jy3, J_Jy3, S_Jy3, Nxy, Npr);

Jc1 = sparse(I_Jc1, J_Jc1, S_Jc1, Npr, Nxy);
Jc2 = sparse(I_Jc2, J_Jc2, S_Jc2, Npr, Nxy);
% Jc3 remains zero

%% 6) Impose Dirichlet BCs on velocity (vectorized)
bV = aBC(:);  % ensure column

% Jx1
Jx1(bV,:) = 0;
Jx2(bV,:) = 0;
Jx3(bV,:) = 0;
Jx1(sub2ind(size(Jx1), bV, bV)) = 1;

% Jx2(bV,:) = 0;
% Jx2(sub2ind(size(Jx2), bV, bV)) = 1;

% Jx3(bV,:) = 0;  % entire row => 0

% Jy1
Jy1(bV,:) = 0;
Jy2(bV,:) = 0;
Jy3(bV,:) = 0;
% Jy1(sub2ind(size(Jy1), bV, bV)) = 1;

% Jy2(bV,:) = 0;
Jy2(sub2ind(size(Jy2), bV, bV)) = 1;
% 
% Jy3(bV,:) = 0; 

%% 7) Combine sub-blocks
K = [ Jx1, Jx2, Jx3;
      Jy1, Jy2, Jy3;
      Jc1, Jc2, Jc3 ];

%% 8) Fix corner pressure DOF
cornerIndex = corner;  
cornerGlobal = 2*Nxy + cornerIndex;
K(cornerGlobal, :) = 0;
K(cornerGlobal, cornerGlobal) = 1;

end % end main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function for shape derivatives at a single Gauss point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dNxAll, dNyAll, detJ] = q2ShapeDerivatives_AllNodes( ...
    xcoords, ycoords, ...
    quad1_vals, quad2_vals, ...
    dquad1_vals, dquad2_vals )

    % 1) partial derivatives wrt xi,eta
    dNdxi  = dquad1_vals .* quad2_vals;
    dNdeta = quad1_vals   .* dquad2_vals;

    % 2) form the Jacobian
    dX_dxi   = sum( xcoords .* dNdxi );
    dX_deta  = sum( xcoords .* dNdeta );
    dY_dxi   = sum( ycoords .* dNdxi );
    dY_deta  = sum( ycoords .* dNdeta );

    J = [ dX_dxi, dY_dxi ;
          dX_deta, dY_deta ];

    detJ = (dX_dxi*dY_deta - dX_deta*dY_dxi);
    invJ = inv(J);

    J11 = invJ(1,1);  J12 = invJ(1,2);
    J21 = invJ(2,1);  J22 = invJ(2,2);

    % 3) for each node i, build dNxAll(i), dNyAll(i)
    dNxAll = J11*dNdxi + J12*dNdeta;
    dNyAll = J21*dNdxi + J22*dNdeta;
end