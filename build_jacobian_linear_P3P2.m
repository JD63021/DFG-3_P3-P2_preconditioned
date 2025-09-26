function [K, Jx1, Jx2, Jx3, Jy1, Jy2, Jy3, Jc1, Jc2] = build_jacobian_linear_P3P2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, U1, U2, U3, gamma, mu, rho, corner) %#ok<INUSD>
% Linear (u–independent) Jacobian blocks for P3–P2:
%   - viscosity + grad-div in momentum (VV blocks)
%   - pressure gradient (VP) and continuity (PV)
%   - applies velocity Dirichlet rows (identity on Jx1 & Jy2) and pressure pin

% ---------- sizes ----------
Nxy         = length(nodeInfo.velocity.x);   % # P3 vel nodes
Npr         = max(elemInfo.presElements(:)); % # P2 pres nodes
nEl         = size(elemInfo.velElements,1);
nV = 10; nP = 6;

% ---------- prealloc triplets ----------
nP3P3 = 100*nEl;  nP3P2 = 60*nEl;  nP2P3 = 60*nEl;

I_jx1 = zeros(nP3P3,1); J_jx1 = I_jx1; V_jx1 = I_jx1;
I_jx2 = zeros(nP3P3,1); J_jx2 = I_jx2; V_jx2 = I_jx2;
I_jx3 = zeros(nP3P2,1); J_jx3 = I_jx3; V_jx3 = I_jx3;

I_jy1 = zeros(nP3P3,1); J_jy1 = I_jy1; V_jy1 = I_jy1;
I_jy2 = zeros(nP3P3,1); J_jy2 = I_jy2; V_jy2 = I_jy2;
I_jy3 = zeros(nP3P2,1); J_jy3 = I_jy3; V_jy3 = I_jy3;

I_jc1 = zeros(nP2P3,1); J_jc1 = I_jc1; V_jc1 = I_jc1;
I_jc2 = zeros(nP2P3,1); J_jc2 = I_jc2; V_jc2 = I_jc2;

off_jx1 = 0; off_jx2 = 0; off_jx3 = 0;
off_jy1 = 0; off_jy2 = 0; off_jy3 = 0;
off_jc1 = 0; off_jc2 = 0;

% ---------- shapes ----------
[Nvel,dNxiVel,dNetaVel,wVel,~] = precomputeShapeFunctionsP3_Tri();
[Npre,~,~,~,~]                 = precomputeShapeFunctionsP2_Tri();
numGp = length(wVel);

% ---------- element loop ----------
for e = 1:nEl
    Kvel = elemInfo.velElements(e,:);   % 10
    Kpr  = elemInfo.presElements(e,:);  %  6

    xV = nodeInfo.velocity.x(Kvel);  yV = nodeInfo.velocity.y(Kvel);

    Kx1e=zeros(nV,nV); Kx2e=zeros(nV,nV); Kx3e=zeros(nV,nP);
    Ky1e=zeros(nV,nV); Ky2e=zeros(nV,nV); Ky3e=zeros(nV,nP);
    Kc1e=zeros(nP,nV); Kc2e=zeros(nP,nV);

    % precompute per-GP geometry
    dNxAll = zeros(nV,numGp); dNyAll = zeros(nV,numGp);
    NAll   = zeros(nV,numGp); detJg  = zeros(numGp,1);
    NpAll  = zeros(nP,numGp);

    for gp=1:numGp
        NAll(:,gp) = Nvel(:,gp);
        [dNx,dNy,detJ] = p3ShapeDerivativesAllNodes(xV,yV,dNxiVel(:,gp),dNetaVel(:,gp));
        dNxAll(:,gp)=dNx; dNyAll(:,gp)=dNy; detJg(gp)=detJ;
        if gp<=size(Npre,2), NpAll(:,gp)=Npre(:,gp); end
    end

    % ---- VV blocks: viscosity + grad-div only (NO convection here) ----
    for i=1:nV
        for k=1:nV
            jxm_x=0; jxm_y=0; jym_x=0; jym_y=0;
            for gp=1:numGp
                wt=wVel(gp); detJ=detJg(gp); cf=wt*detJ;
                dNxi=dNxAll(i,gp); dNyi=dNyAll(i,gp);
                dNxk=dNxAll(k,gp); dNyk=dNyAll(k,gp);

                jxm_x = jxm_x + ( mu*(2*dNxi*dNxk + dNyi*dNyk) + 0        )*cf;
                jxm_y = jxm_y + ( mu*(dNyi*dNxk)                 + 0        )*cf;
                jym_y = jym_y + ( mu*(dNxi*dNxk + 2*dNyi*dNyk) + 0        )*cf;
                jym_x = jym_x + ( mu*(dNxi*dNyk)                 + 0        )*cf;
            end
            Kx1e(i,k)=jxm_x; Kx2e(i,k)=jxm_y;
            Ky1e(i,k)=jym_x; Ky2e(i,k)=jym_y;
        end
    end

    % ---- VP: pressure gradient ----
    for i=1:nV
        for p=1:nP
            jx=0; jy=0;
            for gp=1:numGp
                wt=wVel(gp); detJ=detJg(gp); cf=wt*detJ;
                dNxi=dNxAll(i,gp); dNyi=dNyAll(i,gp);
                Np  =NpAll(p,gp);
                jx = jx - dNxi*Np*cf;
                jy = jy - dNyi*Np*cf;
            end
            Kx3e(i,p)=jx; Ky3e(i,p)=jy;
        end
    end

    % ---- PV: continuity ----
    for p=1:nP
        for k=1:nV
            jcx=0; jcy=0;
            for gp=1:numGp
                wt=wVel(gp); detJ=detJg(gp); cf=wt*detJ;
                Np  =NpAll(p,gp);
                dNxk=dNxAll(k,gp); dNyk=dNyAll(k,gp);
                jcx = jcx - dNxk*Np*cf;
                jcy = jcy - dNyk*Np*cf;
            end
            Kc1e(p,k)=jcx; Kc2e(p,k)=jcy;
        end
    end

    % ---- scatter to triplets ----
    [rv,cv] = ndgrid(Kvel,Kvel); rng = 1:(nV*nV);
    ix = off_jx1 + rng; I_jx1(ix)=rv(:); J_jx1(ix)=cv(:); V_jx1(ix)=Kx1e(:); off_jx1=ix(end);
    ix = off_jx2 + rng; I_jx2(ix)=rv(:); J_jx2(ix)=cv(:); V_jx2(ix)=Kx2e(:); off_jx2=ix(end);
    ix = off_jy1 + rng; I_jy1(ix)=rv(:); J_jy1(ix)=cv(:); V_jy1(ix)=Ky1e(:); off_jy1=ix(end);
    ix = off_jy2 + rng; I_jy2(ix)=rv(:); J_jy2(ix)=cv(:); V_jy2(ix)=Ky2e(:); off_jy2=ix(end);

    [rvp,cvp] = ndgrid(Kvel,Kpr); rng = 1:(nV*nP);
    ix = off_jx3 + rng; I_jx3(ix)=rvp(:); J_jx3(ix)=cvp(:); V_jx3(ix)=Kx3e(:); off_jx3=ix(end);
    ix = off_jy3 + rng; I_jy3(ix)=rvp(:); J_jy3(ix)=cvp(:); V_jy3(ix)=Ky3e(:); off_jy3=ix(end);

    [rpv,cpv] = ndgrid(Kpr,Kvel); rng = 1:(nP*nV);
    ix = off_jc1 + rng; I_jc1(ix)=rpv(:); J_jc1(ix)=cpv(:); V_jc1(ix)=Kc1e(:); off_jc1=ix(end);
    ix = off_jc2 + rng; I_jc2(ix)=rpv(:); J_jc2(ix)=cpv(:); V_jc2(ix)=Kc2e(:); off_jc2=ix(end);
end

% ---------- assemble blocks ----------
Jx1 = sparse(I_jx1,J_jx1,V_jx1,Nxy,Nxy);
Jx2 = sparse(I_jx2,J_jx2,V_jx2,Nxy,Nxy);
Jx3 = sparse(I_jx3,J_jx3,V_jx3,Nxy,Npr);

Jy1 = sparse(I_jy1,J_jy1,V_jy1,Nxy,Nxy);
Jy2 = sparse(I_jy2,J_jy2,V_jy2,Nxy,Nxy);
Jy3 = sparse(I_jy3,J_jy3,V_jy3,Nxy,Npr);

Jc1 = sparse(I_jc1,J_jc1,V_jc1,Npr,Nxy);
Jc2 = sparse(I_jc2,J_jc2,V_jc2,Npr,Nxy);

% pressure block for pin
Jc3 = spalloc(Npr,Npr,1);

% ---------- apply velocity BC rows (identity on Jx1 & Jy2) ----------
b = boundaryInfo.allVelNodes(:);
if ~isempty(b)
    % x-momentum
    Jx1(b,:) = 0; Jx2(b,:) = 0; Jx3(b,:) = 0;
    Jx1(sub2ind(size(Jx1),b,b)) = 1;

    % y-momentum
    Jy1(b,:) = 0; Jy2(b,:) = 0; Jy3(b,:) = 0;
    Jy2(sub2ind(size(Jy2),b,b)) = 1;
end

% ---------- pressure pin (row-zero + diag=1) ----------
cornerIndex = corner;
Jc3(cornerIndex,:) = 0;
Jc3(cornerIndex,cornerIndex) = 1;

% ---------- global K ----------
K = [ Jx1, Jx2, Jx3;
      Jy1, Jy2, Jy3;
      Jc1, Jc2, Jc3 ];
end

% ===== helpers =====
function [dNxAll,dNyAll,detJ] = p3ShapeDerivativesAllNodes(xcoords,ycoords,dNxi,dNeta)
dX_dxi  = sum(xcoords.*dNxi);
dX_deta = sum(xcoords.*dNeta);
dY_dxi  = sum(ycoords.*dNxi);
dY_deta = sum(ycoords.*dNeta);
J = [dX_dxi, dY_dxi; dX_deta, dY_deta];
detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
invJ = inv(J);
dNxAll = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
dNyAll = invJ(2,1)*dNxi + invJ(2,2)*dNeta;
end
