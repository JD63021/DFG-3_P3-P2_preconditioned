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
