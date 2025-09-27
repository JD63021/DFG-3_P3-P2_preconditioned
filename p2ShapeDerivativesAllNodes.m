
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated helper function for shape derivatives at a single Gauss point
% in a P2 triangular element. The same formula applies in principle:
%   x = sum_i xcoords(i) * N_i
% with derivatives wrt (xi,eta).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dNxAll, dNyAll, detJ] = p2ShapeDerivativesAllNodes( ...
    xcoords, ycoords, dNxi, dNeta)
% xcoords, ycoords: the physical coords [x1..x6], [y1..y6] for the 6 local nodes
% dNxi, dNeta: partial derivatives of the shape functions wrt the reference coords
%   each is a 6×1 vector (one derivative for each local node).
%
% Output:
%   dNxAll, dNyAll: the derivatives of shape wrt physical x,y (6×1)
%   detJ: the determinant of the Jacobian

    % Jacobian from reference (xi,eta) to physical (x,y):
    dX_dxi  = sum(xcoords .* dNxi);
    dX_deta = sum(xcoords .* dNeta);
    dY_dxi  = sum(ycoords .* dNxi);
    dY_deta = sum(ycoords .* dNeta);

    J = [ dX_dxi,  dY_dxi ;
          dX_deta, dY_deta ];

    detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
    invJ = inv(J);

    % Transform local derivatives into physical derivatives
    % for each shape function i:
    %   dN/dx = invJ(1,1)*dN/dxi + invJ(1,2)*dN/deta
    %   dN/dy = invJ(2,1)*dN/dxi + invJ(2,2)*dN/deta
    dNxAll = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
    dNyAll = invJ(2,1)*dNxi + invJ(2,2)*dNeta;
end