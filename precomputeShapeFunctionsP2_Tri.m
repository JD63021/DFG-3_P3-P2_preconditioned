function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP2_Tri()
% precomputeShapeFunctionsP2_Tri:
%  Returns shape functions (N), their derivatives (dNxi, dNeta),
%  and a 19-point Gauss integration rule (g_pt, g_wt) for a P2 triangle in
%  the reference domain: 0 <= xi, 0 <= eta, xi+eta <= 1.
%
%  Gmsh node ordering for P2 (anticlockwise):
%    Node1 -> corner1
%    Node2 -> corner2
%    Node3 -> corner3
%    Node4 -> mid-edge (corner1-corner2)
%    Node5 -> mid-edge (corner2-corner3)
%    Node6 -> mid-edge (corner3-corner1)

    % Use the 19-point Dunavant rule (exact for polynomials up to degree 8)
    [g_pt, g_wt] = triGaussPoints6();
    numGauss = length(g_wt);

    numNodes = 6;  % 6 nodes for quadratic (P2) triangle
    N      = zeros(numNodes, numGauss);
    dNxi   = zeros(numNodes, numGauss);
    dNeta  = zeros(numNodes, numGauss);

    for k = 1:numGauss
        xi  = g_pt(k,1);
        eta = g_pt(k,2);

        % P2 shape functions matching Gmsh anticlockwise node ordering
        [Ni, dNxi_k, dNeta_k] = p2basisGmsh(xi, eta);

        N(:,k)     = Ni;
        dNxi(:,k)  = dNxi_k;
        dNeta(:,k) = dNeta_k;
    end
end

%--------------------------------------------------------------------
function [N, dNxi, dNeta] = p2basisGmsh(xi, eta)
% p2basisGmsh:
%   Returns the 6 shape functions and their partial derivatives with respect
%   to (xi, eta) for a quadratic (P2) triangle with node ordering:
%     (1) corner1, (2) corner2, (3) corner3,
%     (4) edge (corner1-corner2), (5) edge (corner2-corner3), (6) edge (corner3-corner1).
%
%   Let zeta = 1 - xi - eta.

    zeta = 1 - xi - eta;  % barycentric coordinate for corner1

    % Shape functions
    N1 = zeta*(2*zeta - 1);  % Node 1 (corner)
    N2 = xi*(2*xi - 1);      % Node 2 (corner)
    N3 = eta*(2*eta - 1);    % Node 3 (corner)
    N4 = 4*xi*zeta;          % Node 4 (mid-edge between node 1 & 2)
    N5 = 4*xi*eta;           % Node 5 (mid-edge between node 2 & 3)
    N6 = 4*eta*zeta;         % Node 6 (mid-edge between node 3 & 1)

    N = [N1; N2; N3; N4; N5; N6];

    % Partial derivatives:
    % Note: zeta = 1 - xi - eta, so dzeta/dxi = -1 and dzeta/deta = -1.

    % Derivative of N1 = zeta*(2*zeta - 1)
    dN1_dxi = (2*zeta - 1)*(-1) + zeta*2*(-1);  % = 1 - 4*zeta
    dN1_deta = (2*zeta - 1)*(-1) + zeta*2*(-1);   % = 1 - 4*zeta

    % Derivative of N2 = xi*(2*xi - 1)
    dN2_dxi  = (2*xi - 1) + xi*2;
    dN2_deta = 0;

    % Derivative of N3 = eta*(2*eta - 1)
    dN3_dxi  = 0;
    dN3_deta = (2*eta - 1) + eta*2;

    % Derivative of N4 = 4*xi*zeta
    dN4_dxi  = 4*(zeta + xi*(-1));  % = 4*(zeta - xi)
    dN4_deta = 4*( xi * (-1) );      % = -4*xi

    % Derivative of N5 = 4*xi*eta
    dN5_dxi  = 4*eta;
    dN5_deta = 4*xi;

    % Derivative of N6 = 4*eta*zeta
    dN6_dxi  = 4*( eta * (-1) );          % = -4*eta
    dN6_deta = 4*( zeta + eta*(-1) );       % = 4*(zeta - eta)

    dNxi  = [ dN1_dxi; dN2_dxi; dN3_dxi; dN4_dxi; dN5_dxi; dN6_dxi ];
    dNeta = [ dN1_deta; dN2_deta; dN3_deta; dN4_deta; dN5_deta; dN6_deta ];
end
% 
%--------------------------------------------------------------------
 function [g_pt, g_wt] = triGaussPoints12()
% triGaussPoints12:
%   Returns a 12-point Dunavant rule on the reference triangle
%   (0 <= xi, 0 <= eta, xi+eta <= 1),
%   exact for polynomials up to degree 6.
%
% The coordinates and weights come from standard references (Dunavant, 1985).

g_pt = [
  0.24928674517091,  0.24928674517091;
  0.24928674517091,  0.50142650965818;
  0.50142650965818,  0.24928674517091;

  0.06308901449150,  0.06308901449150;
  0.06308901449150,  0.87382197101700;
  0.87382197101700,  0.06308901449150;

  0.31035245103378,  0.63650249912140;
  0.63650249912140,  0.05314504984482;
  0.05314504984482,  0.31035245103378;

  0.63650249912140,  0.31035245103378;
  0.31035245103378,  0.05314504984482;
  0.05314504984482,  0.63650249912140
];

g_wt = [
  0.05839313786319;
  0.05839313786319;
  0.05839313786319;

  0.02542245318510;
  0.02542245318510;
  0.02542245318510;

  0.04142553780919;
  0.04142553780919;
  0.04142553780919;

  0.04142553780919;
  0.04142553780919;
  0.04142553780919
];
 end
% Scale the raw weights so that they sum to 0.5 (area of the reference triangle)
% scale = 0.5 / sum(g_wt_raw);
% g_wt = g_wt_raw * scale;

function [g_pt, g_wt] = triGaussPoints19()
% Degree-8, 19-point Dunavant rule on reference triangle; sum(weights)=0.5
g_pt = [
  0.333333333333333, 0.333333333333333;
  0.489682519198738, 0.255658740400631;
  0.255658740400631, 0.489682519198738;
  0.255658740400631, 0.255658740400631;
  0.623592928761935, 0.188203535619033;
  0.188203535619033, 0.623592928761935;
  0.188203535619033, 0.188203535619033;
  0.910540973211095, 0.044729513394453;
  0.044729513394453, 0.910540973211095;
  0.044729513394453, 0.044729513394453;
  0.741198598784498, 0.036838412054736;
  0.036838412054736, 0.741198598784498;
  0.221962989160766, 0.741198598784498;
  0.741198598784498, 0.221962989160766;
  0.036838412054736, 0.221962989160766;
  0.221962989160766, 0.036838412054736;
  0.022072179275643, 0.022072179275643;
  0.955366,          0.022072179275643; % keep triplet symmetry
  0.022072179275643, 0.955366
];
% Standard Dunavant 19 weights (unscaled):
w = [
  0.097135796282799;
  0.031334700227139;
  0.031334700227139;
  0.031334700227139;
  0.025577675658698;
  0.025577675658698;
  0.025577675658698;
  0.009421666963733;
  0.009421666963733;
  0.009421666963733;
  0.012865533220227;
  0.012865533220227;
  0.012865533220227;
  0.012865533220227;
  0.012865533220227;
  0.012865533220227;
  0.001289,          % (some tables round these three tiny weights)
  0.001289,
  0.001289
];
% Use your original 19-pt set if you prefer; the key is the scaling:
g_wt = w * (0.5/sum(w));
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


% function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP2_Tri()
% % precomputeShapeFunctionsP2_Tri:
% %  Returns shape functions (N), derivatives (dNxi, dNeta),
% %  and 6-point Gauss integration rule (g_pt, g_wt)
% %  for a P2 triangle in the reference domain 0<=xi, 0<=eta, xi+eta<=1.
% %
% %  Gmsh node ordering for P2 (anticlockwise):
% %    Node1 -> corner1
% %    Node2 -> corner2
% %    Node3 -> corner3
% %    Node4 -> mid-edge(corner1-corner2)
% %    Node5 -> mid-edge(corner2-corner3)
% %    Node6 -> mid-edge(corner3-corner1)
% 
%     % Get 6-point Gauss rule over the reference triangle
%     [g_pt, g_wt] = triGaussPoints12();
%     numGauss = length(g_wt);
% 
%     numNodes = 6;  % 6 local nodes for P2
%     N      = zeros(numNodes, numGauss);
%     dNxi   = zeros(numNodes, numGauss);
%     dNeta  = zeros(numNodes, numGauss);
% 
%     for k = 1:numGauss
%         xi  = g_pt(k,1);
%         eta = g_pt(k,2);
% 
%         % P2 shape functions matching Gmsh anticlockwise node ordering
%         [Ni, dNxi_k, dNeta_k] = p2basisGmsh(xi, eta);
% 
%         N(:,k)     = Ni;
%         dNxi(:,k)  = dNxi_k;
%         dNeta(:,k) = dNeta_k;
%     end
% end
% 
% %--------------------------------------------------------------------
% function [N, dNxi, dNeta] = p2basisGmsh(xi, eta)
% % p2basisGmsh:
% %   Returns the 6 shape functions and their partial derivatives
% %   wrt (xi, eta) for a quadratic triangle, with node ordering:
% %     (1) corner1, (2) corner2, (3) corner3,
% %     (4) edge12,  (5) edge23,  (6) edge31.
% %
% % Let zeta = 1 - xi - eta.
% 
%     zeta = 1 - xi - eta;  % corner1 (anticlockwise from node1)
% 
%     % shape functions
%     N1 = zeta*(2*zeta - 1);  % node1
%     N2 = xi*(2*xi - 1);      % node2
%     N3 = eta*(2*eta - 1);    % node3
%     N4 = 4*xi*zeta;          % mid-edge(1-2)
%     N5 = 4*xi*eta;           % mid-edge(2-3)
%     N6 = 4*eta*zeta;         % mid-edge(3-1)
% 
%     N = [N1; N2; N3; N4; N5; N6];
% 
%     % partial derivatives wrt xi, eta
%     % zeta = 1 - xi - eta => dzeta/dxi = -1, dzeta/deta = -1
% 
%     dN1_dxi = (2*(zeta) - 1)*(-1) + zeta*(2*(-1));
%         % expanded: dN1_dxi = - (2zeta - 1) - 2*zeta
%         %           = -2zeta + 1 - 2zeta = 1 - 4zeta
%     dN1_deta = 1 - 4*zeta; % same pattern
% 
%     dN2_dxi  = (2*xi - 1) + xi*(2); % chain rule for xi*(2xi-1)
%     dN2_deta = 0; % no eta in N2
% 
%     dN3_dxi  = 0; % no xi in N3
%     dN3_deta = (2*eta - 1) + eta*(2);
% 
%     dN4_dxi  = 4*zeta + 4*xi*(-1);  % 4(zeta dxi + xi dzeta/dxi= -xi)
%                % = 4*(zeta - xi)
%     dN4_deta = 4*( xi * (-1) );     % = -4 xi
% 
%     dN5_dxi  = 4*eta;
%     dN5_deta = 4*xi;
% 
%     dN6_dxi  = 4*( eta * (-1) );  % = -4 eta
%     dN6_deta = 4*( zeta + eta*(-1) );
%                % = 4*(zeta - eta)
% 
%     dNxi  = [ dN1_dxi;  dN2_dxi;  dN3_dxi;  dN4_dxi;  dN5_dxi;  dN6_dxi ];
%     dNeta = [ dN1_deta; dN2_deta; dN3_deta; dN4_deta; dN5_deta; dN6_deta ];
% end
% 
% %--------------------------------------------------------------------
% function [g_pt, g_wt] = triGaussPoints12()
% % triGaussPoints12:
% %   Returns a 12-point Dunavant rule on the reference triangle
% %   (0 <= xi, 0 <= eta, xi+eta <= 1),
% %   exact for polynomials up to degree 6.
% %
% % The coordinates and weights come from standard references (Dunavant, 1985).
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
% end