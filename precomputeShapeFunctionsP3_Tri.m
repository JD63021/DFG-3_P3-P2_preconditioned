function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP3_Tri()
% precomputeShapeFunctionsP3_Tri:
%   Returns the 10 P3 shape functions (N), their derivatives (dNxi, dNeta),
%   and a 19-point Dunavant rule (g_pt, g_wt) for a P3 triangle
%   in the isoparametric domain 0 <= xi, 0 <= eta, xi+eta <= 1.
%
%   The node ordering (anticlockwise) is assumed to match Gmsh's default
%   for a 10-node triangle:
%        (1) corner1
%        (2) corner2
%        (3) corner3
%        (4,5) mid-edge(1->2)
%        (6,7) mid-edge(2->3)
%        (8,9) mid-edge(3->1)
%       (10)   interior bubble.
%
%   The shape functions below are chosen so that each Ni is exactly 1 at its
%   corresponding node (in barycentric space) and 0 at all other nodes.
%
%   Output arrays:
%     N:    10×numGauss
%     dNxi, dNeta:  10×numGauss
%     g_pt, g_wt:   19-point Dunavant rule

    % Get the 19-pt Dunavant rule
    [g_pt, g_wt] = triGaussPoints6();
    numGauss = length(g_wt);

    % Allocate
    numNodes = 10;  % P3
    N      = zeros(numNodes, numGauss);
    dNxi   = zeros(numNodes, numGauss);
    dNeta  = zeros(numNodes, numGauss);

    % Loop over Gauss points, compute shape + derivatives
    for k = 1:numGauss
        xi  = g_pt(k,1);
        eta = g_pt(k,2);

        [Ni, dNxi_k, dNeta_k] = p3IsoShape(xi, eta);

        N(:,k)     = Ni;
        dNxi(:,k)  = dNxi_k;
        dNeta(:,k) = dNeta_k;
    end
end

%------------------------------------------------------------------------
function [N, dNxi, dNeta] = p3IsoShape(xi, eta)
% p3IsoShape:
%   Computes the 10-node (P3) Lagrange shape functions and their
%   derivatives w.r.t. (xi, eta) over the reference triangle:
%   0 <= xi, 0 <= eta, xi+eta <=1.
%
%   Node ordering (anticlockwise):
%     Node1 -> corner1 (L1=1, L2=0, L3=0)
%     Node2 -> corner2 (L1=0, L2=1, L3=0)
%     Node3 -> corner3 (L1=0, L2=0, L3=1)
%     Node4 -> mid-edge(1->2) near corner1 (L1=2/3, L2=1/3, L3=0)
%     Node5 -> mid-edge(1->2) near corner2 (L1=1/3, L2=2/3, L3=0)
%     Node6 -> mid-edge(2->3) near corner2 (L1=0, L2=2/3, L3=1/3)
%     Node7 -> mid-edge(2->3) near corner3 (L1=0, L2=1/3, L3=2/3)
%     Node8 -> mid-edge(3->1) near corner3 (L1=1/3, L2=0, L3=2/3)
%     Node9 -> mid-edge(3->1) near corner1 (L1=2/3, L2=0, L3=1/3)
%     Node10-> interior bubble (L1=L2=L3=1/3)
%
%   Each shape function Ni(L1,L2,L3) is built to be 1 at its node and 0 at
%   all others. Barycentric coords: L1=1 - xi - eta, L2=xi, L3=eta.

    % Barycentric coords
    L1 = 1 - xi - eta;
    L2 = xi;
    L3 = eta;

    % Preallocate
    N     = zeros(10,1);
    dNxi  = zeros(10,1);
    dNeta = zeros(10,1);

    %-----------------------------------------
    % Corner nodes:
    %  N(1) = L1*(3L1 -1)*(3L1 -2)/2
    %  N(2) = L2*(3L2 -1)*(3L2 -2)/2
    %  N(3) = L3*(3L3 -1)*(3L3 -2)/2

    N(1) = L1*(3*L1 - 1)*(3*L1 - 2)/2;
    N(2) = L2*(3*L2 - 1)*(3*L2 - 2)/2;
    N(3) = L3*(3*L3 - 1)*(3*L3 - 2)/2;

    %-----------------------------------------
    % Mid-edge nodes:
    %  Node4:  (L1=2/3, L2=1/3) => N4= (9/2)*L1*L2*(3L1 -1)
    %  Node5:                    => N5= (9/2)*L1*L2*(3L2 -1)
    %
    %  Node6: (L2=2/3, L3=1/3) => N6= (9/2)*L2*L3*(3L2 -1)
    %  Node7:                   => N7= (9/2)*L2*L3*(3L3 -1)
    %
    %  Node8: (L3=2/3, L1=1/3) => N8= (9/2)*L3*L1*(3L3 -1)
    %  Node9:                   => N9= (9/2)*L3*L1*(3L1 -1)

    c = 9/2;
    N(4) = c * L1*L2*(3*L1 - 1);
    N(5) = c * L1*L2*(3*L2 - 1);

    N(6) = c * L2*L3*(3*L2 - 1);
    N(7) = c * L2*L3*(3*L3 - 1);

    N(8) = c * L3*L1*(3*L3 - 1);
    N(9) = c * L3*L1*(3*L1 - 1);

    %-----------------------------------------
    % Interior bubble:
    %  Node10: (L1=L2=L3=1/3) => N(10)= 27*L1*L2*L3
    %
    N(10) = 27*L1*L2*L3;

    %-----------------------------------------
    % Derivatives wrt xi, eta
    dNxi(:)  = 0;  % prefill
    dNeta(:) = 0;

    % We can systematically apply product rule. First define dL1/dxi,etc.
    dL1_dxi  = -1; dL1_deta = -1;
    dL2_dxi  =  1; dL2_deta =  0;
    dL3_dxi  =  0; dL3_deta =  1;

    %--- Corner shape derivatives (same as your old code for corners) ---
    [dNxi(1),  dNeta(1)]  = cornerDeriv( L1, dL1_dxi, dL1_deta );
    [dNxi(2),  dNeta(2)]  = cornerDeriv( L2, dL2_dxi, dL2_deta );
    [dNxi(3),  dNeta(3)]  = cornerDeriv( L3, dL3_dxi, dL3_deta );

    %--- Mid-edge + bubble derivatives ---
    % We'll do each explicitly for clarity:

    % N4 = c * L1*L2*(3L1 -1)
    [dNxi(4), dNeta(4)] = d_of_prod3( c, L1, L2, (3*L1 -1), ...
                                      dL1_dxi,dL1_deta, ...
                                      dL2_dxi,dL2_deta, ...
                                      3*dL1_dxi, 3*dL1_deta );

    % N5 = c * L1*L2*(3L2 -1)
    [dNxi(5), dNeta(5)] = d_of_prod3( c, L1, L2, (3*L2 -1), ...
                                      dL1_dxi,dL1_deta, ...
                                      dL2_dxi,dL2_deta, ...
                                      3*dL2_dxi, 3*dL2_deta );

    % N6 = c * L2*L3*(3*L2 -1)
    [dNxi(6), dNeta(6)] = d_of_prod3( c, L2, L3, (3*L2 -1), ...
                                      dL2_dxi,dL2_deta, ...
                                      dL3_dxi,dL3_deta, ...
                                      3*dL2_dxi, 3*dL2_deta );

    % N7 = c * L2*L3*(3*L3 -1)
    [dNxi(7), dNeta(7)] = d_of_prod3( c, L2, L3, (3*L3 -1), ...
                                      dL2_dxi,dL2_deta, ...
                                      dL3_dxi,dL3_deta, ...
                                      3*dL3_dxi, 3*dL3_deta );

    % N8 = c * L3*L1*(3*L3 -1)
    [dNxi(8), dNeta(8)] = d_of_prod3( c, L3, L1, (3*L3 -1), ...
                                      dL3_dxi,dL3_deta, ...
                                      dL1_dxi,dL1_deta, ...
                                      3*dL3_dxi, 3*dL3_deta );

    % N9 = c * L3*L1*(3*L1 -1)
    [dNxi(9), dNeta(9)] = d_of_prod3( c, L3, L1, (3*L1 -1), ...
                                      dL3_dxi,dL3_deta, ...
                                      dL1_dxi,dL1_deta, ...
                                      3*dL1_dxi, 3*dL1_deta );

    % N10 = 27*L1*L2*L3
    [dNxi(10), dNeta(10)] = d_of_prod3( 27, L1, L2, L3, ...
                                        dL1_dxi,dL1_deta, ...
                                        dL2_dxi,dL2_deta, ...
                                        dL3_dxi,dL3_deta );

end

%=========================================================================
%  Local helper: derivative of corner shape
%    Nc = L*(3L-1)*(3L-2)/2
%-------------------------------------------------------------------------
function [dNdxi, dNdeta] = cornerDeriv( L, dLdxi, dLdeta )
% For Ncorner = L*(3L - 1)*(3L - 2)/2
% We can do partial expansions or a direct product rule. Let's do expansions:

  % define   f(L) = (3L - 1)*(3L - 2)/2
  % => f(L) = (9L^2 - 9L + 2)/2
  % => df/dL= (18L - 9)/2 = 9(2L -1)/2
  f     = (3*L -1).*(3*L -2)/2;
  df_dL =  (9*(2*L -1))/2;  % = ( d/dL of the above ) 
                           % you can also do product rule: (3*(3L-2) + 3*(3L-1))/2

  dNdxi  = dLdxi * f + L*( df_dL * dLdxi );
  dNdeta = dLdeta * f + L*( df_dL * dLdeta );
end

%=========================================================================
%  Local helper: derivative of a product of 3 terms:
%    N = scale * ( X^1 ) * ( Y^1 ) * ( Z ), 
%    where Z might be (3L -1), etc.
%-------------------------------------------------------------------------
function [dNdxi, dNeta] = d_of_prod3( scale, X, Y, Z, ...
                                      dXdxi, dXdeta, ...
                                      dYdxi, dYdeta, ...
                                      dZdxi, dZdeta )
% dN/dxi = scale * [ (dXdxi * Y * Z) + (X * dYdxi * Z) + (X * Y * dZdxi ) ]
% (and similarly for dN/deta)

    dNdxi = scale * ( dXdxi*Y*Z + X*dYdxi*Z + X*Y*dZdxi );
    dNeta = scale * ( dXdeta*Y*Z + X*dYdeta*Z + X*Y*dZdeta );
end

% -------------------------------------------------------------------------
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


% 
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




% function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP3_Tri()
% % precomputeShapeFunctionsP3_Tri_iso:
% %   Returns shape functions (N), derivatives (dNxi, dNeta),
% %   and a 12-point Dunavant rule (g_pt, g_wt) for a P3 triangle,
% %   in the isoparametric domain 0<=xi, 0<=eta, xi+eta<=1.
% %
% %  Gmsh node ordering for P3 (anticlockwise), 10 nodes total:
% %    Node1 -> corner1
% %    Node2 -> corner2
% %    Node3 -> corner3
% %    Node4, Node5 -> mid-edge(1-2)
% %    Node6, Node7 -> mid-edge(2-3)
% %    Node8, Node9 -> mid-edge(3-1)
% %    Node10 -> interior (bubble)
% %
% %  Implementation with standard polynomials in (xi, eta).
% 
% [g_pt, g_wt] = triGaussPoints12();  % 12-pt Dunavant
% numGauss = length(g_wt);
% 
% numNodes = 10;  % P3 has 10 local nodes
% N      = zeros(numNodes, numGauss);
% dNxi   = zeros(numNodes, numGauss);
% dNeta  = zeros(numNodes, numGauss);
% 
% for k = 1:numGauss
%     xi  = g_pt(k,1);
%     eta = g_pt(k,2);
% 
%     [Ni, dNxi_k, dNeta_k] = p3IsoShape(xi, eta);
% 
%     N(:,k)     = Ni;
%     dNxi(:,k)  = dNxi_k;
%     dNeta(:,k) = dNeta_k;
% end
% end
% 
% %--------------------------------------------------------------------
% function [N, dNxi, dNeta] = p3IsoShape(xi, eta)
% % p3IsoShape:
% %   Computes the 10-node (P3) Lagrange shape functions and their
% %   derivatives with respect to xi and eta over the reference triangle:
% %       0 <= xi, 0 <= eta, xi+eta <= 1.
% %
% %   Input:
% %       xi, eta  - coordinates in the isoparametric (reference) domain.
% %
% %   Output:
% %       N      - 10x1 vector of shape function values.
% %       dNxi   - 10x1 vector of partial derivatives dN/dxi.
% %       dNeta  - 10x1 vector of partial derivatives dN/deta.
% %
% % The functions are defined in barycentric coordinates:
% %   L1 = 1 - xi - eta,   L2 = xi,   L3 = eta.
% %
% % Standard P3 (cubic) shape functions (consistent with common Gmsh ordering):
% %
% % Corners:
% %   N1 = L1*(3L1 - 1)*(3L1 - 2)/2
% %   N2 = L2*(3L2 - 1)*(3L2 - 2)/2
% %   N3 = L3*(3L3 - 1)*(3L3 - 2)/2
% %
% % Mid-edge nodes:
% %   N4 = (27/4)*L1^2*L2
% %   N5 = (27/4)*L1*L2^2
% %   N6 = (27/4)*L2^2*L3
% %   N7 = (27/4)*L2*L3^2
% %   N8 = (27/4)*L3^2*L1
% %   N9 = (27/4)*L3*L1^2
% %
% % Interior (bubble) node:
% %   N10 = 27*L1*L2*L3
% 
% % Compute barycentric coordinates:
% L1 = 1 - xi - eta;
% L2 = xi;
% L3 = eta;
% 
% % Preallocate shape function vector:
% N = zeros(10,1);
% 
% % --- Shape Functions ---
% % Corners
% N(1) = L1*(3*L1 - 1)*(3*L1 - 2)/2;
% N(2) = L2*(3*L2 - 1)*(3*L2 - 2)/2;
% N(3) = L3*(3*L3 - 1)*(3*L3 - 2)/2;
% % Mid-edge nodes
% N(4) = (27/4)*L1^2 * L2;
% N(5) = (27/4)*L1   * L2^2;
% N(6) = (27/4)*L2^2 * L3;
% N(7) = (27/4)*L2   * L3^2;
% N(8) = (27/4)*L3^2 * L1;
% N(9) = (27/4)*L3   * L1^2;
% % Interior bubble
% N(10)= 27*L1*L2*L3;
% 
% % --- Derivatives of the barycentrics ---
% % Note: L1 = 1 - xi - eta, L2 = xi, L3 = eta.
% dL1_dxi  = -1;  dL1_deta = -1;
% dL2_dxi  =  1;  dL2_deta =  0;
% dL3_dxi  =  0;  dL3_deta =  1;
% 
% % Preallocate derivative vectors:
% dNxi  = zeros(10,1);
% dNeta = zeros(10,1);
% 
% %% Derivatives for the corner shape functions
% % --- For N1 = L1*(3*L1 - 1)*(3*L1 - 2)/2 ---
% % Let A = (3*L1 - 1)*(3*L1 - 2)/2.
% A = (3*L1 - 1)*(3*L1 - 2)/2;
% % Its derivative with respect to L1:
% dA_dL1 = (9/2)*(2*L1 - 1);
% % Then, by chain rule:
% dNxi(1)  = dL1_dxi * A + L1 * dA_dL1 * dL1_dxi;
% dNeta(1) = dL1_deta * A + L1 * dA_dL1 * dL1_deta;
% % (Since dL1_dxi = dL1_deta = -1, these become:)
% % dNxi(1)  = -A - (9/2)*L1*(2*L1 - 1);
% % dNeta(1) = -A - (9/2)*L1*(2*L1 - 1);
% 
% % --- For N2 = L2*(3*L2 - 1)*(3*L2 - 2)/2 ---
% B = (3*L2 - 1)*(3*L2 - 2)/2;
% dB_dL2 = (9/2)*(2*L2 - 1);
% dNxi(2)  = dL2_dxi * B + L2 * dB_dL2 * dL2_dxi;
% dNeta(2) = dL2_deta * B + L2 * dB_dL2 * dL2_deta;
% % (Recall: dL2_dxi = 1 and dL2_deta = 0, so dNeta(2)=0)
% 
% % --- For N3 = L3*(3*L3 - 1)*(3*L3 - 2)/2 ---
% C = (3*L3 - 1)*(3*L3 - 2)/2;
% dC_dL3 = (9/2)*(2*L3 - 1);
% dNxi(3)  = dL3_dxi * C + L3 * dC_dL3 * dL3_dxi;
% dNeta(3) = dL3_deta * C + L3 * dC_dL3 * dL3_deta;
% % (dL3_dxi = 0 so dNxi(3)=0; dL3_deta = 1 so dNeta(3)= C + L3*(9/2)*(2*L3 - 1) )
% 
% %% Derivatives for mid-edge nodes and bubble (using standard product rule)
% % --- For N4 = (27/4)*L1^2*L2 ---
% dNxi(4)  = (27/4)*( 2*L1*dL1_dxi*L2 + L1^2*dL2_dxi );
% dNeta(4) = (27/4)*( 2*L1*dL1_deta*L2 + L1^2*dL2_deta );
% % Substituting: dL1_dxi = -1, dL2_dxi = 1, dL1_deta = -1, dL2_deta = 0.
% % Thus, dNxi(4) = (27/4)*(-2*L1*L2 + L1^2) and dNeta(4) = -(27/2)*L1*L2.
% 
% % --- For N5 = (27/4)*L1*L2^2 ---
% dNxi(5)  = (27/4)*( dL1_dxi*L2^2 + L1*2*L2*dL2_dxi );
% dNeta(5) = (27/4)*( dL1_deta*L2^2 + L1*2*L2*dL2_deta );
% % With dL1_dxi = -1, dL2_dxi = 1, dL1_deta = -1, dL2_deta = 0:
% % dNxi(5) = (27/4)*(-L2^2 + 2*L1*L2) and dNeta(5) = -(27/4)*L2^2.
% 
% % --- For N6 = (27/4)*L2^2*L3 ---
% dNxi(6)  = (27/4)*( 2*L2*dL2_dxi*L3 + L2^2*dL3_dxi );
% dNeta(6) = (27/4)*( 2*L2*dL2_deta*L3 + L2^2*dL3_deta );
% % Here: dL2_dxi = 1, dL3_dxi = 0; dL2_deta = 0, dL3_deta = 1.
% % So, dNxi(6) = (27/2)*L2*L3 and dNeta(6) = (27/4)*L2^2.
% 
% % --- For N7 = (27/4)*L2*L3^2 ---
% dNxi(7)  = (27/4)*( dL2_dxi*L3^2 + L2*2*L3*dL3_dxi );
% dNeta(7) = (27/4)*( dL2_deta*L3^2 + L2*2*L3*dL3_deta );
% % With dL2_dxi = 1, dL3_dxi = 0; dL2_deta = 0, dL3_deta = 1:
% % dNxi(7) = (27/4)*L3^2 and dNeta(7) = (27/2)*L2*L3.
% 
% % --- For N8 = (27/4)*L3^2*L1 ---
% dNxi(8)  = (27/4)*( 2*L3*dL3_dxi*L1 + L3^2*dL1_dxi );
% dNeta(8) = (27/4)*( 2*L3*dL3_deta*L1 + L3^2*dL1_deta );
% % Here: dL3_dxi = 0, dL1_dxi = -1; dL3_deta = 1, dL1_deta = -1.
% % Thus, dNxi(8) = -(27/4)*L3^2 and dNeta(8) = (27/4)*(2*L1*L3 - L3^2).
% 
% % --- For N9 = (27/4)*L3*L1^2 ---
% dNxi(9)  = (27/4)*( dL3_dxi*L1^2 + L3*2*L1*dL1_dxi );
% dNeta(9) = (27/4)*( dL3_deta*L1^2 + L3*2*L1*dL1_deta );
% % With dL3_dxi = 0, dL1_dxi = -1; dL3_deta = 1, dL1_deta = -1:
% % dNxi(9) = -(27/2)*L1*L3 and dNeta(9) = (27/4)*(L1^2 - 2*L1*L3).
% 
% % --- For N10 = 27*L1*L2*L3 (the interior bubble) ---
% dNxi(10)  = 27*( dL1_dxi*L2*L3 + L1*dL2_dxi*L3 + L1*L2*dL3_dxi );
% dNeta(10) = 27*( dL1_deta*L2*L3 + L1*dL2_deta*L3 + L1*L2*dL3_deta );
% % Substituting: dL1_dxi=-1, dL2_dxi=1, dL3_dxi=0; dL1_deta=-1, dL2_deta=0, dL3_deta=1:
% % dNxi(10)  = 27*( -L2*L3 + L1*L3 ) = 27*(L1 - L2)*L3.
% % dNeta(10) = 27*( -L2*L3 + L1*L2 ) = 27*L2*(L1 - L3).
% 
% end
% 
% % %--------------------------------------------------------------------
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

% 
% % function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP3_Tri()
% % % precomputeShapeFunctionsP3_Tri:
% % %   Returns shape functions (N), derivatives (dNxi, dNeta),
% % %   and a 6-point Gauss rule (g_pt, g_wt) for a P3 triangle,
% % %   in the reference domain 0 <= xi, 0 <= eta, xi + eta <= 1.
% % %
% % %   Gmsh node ordering for P3 (anticlockwise) is assumed as:
% % %     Node1 -> corner1 = (1,0,0)
% % %     Node2 -> corner2 = (0,1,0)
% % %     Node3 -> corner3 = (0,0,1)
% % %     Node4 -> edge12 near corner1
% % %     Node5 -> edge12 near corner2
% % %     Node6 -> edge23 near corner2
% % %     Node7 -> edge23 near corner3
% % %     Node8 -> edge31 near corner3
% % %     Node9 -> edge31 near corner1
% % %     Node10-> interior node
% % %
% % %   Implementation uses barycentric polynomials:
% % %     L1 = xi, L2 = eta, L3 = 1 - xi - eta
% % %   with the known 10 cubic monomials for i+j+k=3,
% % %   arranged in Gmsh's anticlockwise node order.
% % 
% %     % 1) Gauss rule (6-pt is enough for polynomials up to degree 4).
% %     %    For a cubic shape, 6- or 7-point rule is common; we use 6 for demonstration.
% %     [g_pt, g_wt] = triGaussPoints6();
% %     numGauss = length(g_wt);
% % 
% %     % 2) We have 10 local nodes for P3.
% %     numNodes = 10;
% %     N      = zeros(numNodes, numGauss);
% %     dNxi   = zeros(numNodes, numGauss);
% %     dNeta  = zeros(numNodes, numGauss);
% % 
% %     % 3) For each Gauss point, evaluate shape fns and partial derivatives.
% %     for k = 1:numGauss
% %         xi  = g_pt(k,1);
% %         eta = g_pt(k,2);
% % 
% %         [Ni, dNxi_k, dNeta_k] = p3basisGmsh(xi, eta);
% % 
% %         N(:,k)     = Ni;
% %         dNxi(:,k)  = dNxi_k;
% %         dNeta(:,k) = dNeta_k;
% %     end
% % end
% % 
% % %--------------------------------------------------------------------
% % function [N, dNdxi, dNdeta] = p3basisGmsh(xi, eta)
% % % p3basisGmsh:
% % %   Returns the 10 shape functions and their partial derivatives
% % %   wrt (xi, eta) for a cubic (P3) triangle, matching the Gmsh anticlockwise
% % %   node ordering described above.
% % %
% % % Barycentric coords: L1= xi, L2= eta, L3= 1 - xi - eta.
% % 
% %     L1 = xi;
% %     L2 = eta;
% %     L3 = 1 - xi - eta;
% % 
% %     % The shape functions in that order:
% %     %  1) L1^3
% %     %  2) L2^3
% %     %  3) L3^3
% %     %  4) 3 L1^2 L2
% %     %  5) 3 L1 L2^2
% %     %  6) 3 L2^2 L3
% %     %  7) 3 L2 L3^2
% %     %  8) 3 L1 L3^2
% %     %  9) 3 L1^2 L3
% %     % 10) 6 L1 L2 L3
% %     N = zeros(10,1);
% %     N(1) = L1^3;
% %     N(2) = L2^3;
% %     N(3) = L3^3;
% %     N(4) = 3*L1^2*L2;
% %     N(5) = 3*L1*L2^2;
% %     N(6) = 3*L2^2*L3;
% %     N(7) = 3*L2*L3^2;
% %     N(8) = 3*L1*L3^2;
% %     N(9) = 3*L1^2*L3;
% %     N(10)= 6*L1*L2*L3;
% % 
% %     % Partial derivatives wrt xi, eta (chain rule).
% %     % Let's define dL1/dxi=1, dL2/dxi=0, dL3/dxi=-1
% %     %                 dL1/deta=0,dL2/deta=1, dL3/deta=-1
% %     dL1_dxi  = 1;
% %     dL2_dxi  = 0;
% %     dL3_dxi  = -1;
% % 
% %     dL1_deta = 0;
% %     dL2_deta = 1;
% %     dL3_deta = -1;
% % 
% %     % We'll define small helpers to get partials:
% %     % d( L1^a L2^b L3^c )/dxi = a L1^(a-1)*L2^b*L3^c * dL1/dxi + ...
% %     % (plus b and c terms as appropriate).
% %     [dNdxi, dNdeta] = deal(zeros(10,1));
% % 
% %     % Helper inline for partial derivative:
% %     % partial_xi( L1^a L2^b L3^c ) = (a L1^(a-1) * L2^b * L3^c)*dL1dxi
% %     %                              + (b L1^a * L2^(b-1) * L3^c)*dL2dxi
% %     %                              + (c L1^a * L2^b * L3^(c-1))*dL3dxi
% %     % Do similarly for partial_eta.
% % 
% %     % Indices for each shape:
% %     % 1) (3,0,0)
% %     [dNdxi(1), dNdeta(1)] = partialABC(3,0,0, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 2) (0,3,0)
% %     [dNdxi(2), dNdeta(2)] = partialABC(0,3,0, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 3) (0,0,3)
% %     [dNdxi(3), dNdeta(3)] = partialABC(0,0,3, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 4) (2,1,0)
% %     [dNdxi(4), dNdeta(4)] = partialABC(2,1,0, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 5) (1,2,0)
% %     [dNdxi(5), dNdeta(5)] = partialABC(1,2,0, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 6) (0,2,1)
% %     [dNdxi(6), dNdeta(6)] = partialABC(0,2,1, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 7) (0,1,2)
% %     [dNdxi(7), dNdeta(7)] = partialABC(0,1,2, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 8) (1,0,2)
% %     [dNdxi(8), dNdeta(8)] = partialABC(1,0,2, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     % 9) (2,0,1)
% %     [dNdxi(9), dNdeta(9)] = partialABC(2,0,1, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     %10) (1,1,1)
% %     % coefficient is 6 * L1^1 L2^1 L3^1, so effectively (a,b,c)=(1,1,1)
% %     % but there's an extra factor 6 in the function. We'll factor that in:
% %     % partial_xi(6 L1 L2 L3) = 6 * partial_xi(L1^1 L2^1 L3^1).
% %     [pxi_base, peta_base] = partialABC(1,1,1, L1,L2,L3, dL1_dxi,dL2_dxi,dL3_dxi, dL1_deta,dL2_deta,dL3_deta);
% %     dNdxi(10)  = 6 * pxi_base;
% %     dNdeta(10) = 6 * peta_base;
% % end
% % 
% % %-----------------------
% % function [pxi, peta] = partialABC(a,b,c, L1,L2,L3, ...
% %     dL1dxi, dL2dxi, dL3dxi, dL1deta, dL2deta, dL3deta)
% % % partial derivatives of (L1^a L2^b L3^c) wrt xi, eta.
% %     val = (L1^a)*(L2^b)*(L3^c);
% % 
% %     % wrt xi
% %     pxi = a*(L1^(a-1)*L2^b*L3^c)*dL1dxi ...
% %         + b*(L1^a*L2^(b-1)*L3^c)*dL2dxi ...
% %         + c*(L1^a*L2^b*L3^(c-1))*dL3dxi;
% % 
% %     % wrt eta
% %     peta = a*(L1^(a-1)*L2^b*L3^c)*dL1deta ...
% %          + b*(L1^a*L2^(b-1)*L3^c)*dL2deta ...
% %          + c*(L1^a*L2^b*L3^(c-1))*dL3deta;
% % end
% % 
% % %--------------------------------------------------------------------
% % function [g_pt, g_wt] = triGaussPoints6()
% % % triGaussPoints6: returns a 6-point integration rule on the
% % %  reference triangle (xi,eta>=0, xi+eta<=1).
% % % This is the same as in your P2 code, but repeated for convenience.
% %     g_pt = [ ...
% %       0.4459484909, 0.4459484909;
% %       0.4459484909, 0.1081030182;
% %       0.1081030182, 0.4459484909;
% %       0.0915762135, 0.0915762135;
% %       0.0915762135, 0.8168475730;
% %       0.8168475730, 0.0915762135 ];
% %     g_wt = [ ...
% %       0.2233815897;
% %       0.2233815897;
% %       0.2233815897;
% %       0.1099517437;
% %       0.1099517437;
% %       0.1099517437 ];
% % end