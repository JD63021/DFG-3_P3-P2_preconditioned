function testCase = fpc()
% fpc - Test case configuration for the flow past cylinder problem.
%
% This file specifies:
%   (i) the gmsh files to run (now two files: one for velocity and one for pressure),
%   (ii) any velocity boundary flags to exclude (if desired),
%   (iii) pressure boundary flags and values (if desired), and 
%   (iv) the boundary flags to use for enforcing Dirichlet conditions.
%
% For fpc, we choose:
%   - Top boundary flag = 2.
%   - Typical side and bottom boundary flags = [1, 4, 5, 6, 7, 8].
%
% Note: fpc2q (velocity grid) and fpc2l (pressure grid) replace the older fpc2r.

    testCase.gmshVelFile  = 'fpct2p3.m';  % gmsh file for velocity grid (Q2)
    testCase.gmshPresFile = 'fpct2p2.m';  % gmsh file for pressure grid (Q1)
    
    testCase.excludedVelFlags = [2 9 10 11 12];             % velocity boundary flags to exclude
    testCase.pBoundFlags      = [];                % (no pressure BCs in this case)
    testCase.pBoundVals       = [];
    testCase.boundaryFlags.inlet = 4;              % top boundary flag (for inlet)
    testCase.boundaryFlags.wall  = [1 2 3 5 6 7 8 ]; % side and bottom boundaries

    % Specify the inlet velocity profile (function handle).
    % This function should accept: (t, y, H) and return the x-velocity.
    testCase.inletProfile   = @(t, y, H) 4*1.5*sin(pi*t/8)*(y*(H-y))/(H^2);
    testCase.inletProfileSS = @(t, y, H) 4*0.3*(y*(H-y))/(H^2);
    testCase.corner = 2;%113;
end


% function testCase = fpc()
% % fpc - Test case configuration for the flwo past cylinder problem.
% %
% % This file specifies:
% %   (i) the gmsh file to run (here, 'fpc.m'),
% %   (ii) any velocity boundary flags to exclude (if desired),
% %   (iii) pressure boundary flags and values (if desired), and 
% %   (iv) the boundary flags to use for enforcing Dirichlet conditions.
% %
% % For fpc, we choose:
% %   - Top boundary flag = 2.
% %   - Typical Side and bottom boundary flags = [1, 4, 5, 6, 7, 8].
% % fpc2r 2- inlet  3- outlet corner-2 
% % fpc11 4- inlet , 2 -outlet corner -2
% % fpc12 4- inlet , 2 -outlet corner -2
% % fpc16 4- inlet , 2 -outlet corner -205
% % cyl2  19-inlet, 28-33 do nothing (33 outlet), 20-27 walls
% % fpc 18: 13- inlet, 24 - outlet, 15, 21- do nothing,  cylinder walls 17,18, other walls- 14, 16 ,19 ,20 ,22, 23, corner 4
% % fpc_r: 7 inlet, 9 outlet, 8, 10 walls, 11 12 cylinder walls 4 corner
% 
%     testCase.gmshFile = 'fpc2r.m';           % Use the fpc gmsh file
%     testCase.excludedVelFlags = [3];             % (none to exclude)
%     testCase.pBoundFlags = [];                  % (no pressure BCs in this case)
%     testCase.pBoundVals  = [];
%     testCase.boundaryFlags.inlet = 2;             % top boundary flag (for inlet)
%     testCase.boundaryFlags.wall = [1, 4, 5, 6, 7, 8]; % side and bottom boundaries
% 
%     % Specify the inlet velocity profile (function handle).
%     % This function should accept: (t, y, H) and return the x-velocity.
%     testCase.inletProfile = @(t, y, H)4*0.3*(y*(H-y))/(H^2);%  4*0.3*((y+0.2)*(0.21-y))/(0.41^2);
%     testCase.inletProfileSS = @(t, y, H) 4*0.3*(y*(H-y))/(H^2);% 4*0.3*((y+0.2)*(0.21-y))/(0.41^2);
%     testCase.corner = 2;
% end
