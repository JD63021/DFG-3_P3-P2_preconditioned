function uv_new = update_bc(uv, boundaryInfo, nodeInfo, Nxy, t, corner, bcFlags,inletProfile)
% update_bc applies Dirichlet boundary conditions for the driven cavity.
%
% Inputs:
%   uv           : solution vector [2*Nxy + Npr]
%   boundaryInfo : structure with boundary info (fields like 'flag_2', etc.)
%   nodeInfo     : structure with node coordinates.
%   Nxy          : number of velocity DOFs per component.
%   t            : current time.
%   corner       : index for pressure correction.
%   bcFlags      : structure with fields:
%                    .top   - flag number for top boundary.
%                    .sides - array of flag numbers for side/bottom boundaries.
%
% Output:
%   uv_new       : solution vector with boundary conditions applied.

uv_new = uv;
H = max(nodeInfo.velocity.y) - min(nodeInfo.velocity.y);

% --- Apply top boundary condition (Dirichlet inlet)
topFlag = bcFlags.inlet;
if isfield(boundaryInfo, ['flag_' num2str(topFlag)])
    inletNodes = boundaryInfo.(['flag_' num2str(topFlag)])(:).';
else
    inletNodes = [];
end
for nodeID = inletNodes
    rowX = globalRow(nodeID, 'x', Nxy);
    rowY = globalRow(nodeID, 'y', Nxy);
    Uin = inletProfile(t, nodeInfo.velocity.y(nodeID), H);
    uv_new(rowX) = Uin;
    uv_new(rowY) = 0;
end

% --- Apply side and bottom boundary conditions (zero velocity)
sides = bcFlags.wall;
sideNodes = [];
for i = 1:length(sides)
    flagStr = ['flag_' num2str(sides(i))];
    if isfield(boundaryInfo, flagStr)
        sideNodes = [sideNodes; boundaryInfo.(flagStr)(:)];
    end
end
sideNodes = unique(sideNodes);
for nodeID = sideNodes'
    rowX = globalRow(nodeID, 'x', Nxy);
    rowY = globalRow(nodeID, 'y', Nxy);
    uv_new(rowX) = 0;
    uv_new(rowY) = 0;
end

% --- Fix corner pressure DOF
cornerGlobal = 2*Nxy + corner;
uv_new(cornerGlobal) = 0;
end

function row = globalRow(vNode, comp, Nxy)
% globalRow computes the global index for a given velocity node and component.
switch comp
    case 'x'
        row = vNode;
    case 'y'
        row = vNode + Nxy;
    otherwise
        error('Invalid component. Use "x" or "y".');
end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % A. UPDATED update_bc.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function uv_new = update_bc(uv, boundaryInfo,nodeInfo, Nxy,t,corner)
% %UPDATE_BC  Applies Dirichlet BC for a lid-driven cavity:
% %   - top boundary => x-velocity=1, y-velocity=0
% %   - left, right, bottom => x=0, y=0
% %
% % Inputs:
% %   uv           : vector [2*Nxy + Npr], or at least [2*Nxy]
% %   boundaryInfo : from mesh2_gmsh, with fields like flag_5, flag_6, etc.
% %   Nxy          : number of Q2 velocity nodes
% %
% % Output:
% %   uv_new       : same vector but with velocity BC enforced
% %
% 
% uv_new = uv;
% inletNodes = boundaryInfo.flag_2(:).';
% % Umax   = 1.5;
% % H= 0.41;
% H = max(nodeInfo.velocity.y) - min(nodeInfo.velocity.y);
% %% (A) Top boundary => x=1, y=0
% % topNodes = boundaryInfo.flag_2(:)';  % 8 => top
% % 4*1.5*sin(pi*t/8)
% for nodeID = inletNodes
%     rowX = globalRow(nodeID, 'x', Nxy);  % x-velocity DOF
%     rowY = globalRow(nodeID, 'y', Nxy);  % y-velocity DOF
% 
%     yy   = nodeInfo.velocity.y(nodeID);  % the physical y-coordinate
%     % Parabolic formula (adjust as you wish):
%     % Uin = 1;
%     Uin  = 4*1.5*sin(pi*t/8)*( yy*(H - yy) ) / (H^2);
%     % uv_new(rowX) = 1.5*sin(pi*t/2);
%     uv_new(rowX) = Uin;% x=1
%     uv_new(rowY) = 0;
%     % uv_new(rowY) = 0.3*sin(pi*t/0.5);   % y=0
% end
% 
% %% (B) Left/Right/Bottom => x=0, y=0
% LRB = unique([boundaryInfo.flag_1(:); boundaryInfo.flag_5(:); boundaryInfo.flag_4(:);boundaryInfo.flag_6(:); boundaryInfo.flag_7(:);boundaryInfo.flag_8(:);]);
% for nodeID = LRB(:)'
%     rowX = globalRow(nodeID, 'x', Nxy);
%     rowY = globalRow(nodeID, 'y', Nxy);
%     uv_new(rowX) = 0;
%     uv_new(rowY) = 0;
% end
% 
% % R = unique([boundaryInfo.flag_5(:)]);
% % for nodeID = R(:)'
% %     rowX = globalRow(nodeID, 'x', Nxy);
% %     rowY = globalRow(nodeID, 'y', Nxy);
% %     uv_new(rowX) = 0;
% %     uv_new(rowY) = 1;
% % end
% 
% 
% 
% %% (C) Fix corner pressure if you like:
% cornerIndex =  corner; %or find actual "corner pressure node ID"
% cornerGlobal = 2*Nxy + cornerIndex;
% uv_new(cornerGlobal) = 1;
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % A small helper to convert (velocity nodeID, 'x'/'y', Nxy) => global row index
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function row = globalRow(vNode, comp, Nxy)
% switch comp
%     case 'x'
%         row = vNode;            % x-velocity DOFs = [1..Nxy]
%     case 'y'
%         row = vNode + Nxy;      % y-velocity DOFs = [Nxy+1..2*Nxy]
%     otherwise
%         error('Invalid velocity component (use ''x'' or ''y'')');
% end
% end
% 
