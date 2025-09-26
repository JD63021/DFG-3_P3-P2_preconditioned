function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh( ...
    gmshVelFile, gmshPresFile, ...
    excludedVelFlags, pBoundFlags, pBoundVals )
% mesh_p3p2_gmsh:
%   Reads two Gmsh files for a P3–P2 triangular approach:
%     (1) Velocity mesh => P3 => TRIANGLES10, boundary => LINES4
%     (2) Pressure mesh => P2 => TRIANGLES6, boundary => LINES3
%
%   Builds:
%     nodeInfo.velocity (P3) + nodeInfo.pressure (P2),
%     elemInfo.velElements (10 nodes) + elemInfo.presElements (6 nodes),
%     boundaryInfo.velLine4Elements.flag_<X> => Nx4 lines from velocity mesh,
%     boundaryInfo.presLine3Elements.flag_<X> => Nx3 lines from pressure mesh.
%
%   excludedVelFlags => optional set of boundary flags to ignore
%   pBoundFlags, pBoundVals => optional array of flags and boundary values for pressure
%
%   Returns nodeInfo, elemInfo, boundaryInfo structures, analogous to your older code.

if nargin < 3, excludedVelFlags = []; end
if nargin < 4, pBoundFlags = []; end
if nargin < 5, pBoundVals = []; end

%% (A) Read Velocity Mesh => P3 => TRIANGLES10
if ischar(gmshVelFile)
    run(gmshVelFile);  % populates variable 'msh'
    msh_vel = msh;
else
    error('Velocity mesh file must be a .m string from Gmsh.');
end

if ~isfield(msh_vel, 'TRIANGLES10')
    error('Expected TRIANGLES10 in the velocity mesh (P3).');
end
triElementsVel = msh_vel.TRIANGLES10(:,1:10);
numElemsVel    = size(triElementsVel,1);

allVelCoords = msh_vel.POS;  % Nx2 or Nx3
numVelNodes  = size(allVelCoords,1);
vel_nodes.id = (1:numVelNodes).';
vel_nodes.x  = allVelCoords(:,1);
vel_nodes.y  = allVelCoords(:,2);

%% (B) Read Pressure Mesh => P2 => TRIANGLES6
if ischar(gmshPresFile)
    run(gmshPresFile);
    msh_pres = msh;
else
    error('Pressure mesh file must be a .m string from Gmsh.');
end

if ~isfield(msh_pres, 'TRIANGLES6')
    error('Expected TRIANGLES6 in the pressure mesh (P2).');
end
triElementsPres = msh_pres.TRIANGLES6(:,1:6);
numElemsPres    = size(triElementsPres,1);

allPresCoords = msh_pres.POS;
numPresNodes  = size(allPresCoords,1);
pres_nodes.id = (1:numPresNodes).';
pres_nodes.x  = allPresCoords(:,1);
pres_nodes.y  = allPresCoords(:,2);

%% (C) Store element connectivity
elemInfo.velElements  = triElementsVel;  % 10 nodes each
elemInfo.presElements = triElementsPres; % 6 nodes each

%% (D) Fill nodeInfo
nodeInfo.velocity = vel_nodes;
nodeInfo.pressure = pres_nodes;

%% (E) Build boundaryInfo

boundaryInfo = struct();

% 1) Velocity boundaries => LINES4 => Nx5 [nodeA, nodeB, mid1, mid2, flag]
%    or sometimes the order is [nA, nM1, nM2, nB, flag].
%    Check exactly how Gmsh emits them for P3 lines. We'll assume the last column is the flag.
if isfield(msh_vel,'LINES4') && ~isempty(msh_vel.LINES4)
    lines4_vel = msh_vel.LINES4; % Nx5 => [nodeA nodeM1 nodeM2 nodeB flag]
    uniqueVFlags = unique(lines4_vel(:,5));
    velLine4Struct = struct();
    for iF = 1:numel(uniqueVFlags)
        gFlag = uniqueVFlags(iF);
        if ismember(gFlag, excludedVelFlags), continue; end
        mask = (lines4_vel(:,5) == gFlag);
        theseLines = lines4_vel(mask, 1:4);  % columns 1..4 => node IDs
        fName = sprintf('flag_%d', gFlag);
        velLine4Struct.(fName) = theseLines;
    end
    boundaryInfo.velLine4Elements = velLine4Struct;
else
    warning('No LINES4 found in the velocity mesh => boundaryInfo.velLine4Elements = empty.');
    boundaryInfo.velLine4Elements = struct();
end

% 2) Pressure boundaries => LINES3 => Nx4 [nodeA nodeMid nodeB flag]
%    or if the pressure boundary is truly P2, Gmsh might do that. 
%    If your Gmsh uses LINES for linear edges, adapt accordingly.
if isfield(msh_pres,'LINES3') && ~isempty(msh_pres.LINES3)
    lines3_pres = msh_pres.LINES3; 
    uniquePFlags = unique(lines3_pres(:,4));
    presLine3Struct = struct();
    for iF = 1:numel(uniquePFlags)
        pFlag = uniquePFlags(iF);
        maskP = (lines3_pres(:,4) == pFlag);
        theseLinesP = lines3_pres(maskP, 1:3);  % nodeA, nodeMid, nodeB
        fNameP = sprintf('flag_%d', pFlag);
        presLine3Struct.(fNameP) = theseLinesP;
    end
    boundaryInfo.presLine3Elements = presLine3Struct;
else
    warning('No LINES3 found in pressure mesh => boundaryInfo.presLine3Elements = empty.');
    boundaryInfo.presLine3Elements = struct();
end

%% (F) Possibly store velocity boundary node sets
boundaryInfo.allVelNodes = [];
if isfield(boundaryInfo,'velLine4Elements')
    fnames = fieldnames(boundaryInfo.velLine4Elements);
    for iF = 1:numel(fnames)
        thisField = fnames{iF};
        lineEls   = boundaryInfo.velLine4Elements.(thisField);
        if isempty(lineEls), continue; end
        boundaryNodes = unique(lineEls(:));
        boundaryInfo.(thisField) = boundaryNodes;  % store them as well
        boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
    end
    boundaryInfo.allVelNodes = unique(boundaryInfo.allVelNodes);
end

%% (G) Pressure BC if needed
boundaryInfo.pressureConditions = struct('flag',{},'nodes',{},'value',{});
if ~isempty(pBoundFlags)
    for iPf = 1:numel(pBoundFlags)
        thisFlag  = pBoundFlags(iPf);
        thisValue = pBoundVals(iPf);
        boundaryInfo.pressureConditions(iPf).flag  = thisFlag;
        boundaryInfo.pressureConditions(iPf).nodes = [];
        boundaryInfo.pressureConditions(iPf).value = thisValue;
    end
end

end


% function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh( ...
%     gmshVelFile, gmshPresFile, ...
%     excludedVelFlags, pBoundFlags, pBoundVals)
% % mesh6_gmsh reads two Gmsh files for a P2–P1 triangular approach:
% %   (1) Velocity mesh => P2 => TRIANGLES6, boundary => LINES3
% %   (2) Pressure mesh => P1 => TRIANGLES,  boundary => LINES
% %
% % And it builds:
% %   nodeInfo.velocity (P2) + nodeInfo.pressure (P1),
% %   elemInfo.velElements (6 nodes) + elemInfo.presElements (3 nodes),
% %   boundaryInfo.velLine3Elements.flag_<X> => Nx3 lines from velocity mesh,
% %   boundaryInfo.presLine2Elements.flag_<X> => Nx2 lines from pressure mesh.
% %
% % If 'excludedVelFlags' is specified, those velocity boundary segments are ignored.
% % If pBoundFlags/vals are specified, the code sets pressure boundary info as well.
% 
% if ~exist('excludedVelFlags','var') || isempty(excludedVelFlags)
%     excludedVelFlags = [];
% end
% if ~exist('pBoundFlags','var') || isempty(pBoundFlags)
%     pBoundFlags = [];
% end
% if ~exist('pBoundVals','var') || isempty(pBoundVals)
%     pBoundVals = [];
% end
% 
% %% (A) Read Velocity Mesh => P2 => TRIANGLES6
% if ischar(gmshVelFile)
%     run(gmshVelFile);  % loads struct 'msh' into workspace
%     msh_vel = msh;
% else
%     error('Velocity mesh file not provided as a string.');
% end
% 
% if ~isfield(msh_vel,'TRIANGLES6')
%     error('Expected TRIANGLES6 in the velocity mesh (P2).');
% end
% 
% triElementsVel = msh_vel.TRIANGLES6(:,1:6);
% numElemsVel    = size(triElementsVel,1);
% 
% allVelCoords = msh_vel.POS;    % Nx3 or Nx2
% numVelNodes  = size(allVelCoords,1);
% 
% vel_nodes.id = (1:numVelNodes).';
% vel_nodes.x  = allVelCoords(:,1);
% vel_nodes.y  = allVelCoords(:,2);
% 
% %% (B) Read Pressure Mesh => P1 => TRIANGLES
% if ischar(gmshPresFile)
%     run(gmshPresFile);
%     msh_pres = msh;
% else
%     error('Pressure mesh file not provided as a string.');
% end
% 
% if ~isfield(msh_pres,'TRIANGLES')
%     error('Expected TRIANGLES in the pressure mesh (P1).');
% end
% 
% triElementsPres = msh_pres.TRIANGLES(:,1:3);
% numElemsPres    = size(triElementsPres,1);
% 
% allPresCoords = msh_pres.POS;
% numPresNodes  = size(allPresCoords,1);
% 
% pres_nodes.id = (1:numPresNodes).';
% pres_nodes.x  = allPresCoords(:,1);
% pres_nodes.y  = allPresCoords(:,2);
% 
% %% (C) Store element connectivity
% elemInfo.velElements  = triElementsVel;  % 6 nodes each
% elemInfo.presElements = triElementsPres; % 3 nodes each
% 
% %% (D) Fill nodeInfo
% nodeInfo.velocity = vel_nodes;
% nodeInfo.pressure = pres_nodes;
% 
% %% (E) Build boundaryInfo
% boundaryInfo = struct();
% 
% % 1) Velocity boundaries => LINES3 => Nx4 [nodeA nodeMid nodeB flag]
% if ~isfield(msh_vel,'LINES3') || isempty(msh_vel.LINES3)
%     warning('No LINES3 found in velocity mesh => boundaryInfo.velLine3Elements = empty.');
%     boundaryInfo.velLine3Elements = struct();
% else
%     lines3_vel = msh_vel.LINES3; % Nx4 => [nA nM nB flag]
%     uniqueVFlags = unique(lines3_vel(:,4));
%     velLine3Struct = struct();
%     for iF = 1:numel(uniqueVFlags)
%         gFlag = uniqueVFlags(iF);
%         if ismember(gFlag, excludedVelFlags)
%             continue;
%         end
%         mask = (lines3_vel(:,4) == gFlag);
%         theseLines = lines3_vel(mask, 1:3);  % columns 1..3 => node IDs
%         fieldName  = sprintf('flag_%d', gFlag);
%         velLine3Struct.(fieldName) = theseLines;
%     end
%     boundaryInfo.velLine3Elements = velLine3Struct;
% end
% 
% % 2) Pressure boundaries => LINES => Nx3 [nodeA nodeB flag]
% if ~isfield(msh_pres,'LINES') || isempty(msh_pres.LINES)
%     warning('No LINES found in pressure mesh => boundaryInfo.presLine2Elements = empty.');
%     boundaryInfo.presLine2Elements = struct();
% else
%     lines2_pres = msh_pres.LINES; % Nx3 => [nA nB flag]
%     uniquePFlags = unique(lines2_pres(:,3));
%     presLine2Struct = struct();
%     for iF = 1:numel(uniquePFlags)
%         pFlag = uniquePFlags(iF);
%         maskP = (lines2_pres(:,3) == pFlag);
%         theseLinesP = lines2_pres(maskP,1:2);
%         fNameP = sprintf('flag_%d', pFlag);
%         presLine2Struct.(fNameP) = theseLinesP;
%     end
%     boundaryInfo.presLine2Elements = presLine2Struct;
% end
% 
% %% (F) Optionally store velocity boundary node sets
% boundaryInfo.allVelNodes = [];
% if isfield(boundaryInfo,'velLine3Elements')
%     fnames = fieldnames(boundaryInfo.velLine3Elements);
%     for iF = 1:numel(fnames)
%         thisField = fnames{iF}; % e.g. 'flag_5'
%         lineEls   = boundaryInfo.velLine3Elements.(thisField);
%         if isempty(lineEls), continue; end
%         boundaryNodes = unique(lineEls(:));
%         boundaryInfo.(thisField) = boundaryNodes;  % store them
%         boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
%     end
%     boundaryInfo.allVelNodes = unique(boundaryInfo.allVelNodes);
% end
% 
% %% (G) Pressure boundary conditions (if requested)
% boundaryInfo.pressureConditions = struct('flag',{},'nodes',{},'value',{});
% if ~isempty(pBoundFlags)
%     for iPf = 1:numel(pBoundFlags)
%         thisFlag  = pBoundFlags(iPf);
%         thisValue = pBoundVals(iPf);
%         boundaryInfo.pressureConditions(iPf).flag  = thisFlag;
%         boundaryInfo.pressureConditions(iPf).nodes = [];  % or find from LINES
%         boundaryInfo.pressureConditions(iPf).value = thisValue;
%     end
% end
% 
% end
% 
