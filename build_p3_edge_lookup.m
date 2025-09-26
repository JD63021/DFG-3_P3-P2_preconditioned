function boundaryInfo = build_p3_edge_lookup(elemInfo, boundaryInfo)
% Build O(1) lookup from a boundary LINES4 segment (its two corner nodes)
% to (element id, local side id) for P3 triangles (10 nodes).
%
% Output:
%   boundaryInfo.edgeMap.flag_<F>  =  [elemIdx, sideIdx] for each LINES4 row
%
% Assumptions:
%   - elemInfo.velElements is (nEl x 10), with the 3 triangle vertices in cols 1..3
%   - boundaryInfo.velLine4Elements.flag_<F> is (nSeg x 4): [A M1 M2 B]
%
% This eliminates repeated ismember-based searches in force_int.

T = elemInfo.velElements;         % nEl x 10 (P3)
nEl = size(T,1);
corner = T(:,1:3);                % vertices only
pairs  = [1 2; 2 3; 3 1];         % local sides

% Keying: pack a sorted (i,j) node pair into a uint64 key = i*BASE + j
maxNode = max(T(:));
BASE = uint64(maxNode) + 1;

% Build map (i,j) --> (elem, side) packed as uint32: val = (elem-1)*3 + side
K = containers.Map('KeyType','uint64', 'ValueType','uint32');
for e = 1:nEl
    v = corner(e,:);
    for s = 1:3
        i = v(pairs(s,1));  j = v(pairs(s,2));
        if i > j, tmp=i; i=j; j=tmp; end
        key = uint64(i)*BASE + uint64(j);
        K(key) = uint32((e-1)*3 + s);
    end
end

% For every boundary flag, translate its LINES4 into [elem,side] rows.
edgeMap = struct();
if isfield(boundaryInfo,'velLine4Elements')
    flds = fieldnames(boundaryInfo.velLine4Elements);
    for k = 1:numel(flds)
        fn = flds{k};                      % e.g., 'flag_5'
        L  = boundaryInfo.velLine4Elements.(fn);   % nSeg x 4 [A M1 M2 B]
        if isempty(L), continue; end
        A = L(:,1); B = L(:,4);
        i = min(A,B);  j = max(A,B);
        keys = uint64(i)*BASE + uint64(j);

        val = zeros(numel(keys),1,'uint32');
        for r = 1:numel(keys)
            val(r) = K(keys(r));          % must exist if mesh is consistent
        end
        elem  = floor(double(val-1)/3) + 1;
        side  = mod(double(val-1),3) + 1;
        edgeMap.(fn) = [elem(:), side(:)];
    end
end

boundaryInfo.edgeMap = edgeMap;
boundaryInfo._edgeMapBASE = BASE;   % stash for possible checks
end
