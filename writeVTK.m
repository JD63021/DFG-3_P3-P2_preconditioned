function writeVTK(x, y, u, filename)
% writeVTK writes out a legacy VTK PolyData file containing:
%   - A set of points at (x(i), y(i), 0.0)
%   - Each point is stored as a vertex cell
%   - A scalar field "u_value"
%
% x, y      : Arrays of coordinates for each point (1D arrays of same length)
% u         : Scalar field (1D array, same length as x,y)
% filename  : Output filename (e.g. 'solution_0000.vtk')

    fid = fopen(filename, 'w');
    if fid < 0
        error('Cannot open file %s for writing.', filename);
    end

    % --- VTK File Header (PolyData) ---
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK output from MATLAB\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n\n');

    numPoints = length(x);

    % --- Write the Point Coordinates ---
    fprintf(fid, 'POINTS %d float\n', numPoints);
    for i = 1:numPoints
        % z = 0.0 to stay in XY plane
        fprintf(fid, '%f %f 0.0\n', x(i), y(i));
    end
    fprintf(fid, '\n');

    % --- Define Each Point as a Separate VERTEX ---
    % VERTICES <number_of_vertices> <total_number_of_indices>
    % Each line: [1  point_index]
    fprintf(fid, 'VERTICES %d %d\n', numPoints, 2*numPoints);
    for i = 0:numPoints-1
        fprintf(fid, '1 %d\n', i);
    end
    fprintf(fid, '\n');

    % --- Write the Scalar Data ---
    fprintf(fid, 'POINT_DATA %d\n', numPoints);
    fprintf(fid, 'SCALARS u_value float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i = 1:numPoints
        fprintf(fid, '%f\n', u(i));
    end

    fclose(fid);
    fprintf('VTK file "%s" successfully written (PolyData format).\n', filename);
end
