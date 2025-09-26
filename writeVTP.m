function writeVTP(x, y, scalars, filename)
    % writeVTP writes an XML VTK PolyData (VTP) file.
    %   x, y      : vectors of point coordinates (z is set to 0)
    %   scalars   : vector of scalar data for each point (e.g. velocity magnitude)
    %   filename  : output file name (should end with .vtp)
    
    num_points = length(x);
    
    % Create connectivity and offsets for vertices (one vertex per point)
    connectivity = 0:(num_points-1);
    offsets = 1:num_points;
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file %s for writing.', filename);
    end
    
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '  <PolyData>\n');
    fprintf(fid, '    <Piece NumberOfPoints="%d" NumberOfVerts="%d">\n', num_points, num_points);
    
    %% Write Points
    fprintf(fid, '      <Points>\n');
    fprintf(fid, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
    for i = 1:num_points
        fprintf(fid, '          %f %f %f\n', x(i), y(i), 0);
    end
    fprintf(fid, '        </DataArray>\n');
    fprintf(fid, '      </Points>\n');
    
    %% Write Vertices
    fprintf(fid, '      <Verts>\n');
    % Connectivity
    fprintf(fid, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fid, '          ');
    fprintf(fid, '%d ', connectivity);
    fprintf(fid, '\n');
    fprintf(fid, '        </DataArray>\n');
    % Offsets
    fprintf(fid, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid, '          ');
    fprintf(fid, '%d ', offsets);
    fprintf(fid, '\n');
    fprintf(fid, '        </DataArray>\n');
    fprintf(fid, '      </Verts>\n');
    
    %% Write PointData (scalar field)
    fprintf(fid, '      <PointData Scalars="u_value">\n');
    fprintf(fid, '        <DataArray type="Float32" Name="u_value" format="ascii">\n');
    for i = 1:num_points
        fprintf(fid, '          %f\n', scalars(i));
    end
    fprintf(fid, '        </DataArray>\n');
    fprintf(fid, '      </PointData>\n');
    
    fprintf(fid, '    </Piece>\n');
    fprintf(fid, '  </PolyData>\n');
    fprintf(fid, '</VTKFile>\n');
    
    fclose(fid);
end