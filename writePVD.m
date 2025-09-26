function writePVD(timeSteps, baseFilename, pvdFilename)
% writePVD creates a ParaView Data (PVD) file that references a collection
% of VTK (legacy) files, each associated with a given time step.
%
%   timeSteps   : Vector of time values for each snapshot
%   baseFilename: Filename pattern for the .vtk files, e.g., 'solution_%04d.vtk'
%   pvdFilename : Output PVD filename (e.g., 'solutions.pvd')

    fid = fopen(pvdFilename, 'w');
    if fid == -1
        error('Cannot open file %s for writing.', pvdFilename);
    end

    % XML header
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n');
    fprintf(fid, '  <Collection>\n');

    % Each snapshot as a DataSet
    for i = 1:length(timeSteps)
        thisFile = sprintf(baseFilename, i);
        fprintf(fid, ...
            '    <DataSet timestep="%.6g" part="0" file="%s"/>\n', ...
            timeSteps(i), thisFile);
    end

    fprintf(fid, '  </Collection>\n');
    fprintf(fid, '</VTKFile>\n');
    fclose(fid);

    fprintf('Wrote PVD file: %s\n', pvdFilename);
end
