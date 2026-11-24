function ex_writecsv(file, header, M)
%EX_WRITECSV  Write a numeric matrix with a header row (MATLAB and Octave).
fid = fopen(file, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fmt = [repmat('%.6g,', 1, size(M, 2) - 1), '%.6g\n'];
fprintf(fid, fmt, M');
fclose(fid);
end
