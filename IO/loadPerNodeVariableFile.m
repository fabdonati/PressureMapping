function filecontents = loadPerNodeVariableFile(pathtofile,nr_of_points,dimensionOfValues)
% http://vis.lbl.gov/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf        

% Default options:
    if nargin<3
        dimensionOfValues=1;
    end
    
    fid = fopen(pathtofile,'r');
    [~,~,mf] = fopen(fid);
    filecontents.description = fread(fid,[80 1],'*char')'; % description lines                           2*80 chars

    while true && ~feof(fid)
        part   = fread(fid,[1 80],'*char');              % [part                                               80 chars
        if ~isempty(part) && strcmp(part(1:4),'part')  % loop over all parts of the geometry; surface mesh and elements
            % STEP 1: LOAD ALL DATA FOR THE CURRENT ELEMENT
            castingtype = '*int32';
            castingendi = 'b';
            part_number    = fread(fid,                 1, castingtype,castingendi); % #                                                    1 int
            part_descr     = fread(fid,            [1 80], '*char'      ); % description line                                    80 chars
            if strcmp(part_descr(1:5),'block')
                filecontents.data  = fread(fid,[nr_of_points dimensionOfValues], '*float' ,castingendi);
            end
        end
    end

