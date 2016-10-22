function geofilecontents = loadGeoFile(pathtofile,istimeresolved_simulation)
% http://vis.lbl.gov/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf        
% http://vis.lbl.gov/NERSC/Software/ensight/docs8/UserManual.pdf
        % load geometry file; see below for details on the file format
        fid = fopen(pathtofile,'r');
        [~,~,mf] = fopen(fid);
        geofilecontents.type   = fread(fid,[1 80],'*char'); % C Binary                                            80 chars
        if istimeresolved_simulation == 1
            [~]                         = fread(fid,[1 80],'*char') ; % BEGIN TIME STEP                               80 chars % IGNORE
        end
        geofilecontents.description = fread(fid,[80 2],'*char')'; % description lines                           2*80 chars
        
        [~]                         = fread(fid,[1 80],'*char'); % node id <off/given/assign/ignore>                   80 chars % IGNORE
        [~]                         = fread(fid,[1 80],'*char'); % element id <off/given/assign/ignore>                80 chars % IGNORE
        
        if ~strcmp(geofilecontents.type(1:8),'C Binary'), error('This code only works for C binary data'); end
        
        index = 1;
        while true && ~feof(fid)
            part   = fread(fid,[1 80],'*char');              % [part                                               80 chars
            if ~isempty(part) && strcmp(part(1:4),'part')  % loop over all parts of the geometry; surface mesh and elements
                
                % STEP 1: LOAD ALL DATA FOR THE CURRENT ELEMENT
                castingtype = '*int32';
                castingendi = 'b';
                part_number    = fread(fid,                 1, castingtype,castingendi); % #                                                    1 int
                part_descr     = fread(fid,            [1 80], '*char'      ); % description line                                    80 chars
                coords            = fread(fid,            [1 80], '*char'      ); % coordinates                                         80 chars % IGNORE
                if strcmp(coords(1:5),'block')
                    % this is now followed by:
                    % i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1) 3 ints
                    % [imin imax jmin jmax kmin kmax] # if range used: 6 ints
                    geofilecontents.dimensions = fread(fid,                 3, castingtype,castingendi);
                    geofilecontents.extents = fread(fid,                 6, 'double',castingendi);
                    geofilecontents.nr_of_points = prod(geofilecontents.dimensions);
                    geofilecontents.points  = fread(fid,[geofilecontents.nr_of_points 3], '*float' ,castingendi);
                else
                    nr_of_vertices = fread(fid,                 1, castingtype,castingendi); % nn                                                   1 int
                    vertices_xyz   = fread(fid,[nr_of_vertices 3], '*float' ,castingendi); % x y z coordinates // index is skipped                
                    faces_type     = strcat(fread(fid,            [1 80], '*char'      )); % get the type of element
                    nr_of_faces    = fread(fid,                 1, castingtype,castingendi); % ne                                                   1 int
                    switch faces_type(1:5)  % compare first 5 characters in faces_type
                        case {'tetra'}, cornerpts = 4;    % tetraeder; 4 cornerpoints
                        case 'tria3', cornerpts = 3;    % triangle;  3 cornerpoints
                        case 'penta', cornerpts = 6;
                        case 'quad4', cornerpts = 4;
                        otherwise, error('Unknown elementtype, only tetra4 and tria3 are implemented.');
                    end
                    element_faces = fread(fid,[cornerpts nr_of_faces],castingtype,castingendi)';

                    % STEP 2: SAVE DATA FOR CURRENT ELEMENT IN STRUCT
                    geofilecontents.(part_descr).vertices     =  vertices_xyz;
                    geofilecontents.(part_descr).(['faces_' faces_type]) = element_faces;
                    geofilecontents.(part_descr).part_nr      =   part_number;

                    index = index + 1; % next element
                end
            elseif ~isempty(part) && strcmp(part(1:8),'END TIME') % disp('END TIME STEP reached; check if this the end of the file?')
                curpos  = ftell(fid); %save current position
                fread(fid,1); %try reading 1 more line
                if feof(fid) % check if eof is reached
                    geofilecontents.nr_parts = index; % save nr of parts
                    break; %eof reached; break while loop;
                else % No; return to previous curpos
                    if (fseek(fid,curpos,'bof') ~= 0)
                        error('could not go to defined position...')
                    end
                end
            elseif length(part) >= 6 && strcmp(part(1:6),'penta6') % sometimes multiple element face_types are used; load PENTA6 elements too!
                faces_type = strcat(part);
                cornerpts = 6;
                nr_of_pentas    = fread(fid,                 1, castingtype,castingendi);
                geofilecontents.(part_descr).(['faces_' faces_type]) = fread(fid,[cornerpts nr_of_pentas],castingtype,castingendi)';

            elseif length(part) >= 5 && strcmp(part(1:5),'quad4') % sometimes multiple element face_types are used; load QUAD4 elements too!
                faces_type = strcat(part);
                cornerpts = 4;
                nr_of_quads    = fread(fid,                 1, castingtype,castingendi);
                geofilecontents.(part_descr).(['faces_' faces_type]) = fread(fid,[cornerpts nr_of_quads],castingtype,castingendi)';

            elseif isempty(part)
                [~] = fread(fid,1,'*char');
                if feof(fid)
                    break;
                else
                    error('Unexpected error; did not reach end of file after empty part');
                end

            else
                warning('This statement should not be reached. Most likely the GEO file was saved as multiple timesteps.. ignoring next timesteps.');
                break; % break out of while loop....
            end
        end
        fclose(fid); % close geo file
        
        %%% SOURCE OF GEO FILE FORMAT SPECIFIER BELOW
        %%% source: https://visualization.hpc.mil/wiki/EnSight_Gold_Geometry_Casefile_-_C_Binary_Example
        % All Data is plainly assumed to be exported as C binary Little Endian!
        %
        % C Binary                                            80 chars
        %+% ADDITIONAL_PARAMETER BEGIN TIME STEP                80 chars     % apparently added in fluent12 exported encas/geo files
        % description line 1                                  80 chars
        % description line 2                                  80 chars
        % node id <off/given/assign/ignore>                   80 chars
        % element id <off/given/assign/ignore>                80 chars
        %#% [extents                                            80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% xmin xmax ymin ymax zmin zmax]                       6 floats    % apparently skipped in fluent12 exported encas/ geo files
        
        % part                                                80 chars     % Multiple parts may exist within 1 file; loop over these until end of file
        % #                                                    1 int       % Index of file
        % description line                                    80 chars     % Name of current part
        % coordinates                                         80 chars
        % nn                                                   1 int       % Count of xyz coordinates
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        % x_n1 x_n2 ... x_nn                                  nn floats
        % y_n1 y_n2 ... y_nn                                  nn floats    % Loaded at once into vertices
        % z_n1 z_n2 ... z_nn                                  nn floats
        % element type                                        80 chars     % Determines np (nr of cornerpoints
        % ne                                                   1 int       % Count of faces
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        % n1_e1 n2_e1 ...                                     np_e1
        % n1_e2 n2_e2 ...                                     np_e2
        % .
        % .
        % n1_ne n2_ne ... np_ne                            ne*np ints
        % element type                                        80 chars     % apparently skipped in fluent12 exported encas/ geo files
        % .
        % .
        % part                                                80 chars     % Go BACK to PART above and repeat all elements
        
        %#% .                                                                % apparently skipped in fluent12 exported encas/ geo files
        %#% .                                                                % apparently skipped in fluent12 exported encas/ geo files
        %#% part                                                80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% #                                                    1 int       % apparently skipped in fluent12 exported encas/ geo files
        %#% description line                                    80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% block [iblanked] [with_ghost] [range]               80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1)           3 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [imin imax jmin jmax kmin kmax] # if range used:     6 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% nn = (imax-imin+1)* (jmax-jmin+1)* (kmax-kmin+1)                 % apparently skipped in fluent12 exported encas/ geo files
        %#% ne = (imax-imin)*(jmax-jmin)*(kmax-kmin)                         % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% x_n1 x_n2 ... x_nn                                  nn floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% y_n1 y_n2 ... y_nn                                  nn floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% z_n1 z_n2 ... z_nn                                  nn floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% [ib_n1 ib_n2 ... ib_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [ghost_flags]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [gf_e1 gf_e2 ... gf_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [node_ids]                                          80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [element_ids]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% part                                                80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% #                                                    1 int       % apparently skipped in fluent12 exported encas/ geo files
        %#% description line                                    80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% block rectilinear [iblanked] [with_ghost] [range]   80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1)           3 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [imin imax jmin jmax kmin kmax] # if range used:     6 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% nn = (imax-imin+1)* (jmax-jmin+1)* (kmax-kmin+1)                 % apparently skipped in fluent12 exported encas/ geo files
        %#% ne = (imax-imin)*(jmax-jmin)*(kmax-kmin)                         % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% x_1 x_2 ... x_i                                      i floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% y_1 y_2 ... y_j                                      j floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% z_1 z_2 ... z_k                                      k floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% [ib_n1 ib_n2 ... ib_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [ghost_flags]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [gf_e1 gf_e2 ... gf_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [node_ids]                                          80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [element_ids]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% part                                                80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% #                                                    1 int       % apparently skipped in fluent12 exported encas/ geo files
        %#% description line                                    80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% block uniform [iblanked] [with_ghost] [range]       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1)           3 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [imin imax jmin jmax kmin kmax] # if range used:     6 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% nn = (imax-imin+1)* (jmax-jmin+1)* (kmax-kmin+1)                 % apparently skipped in fluent12 exported encas/ geo files
        %#% ne = (imax-imin)*(jmax-jmin)*(kmax-kmin)                         % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% x_origin y_origin z_origin                           3 floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% x_delta y_delta z_delta                              3 floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% [ib_n1 ib_n2 ... ib_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [ghost_flags]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [gf_e1 gf_e2 ... gf_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [node_ids]                                          80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [element_ids]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
    end                            % load .geo file

    function   scalarfilecontents = loadScalarFile(pathtofile,geometry_struct,istimeresolved_simulation)
        % load scalar file; see below for details on the file format
        version = ver('matlab'); 
        if str2double(version.Version) > 8
            narginchk(3,3); nargoutchk(1,1);
        end
        fid = fopen(pathtofile,'r');
        timepoint = 1;
        
        if istimeresolved_simulation == 1
            [first_line]                   = fread(fid,[1 80],'*char' ); % BEGIN TIME STEP % IGNORE
        end

        scalarfilecontents.description = fread(fid,[1 80],'*char' ); % description line 1          80 chars
        
        while ~feof(fid)
            part                       = fread(fid,[1 80],'*char');  % part                        80 chars
            
            if ~isempty(part) && strcmp(part(1:4),'part')
                
                %%% STEP 1: LOAD ALL DATA FOR THE CURRENT ELEMENT
                part_nr                = fread(fid,     1,castingtype,castingendi );    % #                            1 int
                [~]                    = fread(fid,[1 80],'*char');
                [nrofelements_toread,elementdescr]    = getCountAndNameFromGeo( geometry_struct, part_nr );
                read_data              = fread(fid,nrofelements_toread,'*float',castingendi );
                
                %%% STEP 2: SAVE ALL DATA FOR THE CURRENT ELEMENT
                scalarfilecontents.(elementdescr).(['timestep' num2str(timepoint)]) = read_data;
                scalarfilecontents.(elementdescr).part_nr                           =   part_nr;
                
            elseif ~isempty(part) && strcmp(part(1:13),'END TIME STEP') % Check if NEXT TIME STEP IS STARTING OR IF FILE ENDS
                BEGIN_TIME_STEP = fread(fid,[1 80],'*char');      % Try reading a 80char line
                if feof(fid) % END OF FILE REACHED; STOP READING
                    break;
                elseif strcmp(BEGIN_TIME_STEP,first_line)   % NEW TIME STEP FOUND, restart while loop after reading the duplicate description
                    [~]       = fread(fid,[1 80],'*char');  % description line 1          80 chars % IGNORE description after saving 1st time
                    timepoint = timepoint + 1;              % Increase Time Step with 1
                end
                
            elseif isempty(part)
                [~] = fread(fid,1,'*char');
                if feof(fid)
                    break;
                else
                    error('Unexpected error; did not reach end of file after empty part')
                end
                
            else % code should not arrive here... throw error!
                disp(['lastfile: ' pathtofile ])
                error(['Unexpected error; Last element read: ' num2str(part_nr) ' (timepoint: ' num2str(timepoint) ')'])
                
            end % if
        end % while
        
        scalarfilecontents.timesteps = timepoint; % set number of loaded timepoints
        fclose(fid); % close file
        
        
        %%% SCALAR - from: http://www-vis.lbl.gov/NERSC/Software/ensight/docs/OnlineHelp/UM-C11.pdf
        %+% ADDED: BEGIN TIME STEP    80 chars %%% ADDED
        % description line 1          80 chars
        % part                        80 chars
        % #                            1 int
        % coordinates                 80 chars
        % s_n1 s_n2 ... s_nn          nn floats
        % part                        80 chars
        % .
        % .
        % part                        80 chars
        % #                            1 int
        % block # nn = i*j*k          80 chars
        % s_n1 s_n2 ... s_nn          nn floats
        
    end         % load .scl# file
   