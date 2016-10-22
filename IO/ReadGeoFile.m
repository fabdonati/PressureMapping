function [Imageheader,points,bError] = ReadGeoFile(directory,GEOfileName,bReadPoints,bBinaryFormat)
% Function to read the EnSight GEO file, and retrieve the necessary
% information to build an image: the origin, spacing and orientation.
% TODOs: 
% - get the correct image size from the "truncated extents" now exported
% accordingly to format Ensight 10.0
%
% Example header:
%  OpenCMISS Exported Encas Model Geometry File
%  EnSight 9.0.3
%  node id assign
%  element id assign
%  extents
%      0.00000    95.45460
%      0.00000    17.89770
%      0.00000    18.00003
%  part
%          1
% Total Volume, Mins/Maxs=(I:1/97,J:1/19,K:1/19)
%  block
%         97        19        19
%
%--------------------------------------------------------------------------
%
% EnSight Model Geometry File
% EnSight 10.0.3
% node id assign
% element id assign
% extents
% -3.33749e+01 4.36251e+01
% -1.75902e+02 8.24317e+01
% -1.51235e+02 1.67099e+02
% part
%          1
% Total Volume, Mins/Maxs=(I:1/192,J:1/156,K:1/40)
% block range
%        192       156        40
%          1       192         1       156         5        40

% Parameters:
EnsightUnits = 'mm'; 
ImageUnits   = 'mm';

bDebug=0;
bFindOrigin = 0;

Imageheader = [];
points = NaN;
bError = 0;
%--------------------------------------------------------------------------
    if nargin<3
        bReadPoints=0;
    end
    if nargin<4
        bBinaryFormat=0;
    end
    if ~exist(fullfile(directory, GEOfileName),'file')
        % Check if the ending is not correct:
        if ~strcmp(GEOfileName(end-3:end),'.geo')
            fprintf(' ... adding .geo to %s\n',GEOfileName);
            GEOfileName = [GEOfileName '.geo'];            
        end
        if ~exist(fullfile(directory, GEOfileName),'file')
            % Check if a novel ending is not correct:
            GEOfileName = [GEOfileName '0000'];
        end
    end
    filename = fullfile(directory, GEOfileName);
    fid = fopen(filename,'r');
    [a b c d] = fopen(fid);
    if fid==-1
        fprintf('ERROR in ReadGeoFile! GEO file could not be opened: %s\n',filename);
        bError = 1;
        return;
    end
    % Check if this is a binary or an ASCII format:        
    geofilecontents.type   = fread(fid,[1 80],'*char'); %
    if strcmp(deblank2(geofilecontents.type),'C Binary') || bBinaryFormat
        bBinaryFormat = 1;
        % This is a binary!
        fclose(fid);
        GeoData = loadGeoFile(filename,0);
    else
        fclose(fid);
        fid = fopen(filename,'r');
        
        Imageheader.Description = 'Header constructed from an Ensight GEO file';
        % Units in the PPE workflow use m/s!
        Imageheader.VelocityUnits = 'm/s';

        bExtentsFound = 0;    
        bBlockFound = 0;
        % When reading the ensight values of points, get rid of one line
        % after the keyword "block"
        DiscardLines = 1;

        % Get the extents of the volume in physical units:
        bExtentsFound = findWordInFile(fid,'extents',filename);    
        if(bExtentsFound)
            Extent = fscanf(fid,'%f',6);
            Extents(1:3,1) = Extent(1:2:5);
            Extents(1:3,2) = Extent(2:2:6);
        else
            fprintf('ERROR! No extents found in file %s\n',filename);
            Imageheader = NaN;
            bError =1;
            return;
        end
        % Get the size of the volume in number of samples (voxels):
        [bBlockFound, tline] = findWordInFile(fid,'block',filename);  
        if(bBlockFound)
            ImSize = fscanf(fid,'%i',3);
            spacing = Extents(1:3,2)./ImSize;
            % The format 10 of Ensight has also the block range:
            if strcmp(deblank2(tline),'block range')
                % When reading the ensight values of points, get rid of
                % two lines after the keyword "block"
                DiscardLines = 2;
                ImExtents = fscanf(fid,'%i',6);
                ImSize(1) = ImExtents(2) - ImExtents(1) + 1;
                ImSize(2) = ImExtents(4) - ImExtents(3) + 1;
                ImSize(3) = ImExtents(6) - ImExtents(5) + 1; 
            end                    
        else
            fprintf('ERROR! No block found in file %s\n',filename);
            Imageheader = NaN;
            bError = 1;
            return;
        end
        bSizeCorrect = CheckSize(Extents(:,2),EnsightUnits);
        if(~bSizeCorrect)
            fprintf('WARNING in ReadGeoFile! The extents of the Ensight file (Max: %f,%f,%f) are not as expected (units: %s)\n',Extents(:,2),EnsightUnits);
        end
        Imageheader.dim = ImSize;
        Imageheader.Extents = Extents;
        Imageheader.origin  = Extents(1:3,1);
        Imageheader.spacing = spacing;
%             if exist('ImExtents','var')
%                 % it may be possible that the origin is shifted:
%                 for iC = 1:3
%                     if ImExtents(2*iC -1) > 1
%                         Offset = (ImExtents(2*iC -1) -1) * Imageheader.spacing(iC);
%                         Imageheader.origin(iC) = Imageheader.origin(iC) + Offset;
%                         Imageheader.Extents(iC,1) = Imageheader.Extents(2*iC -1) + Offset;
%                     end
%                 end
%             end
        fclose(fid);
    end

    if(~bBinaryFormat)
        dimensionOfValues = 3;
        if (bReadPoints)
            optionsRead.dimensionOfValues =dimensionOfValues;
            optionsRead.DiscardLines = DiscardLines;
            [points,status] = ReadEnsightValues(directory,GEOfileName,Imageheader,optionsRead);
            if(status==-1);
                fprintf('ERROR in ReadGeoFile! GEO file looks corrupted, no points read!: %s\n',filename);
                clear points;
                points = NaN;
                RotationMatrix = eye(3);
                bError = 1;
                return;
            end
            Xcoor = points(:,:,:,1); Xcoor = Xcoor(:);
            Ycoor = points(:,:,:,2); Ycoor = Ycoor(:);
            Zcoor = points(:,:,:,3); Zcoor = Zcoor(:);
            D1 = squeeze(points(2,1,1,:) - points(1,1,1,:))';
            D2 = squeeze(points(1,2,1,:) - points(1,1,1,:))';
            D3 = squeeze(points(1,1,2,:) - points(1,1,1,:))';
            spacing(1) = sqrt(D1(1)^2 + D1(2)^2 + D1(3)^2);
            spacing(2) = sqrt(D2(1)^2 + D2(2)^2 + D2(3)^2);
            spacing(3) = sqrt(D3(1)^2 + D3(2)^2 + D3(3)^2);
        end
    else
        % Parse the output gathered from the binary reader:
        points = GeoData.points;
        Xcoor = points(:,1);
        Ycoor = points(:,2);
        Zcoor = points(:,3);
        Extents(1:3,1) = GeoData.extents(1:2:5);
        Extents(1:3,2) = GeoData.extents(2:2:6);
        [D1,D2,D3,spacing] = FindUnitaryVectors(points,bDebug);
        ImSize = GeoData.dimensions;

    end
    % Check the limits of the Bounding Box:
    xMin = min(Xcoor);
    yMin = min(Ycoor);
    zMin = min(Zcoor);
    xMax = max(Xcoor);
    yMax = max(Ycoor);
    zMax = max(Zcoor);
    Extents2 = [ xMin yMin zMin; xMax yMax zMax]';
    if bBinaryFormat
        % Not managed to cast the type in Extents, not possible to check
        % this, simply take Extents2 as valid
        Extents = Extents2;
        Imageheader.Extents = Extents;
        Imageheader.dim = ImSize;
    else
        for iE=1:3
            for jE=1:2
                A =Extents(iE,jE);
                B =Extents2(iE,jE); 
                if A~=B
                    fprintf('ERROR in ReadGeoFile, extents do not match! (%f and %f), coordinate %i.%i\n',A,B,iE,jE);
                    bError = 1;
                end
            end
        end
    end

    nD1 = VectorNormalization(D1);
    nD2 = VectorNormalization(D2);
    nD3 = VectorNormalization(D3);
    RotationMatrix = [nD1;nD2;nD3];
    if min(spacing)==0
        fprintf('ERROR in ReadGeoFile! Spacing calculated from points results in zero!\n');
        fprintf('MOST LIKELY an error during generation of the GEO file\n');
        fprintf('  ... PLEASE DO CHECK THE CONVERSION SETTINGS IN STEP 3!\n');
        bError = 1;
    else
        Imageheader.spacing = spacing;
    end
    % Finally, return units as specified
    if strcmp(EnsightUnits,ImageUnits)
        Scale = 1;
    else
        switch EnsightUnits
            case 'm'
                switch ImageUnits
                    case 'mm'
                        Scale = 1000;
                    case 'cm'
                        Scale = 100;
                    case 'dm'
                        Scale = 10;
                end
            case 'mm'
                switch ImageUnits
                    case 'cm'
                        Scale = 0.1;
                    case 'dm'
                        Scale = 0.01;
                    case 'm'
                        Scale = 0.001;
                end
            otherwise
                fprintf('ERROR!  Scaling option not coded!\n');
        end
    end
    if(bFindOrigin)
        % Find the origin: one arbitrary corner:
        Coords = points(:,:,:,2);
        Coords = Coords(:);
        Icandidates = find(Coords == yMin);
        Coords = points(:,:,:,1); Coords = Coords(:); Xcoor = Coords(Icandidates);
        Coords = points(:,:,:,3); Coords = Coords(:); Zcoor = Coords(Icandidates);
        MinZ = min(Zcoor);
        Icand2 = find(Zcoor==MinZ);
        X = Xcoor(Icand2(1));
        origin = [X,yMin,MinZ]';
    else
        origin = Imageheader.Extents(:,1);
    end
    
    Imageheader.Extents = Imageheader.Extents*Scale;
    Imageheader.origin  = origin*Scale;    
    Imageheader.spacing = Imageheader.spacing*Scale; 
    Imageheader.TransformMatrix = RotationMatrix;
    if(bDebug)&&(bReadPoints)
        figure
        % plot bounding box:
        x = [xMin xMin xMin xMax xMax xMax xMin xMin xMin];
        y = [yMin yMin yMax yMax yMax yMin yMin yMax yMax];
        z = [zMin zMax zMax zMax zMin zMin zMin zMin zMax];
        plot3(x,y,z);
        hold on
        f = squeeze(points(ImSize(1),ImSize(2),ImSize(3),:));
        o = squeeze(points(1,1,1,:));
        plot3(o(1),o(2),o(3),'ro');
        plot3(f(1),f(2),f(3),'ro');
        plotline(o,D1,'r');
        plotline(o,D2,'g');
        plotline(o,D3,'b');        
        plot3(points(:,1,1,1),points(:,1,1,2),points(:,1,1,3),'g.');
        plot3(points(:,end,1,1),points(:,end,1,2),points(:,end,1,3),'g.');
        plot3(points(:,1,end,1),points(:,1,end,2),points(:,1,end,3),'g.');
        plot3(points(:,end,end,1),points(:,end,end,2),points(:,end,end,3),'g.');
        plot3(points(1,:,1,1),points(1,:,1,2),points(1,:,1,3),'k.');
        plot3(points(end,:,1,1),points(end,:,1,2),points(end,:,1,3),'k.');
        plot3(points(1,:,end,1),points(1,:,end,2),points(1,:,end,3),'k.');
        plot3(points(end,:,end,1),points(end,:,end,2),points(end,:,end,3),'k.');
        plot3(squeeze(points(1,1,:,1)),squeeze(points(1,1,:,2)),squeeze(points(1,1,:,3)),'c.')
        plot3(squeeze(points(1,end,:,1)),squeeze(points(1,end,:,2)),squeeze(points(1,end,:,3)),'c.')
        plot3(squeeze(points(end,1,:,1)),squeeze(points(end,1,:,2)),squeeze(points(end,1,:,3)),'c.')
        plot3(squeeze(points(end,end,:,1)),squeeze(points(end,end,:,2)),squeeze(points(end,end,:,3)),'c.')
        plot3(origin(1),origin(2),origin(3),'ro','LineWidth',3);
        axis equal;
    end
    points = points*Scale;  
    bSizeCorrect = CheckSize(Imageheader.Extents(:,2),ImageUnits);
    if(~bSizeCorrect)
        fprintf('WARNING! The extents of the resulting image (Max: %f,%f,%f) are not as expected (units: %s)\n',Imageheader.Extents(:,2),ImageUnits);
    end
    if isstruct(Imageheader)
        Imageheader = ParseHeader(Imageheader);
    end
end

function plotline(origin,vector,colour)
    scaling = 5;
    x = [origin(1) origin(1)+(vector(1)*scaling)];
    y = [origin(2) origin(2)+(vector(2)*scaling)];
    z = [origin(3) origin(3)+(vector(3)*scaling)];
    plot3(x,y,z,[colour '-']);
end

function [bSizeCorrect] = CheckSize(data,Units)
    bSizeCorrect=1;
    MeanData = mean(data(:));
    % Check the magnitud order:
    MinThreshold = 50;
    MaxThreshold = 1500;
    switch Units
        case 'mm', scale = 1;
        case 'cm', scale = 0.1;
        case 'dm', scale = 0.01;
        case 'm', scale = 0.001;
    end
    if MeanData<MinThreshold*scale
        bSizeCorrect=0;
        fprintf('\n +++ WARNING! +++ Data looks too small!\n\n');
    end
    if MeanData>MaxThreshold*scale
        bSizeCorrect=0;
        fprintf(' WARNING! Data looks too big!\n');
    end
end

function [Vout] =VectorNormalization(Vin)
    Vout = Vin / sqrt(sum(Vin.^2));
end
function [n] = normalise(v)
    n = v ./ sqrt(sum(v.^2));
end
function [] = plotvector(P,D,s,color)
    plot3([P(1) P(1)+s*D(1)],[P(2) P(2)+s*D(2)],[P(3) P(3)+s*D(3)],color)
end
function         [D1,D2,D3,spacing] = FindUnitaryVectors(points,bDebug)
    % Estimate the 3D orientation, first with the two first points in
    % the list:
    P0 = points(1,:);
    D1 = (points(2,:) - P0);
    spacing(1) = sqrt(sum(D1.^2));
    D1 = normalise(D1);
    s = 20;
    % Now detect the first point out of the first line:
    bSecondLineFound = 0;
    iP = 2;
    epsilon = 0.01;
    if bDebug, figure('color',[1 1 1]); hold on; end
    while ~bSecondLineFound && iP <= size(points,1)
        iP = iP+ 1;             
        P2 = points(iP,:);
        NewD = P2 - P0;
        if abs(1-abs(dot(D1,normalise(NewD)))) > epsilon
            bSecondLineFound = 1;
            if bDebug, plot3(P2(1),P2(2),P2(3),'r*'); end               
        else
%             % Refine D1 (avoid rounding error of a too close sample)
%             D1 = NewD;
%             spacing(1) = sqrt(sum(D1.^2))/(iP-1);
%             D1 = normalise(D1);
            if bDebug, plot3(P2(1),P2(2),P2(3),'.'); end               
        end
    end
    nPerLine = iP+5;
    if bDebug, 
        s1 = iP;
        plotvector(P0,D1,s1,'r');
    end
    D2 = NewD;
    % Check if they are orthogonal
    if(abs(dot(D1,D2))>0.01)
        t = cross(D1,D2);
        D2 = cross(D1,-t);        
    end
    D2 = normalise(D2);
    spacing(2) = dot(NewD,D2);
    if bDebug, plotvector(P0,D2,s,'g'); end
    % Now detect the first point out of the first plane:
    bSecondPlaneFound = 0;
    epsilon2 = 0.000001;
    while ~bSecondPlaneFound && iP <= size(points,1)
        iP = iP + nPerLine;
        P2 = points(iP,:);
        D3temp = P2 - P0;
        NewD = normalise(D3temp);
        if abs(dot(NewD,normalise(D1))^2 + dot(NewD,normalise(D2))^2 - 1) > epsilon2
            bSecondPlaneFound = 1;
            if bDebug, plot3(P2(1),P2(2),P2(3),'r*'); end               
        else
            if bDebug, plot3(P2(1),P2(2),P2(3),'.'); end               
        end
    end
    % Again, find the orthogonal:
    D3 = cross(D1,D2);
    if dot(D3temp,D3)<1
        D3 = -D3;
    end
    D3 = normalise(D3);
    spacing(3) = dot(D3temp,D3);
    if bDebug, plotvector(P0,D3,s,'b'); end
end