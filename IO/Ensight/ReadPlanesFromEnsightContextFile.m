function [planes] = ReadPlanesFromEnsightContextFile(filename)
% Function to read and extract the location and extents of the planes drawn
% in an Ensight volume, an information which is stored in a context file in
% this way:
%
% part: select_byname_begin
% "(CASE:CMDKEYWORD_CURCASE)Total Volume, Mins/Maxs=(I:1/224,J:1/168,K:1/28)"
% part: select_byname_end
% clip: begin
% part: description Clip_grid3
% clip: type grid
% clip: domain intersect
% clip: grid_pts 40 40
% clip: extents finite
% clip: tool plane
% clip: delta 0.000000e+000 0.000000e+000 0.000000e+000
% clip: plane 1 -5.294734e+000 -2.170538e+001 1.112376e+002
% clip: plane 2 -1.395352e+001 -4.235797e+001 1.443811e+002
% clip: plane 3 1.763874e+001 -6.601831e+001 1.378914e+002
% clip: end
% clip: create
% part: select_byname_begin
% "(CASE:CMDKEYWORD_CURCASE)Total Volume, Mins/Maxs=(I:1/224,J:1/168,K:1/28)"
% part: select_byname_end

    bDebug=0;
    iPlane = 0;
    planes = [];
    fid = fopen(filename,'r');
    if fid==-1
        fprintf('ERROR! File %s could not be opened!\n',filename);
    else
        fprintf('Reading ensight context file %s ...',filename);
        l = 0;
        while(~feof(fid))
            [tline,l] = ReadOneLine(fid,l);
            if(bDebug),fprintf('Line %i processed: %s\n',l,tline);end
            if l==470
               a=1; 
            end
            if(numel(tline)>10)
                if strcmp(tline(1:11),'clip: begin')
                    if(bDebug), fprintf('New plane candidate from line %i, %s\n',l,tline); end
                    % In the next line we have the description
                    %[tline,l] = ReadOneLine(fid,l);
                    [plane,bValidPlane,l] = GetPlaneDefinition(fid,l);
                    if (bValidPlane)
                       iPlane = iPlane+1;
                       planes(iPlane).plane = plane;
                       if(bDebug), fprintf('New plane added (line %i): %s\n',l,planes(iPlane).plane.Description); end
                    end
                end
            end
        end        
        fprintf('\n %i = number of planes read\n',iPlane);
    end
end
function [plane,bValidPlane,l] = GetPlaneDefinition(fid,l)
    bDebug=0;
    bEnd = 0;
    bValidPlane=0;
    Description= '';
    Type = '';
    Delta = NaN;
    Plane = NaN*ones(3,3);
    while(~bEnd)
        [tline,l] = ReadOneLine(fid,l);
        switch tline(1:5)
            case 'part:'
                if strcmp(deblank2(tline(6:17)),'description')
                    Description = deblank2(tline(18:end));
                else
                    if(bDebug),fprintf('Warning! a "part:" piece of information not recognised, line %i: %s!\n',l,tline); end
                end
            case 'clip:'
                switch deblank2(tline(6:9))
                    case 'typ'
                        Type = deblank2(tline(11:end));
                        if strcmp(Type,'grid'), 
                            bValidPlane=1;
                        else                    bValidPlane=0;  end
                    case 'end'
                        bEnd = 1;
                    case 'del'
                        Delta = sscanf(tline(12:end),'%f',3);
                    case 'pla'
                        PlaneNumber = sscanf(tline(12:end),'%i',1);
                        Plane(PlaneNumber,1:3) = sscanf(tline(15:end),'%f',3);                    
                end
            otherwise
                if(bDebug),fprintf('Warning! a piece of information not recognised, line %i: %s!\n',l,tline); end
        end
    end
    plane.Description = Description;
    plane.Type  = Type;
    plane.Delta = Delta;
    plane.Plane = Plane;
end
function [tline,l] = ReadOneLine(fid,l)
    tline = deblank2(fgetl(fid)); 
    l=l+1;
end