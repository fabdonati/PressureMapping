function [a,hd] = read_image_vtk2(name,bl)
% Updated version by P.Lamata
%
% Version control:
% - 04/02/2014: hack to recover from negative label values (an issue when
% you get the endian wrong. 
% - 20/12/2013: default big/little depending on OS: unix takes 'l', else 'b'
% - 29/06/2012: definition of the default big/little endian as 'b' (changed
% from 'l')
% - (20-11-09) Option of read a binary managing the memory (needed for big datasets)
% memory is introduced if 
% - (04-11-09) Removal of the permutation XY

if isunix()
    DefaultBigOrLittleEndian = 'l'; 
else
    DefaultBigOrLittleEndian = 'b'; 
end
if nargin<2
    bl = DefaultBigOrLittleEndian;
else
    % Map numerical convention into little/big character:
    switch bl
        case 0
            bl = 'l';
        case 1
            bl = 'b';
    end
end

bDebug=1;
hd.File = name;
if (bDebug), fprintf(1,'        Reading image %s\n',name); end

bLoadBigBinary=0;


fid = fopen(name,'rb',bl);

p = 1;
head = do_blanks(fgetl(fid));
hd.head = head;
if head(1) ~= '#'
    fclose(fid);
    error('line 1 header not recognised');
    return
end
contents = fgetl(fid);
hd.contents = contents;
datatype = do_blanks(fgetl(fid));
if ~(strcmp(upper(datatype),'BINARY') | strcmp(upper(datatype),'ASCII'))
    fclose(fid);
    error('data type not recognised')
    return
end
if strcmp(upper(datatype),'BINARY')
    dform = 0;
else
    dform = 1;
end
hd.datatype = datatype;
dataset = do_blanks(fgetl(fid));
if isempty(findstr('DATASET',upper(dataset)))
    fclose(fid);
    error('data set not specified')
    return
end
[H,dataset] = strtok(dataset);
dataset = do_blanks(dataset);
hd.dataset = dataset;
if ~strcmp(dataset,'STRUCTURED_POINTS') & ~strcmp(dataset,'RECTILINEAR_GRID')
    error('not an image');
end
data = '';
fin = 0;
hd.origin = [0 0 0];
hd.spacing = [1 1 1];
formatdata = '';
while fin == 0
    q = ftell(fid);
    b = fgetl(fid);
    if isempty(b)
        continue
    end
    p = length(b);
    w = q+p+1;
    b = do_blanks(b);
    [H,T] = strtok(b);
    switch upper(H)
        case 'DIMENSIONS'
            dim = str2num(T);
            hd.dim = dim;
        case 'ORIGIN'
            origin = str2num(T);
            hd.origin = origin;
        case 'SPACING'
            spacing = str2num(T);
            hd.spacing = spacing;
       case 'X_COORDINATES'
            [H, T] = strtok(T);
            xp = str2num(H);
            if dform == 0
                hd.xc = fread(fid,xp,'float');
            else
                hd.xc = fscanf(fid, '%f ', xp);
            end
        case 'Y_COORDINATES'
            [H, T] = strtok(T);
            yp = str2num(H);
            if dform == 0
                hd.yc = fread(fid,yp,'float');
            else
                hd.yc = fscanf(fid, '%f ', yp);
            end
        case 'Z_COORDINATES'
            [H, T] = strtok(T);
            zp = str2num(H);
            if dform == 0
                hd.zc = fread(fid,zp,'float');
            else
                hd.zc = fscanf(fid, '%f ', zp);
            end
          case 'POINT_DATA'
            points = str2num(T);
            hd.points = points;
            indata = 'point_data';
        case 'CELL_DATA'
            cells = str2num(T);
            hd.cells = cells;
            indata = 'cell_data';
        case 'SCALARS'
            [H,T] = strtok(T);
            hd.scalars_name = H;
            [H,T] = strtok(T);
            formatdata = H;
            hd.format = formatdata;
            hd.comp = str2num(T);
            switch formatdata
                case 'unsigned_short'
                    fd = 'uint16';
                case 'short'
                    fd = 'int16';
                case 'float'
                    fd = 'float';
                case 'unsigned_char'
                    fd = 'uchar';
                case 'char'
                    fd = 'char';
                case 'double'
                    fd = 'double';
                otherwise
                    fprintf(['ERROR! format type ' formatdata ' not recognised']);
                    %fclose(fid);
                    
                    fprintf('Debugging info: formatdata is char = %i (what it should be)\n',ischar(formatdata));
                    fprintf('                formatdata is nume = %i\n',isnumeric(formatdata));
                    fprintf('                formatdata is cell = %i\n',isa(formatdata,'cell'));
                    fprintf('                formatdata is struct = %i\n',isa(formatdata,'struct'));
                    fprintf(' ... attempting to recover making fd = %s\n',formatdata);
                    fd = formatdata;
                    %error(['format type ' formatdata ' not recognised']);
            end
%         case 'COLOR_SCALARS'            
%             [H,T] = strtok(T);
%             hd.scalars_name = H;            
%             formatdata = T;
%             hd.format = formatdata;
%             hd.comp = str2num(T);
%             switch formatdata
%                 case 'scalars'                    
%                     fd = 'uchar';
%             end
        case {'LOOKUP_TABLE','COLOR_SCALARS'}
            if strcmp(H,'COLOR_SCALARS')
                [H,T] = strtok(T);
                hd.scalars_name = H;            
                formatdata = H;
                hd.format = formatdata;
                hd.comp = str2num(T);
                switch formatdata
                    case 'scalars'                    
                        fd = 'uchar';
                    otherwise
                        fprintf(['ERROR! format type ' formatdata ' not recognised']);              
                        fd = 'uchar';
                end
            end
            q = ftell(fid);
            if q ~= w
                fseek(fid,w-q,0);
            end
            if isempty(formatdata)
                error('no data format given')
            end
            if strcmp(indata,'point_data')
                a=false((dim(1)-1)*(dim(2)-1)*(dim(3)-1),1);
                if dform == 0
                    if (bLoadBigBinary)
                        mem=memory;
                        Maxn= floor(mem.MaxPossibleArrayBytes/16);
                        if points<Maxn
                            a = fread(fid,points,fd);
                        else
                            nBlocks = ceil(points/Maxn);
                            fprintf('Reading work splitted in %i blocks\n',nBlocks);
                            for iB=1:nBlocks
                                i1=1+(iB-1)*Maxn;
                                if (iB==nBlocks)
                                    i2=points;
                                else
                                    i2=iB*Maxn;
                                end
                                a(i1:i2) = fread(fid,(i2-i1+1),fd);
                            end
                        end
                    else
                        a = fread(fid,points,fd);
                    end
                else
                    booleanload=0;
                    if (booleanload)
                        p1 = floor(points/2);
                        a1 = boolean(fscanf(fid, '%f', p1));
                        p2 = points-p1;
                        a2 = boolean(fscanf(fid, '%f', p2));
                        a = [a1 a2];
                    else
                        a = fscanf(fid, '%f', points);
                    end
                end
                a = reshape(a,dim(1),dim(2),dim(3));
                %a = permute(a,[2 1 3]);
                hd.dim = size(a);
                fin = 1;
            elseif strcmp(indata,'cell_data')
                a = fread(fid,[cells],fd);
                a = reshape(a,dim(1)-1,dim(2)-1,dim(3)-1);
                %a = permute(a,[2 1 3]);
                hd.dim = size(a);
                fin = 1;
            end
        otherwise
            fclose(fid);
            error(['do not recognise ' b]);
    end
end
if (bDebug), fprintf(1,'        Image read. Size=%i,%i,%i. Sampling=%1.1f,%1.1f,%1.1f\n\n',size(a),hd.spacing ); end
fclose(fid);

% Check if the endian is correct: recover to other option if the values are
% only negative (case of an image with labels):

Labels = unique(a);
nMaxLabels = 100;
if numel(Labels) < nMaxLabels
    % this seems to be an image of labels,
    if numel(Labels) > 1
        % This is not an empty image
        if max(Labels) == 0
            % All labels are negative
            fprintf(' WARNING! All grayvalues (labels) were negative, and this is interpreted to be an error caused by wrong choise of endian while reading VTK file\n')
            switch bl
                case 'l', blnew = 'b'; bReadAgain = 1;
                case 'b', blnew = 'l'; bReadAgain = 1;
                otherwise, bReadAgain = 0; fprintf('   ERROR: not recognised endian: %s\n',bl);
            end
            if(bReadAgain)
                fprintf('   ... attempt to recover by changing from %s to %s endian\n',bl,blnew)
                [a,hd] = read_image_vtk2(name,blnew);
            end 
        end
    end
end
        


function b = do_blanks(a)
%
b = deblank(a);
b = b(end:-1:1);
b = deblank(b);
b = b(end:-1:1);
r = findstr(b,'  ');
while ~isempty(r)
    b = strrep(b,'  ',' ');
    r = findstr(b,'  ');
end

    