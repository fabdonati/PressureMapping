function [PH,hPH] = io_ReadMedicalImage(file,options)
% Function that encapsulates the reading of several medical imaging
% formats, and returns a simplified and unified one adapted to the use of
% the code by Pablo Lamata

bViewIsosurface = 0;

[nameHeart formatHeartFile] = RemoveExtensionFromImageName(file);
switch formatHeartFile
    case {'.bm'}
        ImReader = 'PhilipsImage';
    case {'.nii','.NII','.nii.gz'}
        ImReader = 'nii_reader';
    case {'.vtk','.VTK'}
        ImReader = 'ReadVTK2';
    case {'.par','.rec'}
        ImReader = 'ParRec';
    otherwise
        ImReader = 'ReadData3D';
end

if isunix()
    DefaultBigOrLittleEndian = 'l'; 
else
    DefaultBigOrLittleEndian = 'b'; 
end

bEndingSpecified = 0;
bAvoidEndianIssue = 1;
bPhilipsParseDataOnly = 0;

if nargin>=2
    if isfield(options,'bl'), bl = options.bl; bEndingSpecified = 1; end
    if isfield(options,'bViewIsosurface'), bViewIsosurface = options.bViewIsosurface; end
    if isfield(options,'bAvoidEndianIssue'), bAvoidEndianIssue = options.bAvoidEndianIssue; end
    if isfield(options,'bPhilipsParseDataOnly'), bPhilipsParseDataOnly = options.bPhilipsParseDataOnly; end
end
if exist(file,'file');
    switch ImReader
        case 'ParRec'
            info= par_read_header([nameHeart '.par']) ;
            PH  = par_read_volume(info);
            hPH = ParseHeader(info); 
        case 'PhilipsImage'
            fprintf('Reading raw image %s (assuming Philips format) ...\n',file)
            
            %For Philips conversion, the machine format endian changes as
            %follows
            if isunix()
                optionsPhilips.machineformat = 'b'; 
            else
                optionsPhilips.machineformat = 'l'; 
            end
                
            optionsPhilips.bParseOnly = bPhilipsParseDataOnly;
            
            [PH,hPH] = io_ConvertPhilipsImage(file,optionsPhilips);
            hPH = ParseHeader(hPH);
        case 'nii_reader'
            fprintf('Reading image %s (nifti format)... ',file);
            niiImage = load_untouch_nii(file);
            fprintf(' Finished!\n ');
            PH  = niiImage.img;
            hPH = ParseHeader(niiImage); 
        case 'ReadData3D'
            fprintf('Reading image %s (with ReadData3D)... ',file);
            [PH,info]=ReadData3D(file);
            fprintf(' Finished!\n ');
            if numel(find(PH))==0
                % A hack over a strange bug in ReadData3D with GIPL:
                PH = gipl_read_volume(info);
            end
            hPH = ParseHeader(info); 
        case 'ReadVTK2'
            fprintf('Reading image %s (vtk format)... ',file);
            if(bEndingSpecified)
                fprintf(' ending specified: %s',bl);
            else
                bl = DefaultBigOrLittleEndian;
                fprintf(' ending by default from machine: %s',bl);
            end
            [PH,hPH] = read_image_vtk2(file,bl);
            % Little/big endian issues will cause a distortion of the
            % levels of intensity, and normally introducing negative
            % values. The other possibility is to introduce strange
            % exponential values, very tiny or very big. 
            if(bAvoidEndianIssue)
                values = sort(unique(PH)); 
                Range = values(end) - values(1); 
                % Remove the 0 for the log scale:
                if numel(values)>1
                    Extremes = [values(1) values(end)]; 
                    I = find(Extremes==0);
                    if numel(I)>0, 
                        if I==1,
                            Extremes(1) = values(2); 
                        else
                            Extremes(2) = values(end-1);
                        end
                    end
                else
                    if values == 0,
                        % Avoid the Inf in the log, no need to change
                        % endian
                        Extremes = 1;
                    else
                        Extremes = values;
                    end
                end
                LogRange = sum(abs(log(Extremes)));
                if min(PH(:)) < 0 || Range > 1e30 || LogRange > 50
                    fprintf('WARNING! Reading vtk file led to some negative gray levels (labels).\n');
                    fprintf('         Attempt to recover it by changing little/big endinan.\n');
                    switch bl
                        case 'b', options.bl = 'l';
                        case 'l', options.bl = 'b';
                    end
                    options.bAvoidEndianIssue = 0;
                    [PH,hPH] = io_ReadMedicalImage(file,options);
                    fprintf(' labels available:\n'); fprintf('%i, ',unique(PH(:))); fprintf('\n');
                    return;
                end
            end
            fprintf(' Finished!\n ');
            hPH = ParseHeader(hPH);
    end
else
    fprintf('ERROR! image file does not exist: %s\n',file);
    PH = NaN;
    hPH = NaN;
end
if(bViewIsosurface)
    show_segment_surface(PH,hPH.Mv2w);
end