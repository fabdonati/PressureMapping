function parrec2Metafile(parFileName,gatingType)
%input: filename (*.PAR or *.par) with extension! (Philips format)
%output: rootfilename.mhd and rootfilename.raw (metafile format in the
%correct dicom space)

%OPTIONAL input: gatingType; values: 'retrospective' (default) or 'prospective'

%check existence of input files and get rid of the extension if it was provided:
% if (strcmp(rootfilename(end-3:end),'.PAR') || ...
%         strcmp(rootfilename(end-3:end),'.par') || ...
%         strcmp(rootfilename(end-3:end),'.rec') || ...
%         strcmp(rootfilename(end-3:end),'.REC'))
%     rootfilename = rootfilename(1:end-4);
% end

% if exist([rootfilename '.PAR'],'file')
%     inputPAR = [rootfilename '.PAR'];
% elseif  exist([rootfilename '.par'],'file')
%     inputPAR = [rootfilename '.par'];
% else
%     error('Input PAR file does not exist!');
% end

if ~exist(parFileName,'file')
    error('PAR file name on the input does not exist!');
end
rootfilename = parFileName(1:end-4);

%load header and in particular the image orientation and position info from PAR file:
hdr = par_read_header(parFileName);
inputREC = hdr.FilenameREC;
if ~exist(inputREC,'file')
    error(['The REC file with name ''' inputREC ''' does not exist. It may be a problem of capital letters?']);
end


%the info about dataType is taken from PAR/REC reader
switch(hdr.BitDepth)
    case 8
        dataType='char';
    case 16
        dataType='short';
    case 32
        dataType='float';
    case 64
        dataType='double';
end

%read the binary data
data = par_read_volume(hdr);
%extract image position and orientation for 3D data (for 3D+t we need to add later the time components):
[TransformMatrix, Offset] = extractparrecinformationfromfile(parFileName);

%For 3D+t:
if (length(hdr.Dimensions) == 4)
    if ~exist('gatingType','var')
        warning('''gatingType'' is set to ''retrospective''');
        gatingType = 'retrospective';
    end
    
    TransformMatrix3D =  TransformMatrix;
    Offset3D = Offset;
    TransformMatrix = eye(4);
    TransformMatrix(1:3,1:3) = TransformMatrix3D;
    Offset = zeros(1,4);
    Offset(1:3) = Offset3D;
    
    %get trigger times (+ time offset in prospective gating)
    triggerTime = zeros(hdr.Dimensions(3:4)); %triggerTime(1,:) is the first slice, all time steps
    for j = 1 : hdr.Dimensions(3)  %loop over slices
        for k = 1 : hdr.Dimensions(4)  %loop over time Frames
            %FOR 2D CINE parrec, where each whole slice is written
%            triggerTime(j,k) = ...
%            hdr.SliceInformation(k+(j-1)*hdr.Dimensions(4)).TriggerTime;
%         FOR truly 3D data
            triggerTime(j,k) = ...
            hdr.SliceInformation(j+(k-1)*hdr.Dimensions(3)).TriggerTime;
            
        end
    end
%    triggerTime = triggerTime / 1000.0; %milliseconds --> seconds

    %number of data sets in a single PAR-REC file; it is 2 typically if the
    %PAR-REC contains magnitude and phase image.
    noOfDataComp = 1;
    if (mod(hdr.Dimensions(4),2) == 0)
        if (triggerTime(1,1) == triggerTime(1,(hdr.Dimensions(4)/2) + 1))
            warning('Input 3D+t image seems to consist of two type of data, e.g. magnitude and image');
            noOfDataComp = 2;
        end
    end
    nTimeFrames = hdr.Dimensions(4) / noOfDataComp;
    
    timeStep = zeros(hdr.Dimensions(3),1);  %time step for each slice (has sense only in retrospective triggering, it should be constant in the prospective one)
    for j = 1 : hdr.Dimensions(3)
        timeStep(j) = (triggerTime(j,nTimeFrames) - triggerTime(j,1)) / (nTimeFrames - 1);
    end
    meanTimeStep = mean(timeStep);
    
    %For the prospective trigger, the mid-point scheme is applied:
    Offset(4) = triggerTime(1,1);  %timeOffset, should be 0 for retrospective trigger
    if strcmp(gatingType,'prospective')
        warning('For the prospective trigger, the mid-point time scheme is applied');
        Offset(4) = Offset(4) + (meanTimeStep / 2.0);
    end
    
    %SAVE THE MAT and TXT file about the timings
    timeOffset = Offset(4)
    save([rootfilename '_timing.mat'],'triggerTime','timeStep','meanTimeStep','timeOffset');
    disp(['Saving output for Imperial tracking into file ' rootfilename '_timing.txt']);
    f = fopen([rootfilename '_timing.txt'],'wt');
    
    for i = 1 : nTimeFrames
        if (i == 1)
            t = timeOffset;
        else
            t = t + meanTimeStep;
        end
        fprintf(f,'%f  ',t);
    end
    fprintf(f,'\n');
    fclose(f);
    
    meanTimeStep = meanTimeStep / 1000.0; %   in seconds for MHD header
    timeOffset = timeOffset / 1000.0;
    Offset(4) = Offset(4) / 1000.0;
    
    if (noOfDataComp == 2)
        Dimensions = hdr.Dimensions;
        Dimensions(4) = hdr.Dimensions(4) / 2;
        %first component of the data set:
        WriteMhaFile([rootfilename '_1.mhd'], Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
        WriteBinary([rootfilename '_1.raw'],data(:,:,:,1:(hdr.Dimensions(4)/2)),dataType);
        %second component of the data set (non-scaled, original values as in the PAR/REC):
%        WriteMhaFile([rootfilename '_2.mhd'], Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
%        WriteBinary([rootfilename '_2.raw'],data(:,:,:,(hdr.Dimensions(4)/2)+1:end),dataType);
        warning(['Assuming that the 2nd component image is the phase image, it will be scaled to cm / sec and written as ' rootfilename '_2scaled.raw']);
        %V = pixel value in REC file, FP = floating point value, DV = displayed value on console
        %RS = rescale slope,           RI = rescale intercept,    SS = scale slope
        %DV = PV * RS + RI
        data2Scaled = double(data(:,:,:,(hdr.Dimensions(4)/2)+1:end));
        data2Scaled = (data2Scaled*hdr.SliceInformation(end).RescaleSlope + hdr.SliceInformation(end).RescaleIntercept) * sum(hdr.PhaseEncodingVelocity) / (1000.0*pi);
        %simple scaling not using PARREC things except for VENC:
        %DEBUG scaling "Radek":
%         data2ScaledRadek = double(data(:,:,:,(hdr.Dimensions(4)/2)+1:end));
%         data2ScaledRadek = (data2ScaledRadek - (4096/2)) / (4096/2) * sum(hdr.PhaseEncodingVelocity);
         WriteMhaFile([rootfilename '_2scaled.mhd'], Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
         WriteBinary([rootfilename '_2scaled.raw'],data2Scaled,dataType);
%         WriteMhaFile([rootfilename '_2scaledRadek.mhd'], Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
%         WriteBinary([rootfilename '_2scaledRadek.raw'],data2ScaledRadek,dataType);
        return;
    else
        WriteMhaFile([rootfilename '.mhd'], hdr.Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
        WriteBinary([rootfilename '.raw'],data,dataType);
        return;
    end
end

%writing of output for 3D image (without time) or for data set with a
%single component
WriteMhaFile([rootfilename '.mhd'], hdr.Dimensions, hdr.Scales, dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
WriteBinary([rootfilename '.raw'],data,dataType);

end %of main

function WriteBinary(filename,data,dataType)
f = fopen(filename,'wb');
fwrite(f,data,dataType);
fclose(f);
end


%FUNCTION par_read_header
function info =par_read_header(filename)
% Function for reading the header of a Philips Par / Rec  MR V4.* file 
%
% info  = par_read_header(filename);
%
% examples:
% 1,  info=par_read_header()
% 2,  info=par_read_header('volume.par');

if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.par', 'Read par-file');
    filename = [pathname filename];
end

fid=fopen(filename,'rb');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end
info.Filename=filename;


mode = -1; nHC=0; nIC=0; nSC=0;
while(true)
    str=fgetl(fid);
    if ~ischar(str), break, end
    if(isempty(str)), continue, end
    
    if(strfind(str,'= DATA DESCRIPTION FILE =')), mode=0; end
    if(strfind(str,'= GENERAL INFORMATION =')), mode=1; end
    if(strfind(str,'= PIXEL VALUES =')), mode=2; end
    if(strfind(str,'= IMAGE INFORMATION DEFINITION =')), 
        mode=3; fgetl(fid);  str=fgetl(fid); 
        % Skip a line
    end
    if(strfind(str,'= IMAGE INFORMATION =')), mode=4; end
    
    if(strfind(str,'= END OF DATA DESCRIPTION FILE =')), mode=5; end
    
    switch(mode)
        case -1;
        case 0
            nHC=nHC+1; HeaderComment{nHC}=str;
        case 1
            if(str(1)=='.')
                [type data]=General_Information_Line(str);
                switch(type)
                    case 'PatientName'
                        info.(type)=data;
                    case 'ProtocolName'
                        info.(type)=data;
                    case 'ExaminationName'
                        info.(type)=data;
                    case 'ExaminationDateTime'
                        info.(type)=data;
                    case 'SeriesType'
                        info.(type)=data;
                    case 'AcquisitionNr'
                        info.(type)=sscanf(data, '%d')';
                    case 'ReconstructionNr'
                        info.(type)=sscanf(data, '%d')';
                    case 'ScanDuration'
                        info.(type)=sscanf(data, '%lf')';
                    case 'MaxNumberOfCardiacPhases'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfEchoes'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfSlicesLocations'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfDynamics'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfMixes'
                        info.(type)=sscanf(data, '%d')';
                    case 'PatientPosition'
                        info.(type)=data;
                    case 'PreparationDirection'
                        info.(type)=data;
                    case 'Technique'
                        info.(type)=data;
                    case 'ScanResolution'
                        info.(type)=sscanf(data, '%d')';
                    case 'ScanMode'
                        info.(type)=data;
                    case 'RepetitionTime'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Fov'
                        info.(type)=sscanf(data, '%lf')';
                    case 'WaterFatShift'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Angulation'
                        info.(type)=sscanf(data, '%lf')';
                    case 'OffCentre'
                        info.(type)=sscanf(data, '%lf')';
                    case 'FlowCompensation'
                        info.(type)=sscanf(data, '%d')';
                    case 'Presaturation'
                        info.(type)=sscanf(data, '%d')';
                    case 'PhaseEncodingVelocity'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Mtc'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Spir'
                        info.(type)=sscanf(data, '%lf')';
                    case 'EpiFactor'
                        info.(type)=sscanf(data, '%lf')';
                    case 'DynamicScan'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Diffusion'
                        info.(type)=sscanf(data, '%lf')';
                    case 'DiffusionEchoTime'
                        info.(type)=sscanf(data, '%lf')';
                    case 'MaxNumberOfDiffusionValues'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfGradientOrients'
                        info.(type)=sscanf(data, '%d')';
                    case 'NumberOfLabelTypes'
                        info.(type)=sscanf(data, '%d')';
                    case 'HeaderComment'
                    otherwise
                        display(' did not recognize the information type: ');
                        display(type);
                        if (size(type) == 0)
                            display('empty, skipping line');
                        else
                            info.(type)=data;
                        end
                end
            end
        case 2
        case 3
            if(str(1)=='#');
                [type datatype datalength]=Image_Information_Line(str);
                if(~isempty(type))
                    nIC=nIC+1;
                    ImageInformationTags(nIC).Name=type;
                    ImageInformationTags(nIC).DataType=datatype;
                    ImageInformationTags(nIC).NumberOfValues=datalength;
                end
            end
        case 4
            if(str(1)~='#');
                nSC=nSC+1;
                vals=regexp(str, '\s+','split');
                vald=sscanf(str, '%lf')';
                current_loc=0;
                for i=1:length(ImageInformationTags)
                    IIT=ImageInformationTags(i);
                    if(strcmp(IIT.DataType,'string'))
                        SliceInformation(nSC).(IIT.Name)=vals{current_loc+1};
                    else
                        SliceInformation(nSC).(IIT.Name)=vald(current_loc+1:current_loc+IIT.NumberOfValues);
                    end
                    current_loc=current_loc+IIT.NumberOfValues;
                end
            end
            
        case 5
        otherwise
            %disp(str);
    end
end
fclose(fid);
info.HeaderComment=HeaderComment;
info.SliceInformation=SliceInformation;
info.ImageInformationTags=ImageInformationTags;

% Add Dimensions and Voxel Spacing. Warning, based on only 1 slice!
infof=info.SliceInformation(1);
if(isfield(infof,'ReconResolution'))
    if(isfield(info,'MaxNumberOfSlicesLocations'))
        zs(1)=info.MaxNumberOfSlicesLocations;
        zs(2)=length(SliceInformation)/zs(1);
        if((mod(zs(2),1)>0)||zs(2)==1)
            zs=length(SliceInformation);
        end
    else
        zs=length(SliceInformation);
    end
    
    info.Dimensions=[infof.ReconResolution zs];
else
    info.Dimensions=[info.ScanResolution length(SliceInformation)];
end
if(isfield(infof,'PixelSpacing'))
    if(isfield(infof,'SliceThickness')&&isfield(infof,'SliceGap'))
        zs=infof.SliceThickness+infof.SliceGap;
    else
        zs=0;
    end
    info.Scales=[infof.PixelSpacing zs];
else
    info.Scales=[0 0 0];
end


[folder,filen]=fileparts(info.Filename);
info.FilenameREC=fullfile(folder,[filen '.rec']);
   
    
% Add bith depth
if(infof.ImagePixelSize)
    info.BitDepth=infof.ImagePixelSize;
else
    if(exist(info.FilenameREC,'file'))
        file_info=dir(info.FilenameREC);
        bytes=file_info.bytes;
        info.BitDepth=(bytes/prod(info.Dimensions))*8;
    end
end
end


function [type datatype datalength]=Image_Information_Line(str)
s=find(str=='(',1,'last');
if(isempty(s)), s=length(str); end
type=str(1:s-1);  data=str(s+1:end);
type=regexp(type, '\s+|/|_', 'split');
type_clean='';
for i=1:length(type)
    part=type{i};
    part(part=='#')=[]; part(part==' ')=[]; partu=uint8(part);
    if(~isempty(part))
        check=((partu>=97)&(partu<=122))|((partu>=65)&(partu<=90));
        if(check)
            part=lower(part); part(1)=upper(part(1));
            type_clean=[type_clean part];
        else
            break;
        end
    end
end
type=type_clean;
while(~isempty(data)&&data(1)==' '), data=data(2:end); end
while(~isempty(data)&&data(end)==' '), data=data(1:end-1); end
if(~isempty(data))
    data=data(1:end-1);
    s=find(data=='*',1,'first');
    if(isempty(s)),
        datalength=1; datatype=data;
    else
        datalength=str2double(data(1:s-1));  datatype=data(s+1:end);
    end
else
    datalength=0; datatype=''; type='';
end
end

function [type data]=General_Information_Line(str)
s=find(str==':',1,'first');
if(isempty(s)), s=length(str); end
type=str(1:s-1);  data=str(s+1:end);
type=regexp(type, '\s+|/', 'split');
type_clean='';
for i=1:length(type)
    part=type{i};
    part(part=='.')=[]; part(part==' ')=[]; partu=uint8(part);
    if(~isempty(part))
        check=((partu>=97)&(partu<=122))|((partu>=65)&(partu<=90));
        if(check)
            part=lower(part); part(1)=upper(part(1));
            type_clean=[type_clean part];
        else
            break;
        end
    end
end
type=type_clean;
while(~isempty(data)&&data(1)==' '), data=data(2:end); end
while(~isempty(data)&&data(end)==' '), data=data(1:end-1); end
end

function V = par_read_volume(info)
% Function for reading the volume of a Philips Par / Rec  MR V4.* file 
%
% volume = par_read_volume(file-header)
%
% examples:
% 1: info = par_read_header()
%    V = par_read_volume(info);
%    imshow(squeeze(V(:,:,round(end/2),1)),[]);
%
% 2: V = par_read_volume('test.par');

if(~isstruct(info)), info=par_read_header(info); end

% Open file
fid=fopen(info.FilenameREC','rb','ieee-le');
% Skip header
fseek(fid,0,'bof');

datasize=prod(info.Dimensions)*info.BitDepth/8;

% Read the Data
switch(info.BitDepth)
    case 8
        info.DataType='char';
    case 16
        info.DataType='short';
    case 32
        info.DataType='float';
    case 64
        info.DataType='double';
end

switch(info.DataType)
    case 'char'
        V = int8(fread(fid,datasize,'char=>int8'));
    case 'uchar'
        V = uint8(fread(fid,datasize,'uchar=>uint8'));
    case 'short'
        V = int16(fread(fid,datasize,'short=>int16'));
    case 'ushort'
        V = uint16(fread(fid,datasize,'ushort=>uint16'));
    case 'int'
        V = int32(fread(fid,datasize,'int=>int32'));
    case 'uint'
        V = uint32(fread(fid,datasize,'uint=>uint32'));
    case 'float'
        V = single(fread(fid,datasize,'float=>single'));
    case 'double'
        V = double(fread(fid,datasize,'double=>double'));
end

fclose(fid);
V = reshape(V,info.Dimensions);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION Read Par-rec header:



function info = par_read_header(filename)
% Function for reading the header of a Philips Par / Rec  MR V4.* file 
%
% info  = par_read_header(filename);
%
% examples:
% 1,  info=par_read_header()
% 2,  info=par_read_header('volume.par');

if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.par', 'Read par-file');
    filename = [pathname filename];
end

fid=fopen(filename,'rb');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end
info.Filename=filename;


mode = -1; nHC=0; nIC=0; nSC=0;
while(true)
    str=fgetl(fid);
    if ~ischar(str), break, end
    if(isempty(str)), continue, end
    
    if(strfind(str,'= DATA DESCRIPTION FILE =')), mode=0; end
    if(strfind(str,'= GENERAL INFORMATION =')), mode=1; end
    if(strfind(str,'= PIXEL VALUES =')), mode=2; end
    if(strfind(str,'= IMAGE INFORMATION DEFINITION =')), 
        mode=3; fgetl(fid);  str=fgetl(fid); 
        % Skip a line
    end
    if(strfind(str,'= IMAGE INFORMATION =')), mode=4; end
    
    if(strfind(str,'= END OF DATA DESCRIPTION FILE =')), mode=5; end
    
    switch(mode)
        case -1;
        case 0
            nHC=nHC+1; HeaderComment{nHC}=str;
        case 1
            if(str(1)=='.')
                [type data]=General_Information_Line(str);
                switch(type)
                    case 'PatientName'
                        info.(type)=data;
                    case 'ProtocolName'
                        info.(type)=data;
                    case 'ExaminationName'
                        info.(type)=data;
                    case 'ExaminationDateTime'
                        info.(type)=data;
                    case 'SeriesType'
                        info.(type)=data;
                    case 'AcquisitionNr'
                        info.(type)=sscanf(data, '%d')';
                    case 'ReconstructionNr'
                        info.(type)=sscanf(data, '%d')';
                    case 'ScanDuration'
                        info.(type)=sscanf(data, '%lf')';
                    case 'MaxNumberOfCardiacPhases'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfEchoes'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfSlicesLocations'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfDynamics'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfMixes'
                        info.(type)=sscanf(data, '%d')';
                    case 'PatientPosition'
                        info.(type)=data;
                    case 'PreparationDirection'
                        info.(type)=data;
                    case 'Technique'
                        info.(type)=data;
                    case 'ScanResolution'
                        info.(type)=sscanf(data, '%d')';
                    case 'ScanMode'
                        info.(type)=data;
                    case 'RepetitionTime'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Fov'
                        info.(type)=sscanf(data, '%lf')';
                    case 'WaterFatShift'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Angulation'
                        info.(type)=sscanf(data, '%lf')';
                    case 'OffCentre'
                        info.(type)=sscanf(data, '%lf')';
                    case 'FlowCompensation'
                        info.(type)=sscanf(data, '%d')';
                    case 'Presaturation'
                        info.(type)=sscanf(data, '%d')';
                    case 'PhaseEncodingVelocity'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Mtc'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Spir'
                        info.(type)=sscanf(data, '%lf')';
                    case 'EpiFactor'
                        info.(type)=sscanf(data, '%lf')';
                    case 'DynamicScan'
                        info.(type)=sscanf(data, '%lf')';
                    case 'Diffusion'
                        info.(type)=sscanf(data, '%lf')';
                    case 'DiffusionEchoTime'
                        info.(type)=sscanf(data, '%lf')';
                    case 'MaxNumberOfDiffusionValues'
                        info.(type)=sscanf(data, '%d')';
                    case 'MaxNumberOfGradientOrients'
                        info.(type)=sscanf(data, '%d')';
                    case 'NumberOfLabelTypes'
                        info.(type)=sscanf(data, '%d')';
                    case 'HeaderComment'
                    otherwise
                        display(' did not recognize the information type: ');
                        display(type);
                        if (size(type) == 0)
                            display('empty, skipping line');
                        else
                            info.(type)=data;
                        end
                end
            end
        case 2
        case 3
            if(str(1)=='#');
                [type datatype datalength]=Image_Information_Line(str);
                if(~isempty(type))
                    nIC=nIC+1;
                    ImageInformationTags(nIC).Name=type;
                    ImageInformationTags(nIC).DataType=datatype;
                    ImageInformationTags(nIC).NumberOfValues=datalength;
                end
            end
        case 4
            if(str(1)~='#');
                nSC=nSC+1;
                vals=regexp(str, '\s+','split');
                vald=sscanf(str, '%lf')';
                current_loc=0;
                for i=1:length(ImageInformationTags)
                    IIT=ImageInformationTags(i);
                    if(strcmp(IIT.DataType,'string'))
                        SliceInformation(nSC).(IIT.Name)=vals{current_loc+1};
                    else
                        SliceInformation(nSC).(IIT.Name)=vald(current_loc+1:current_loc+IIT.NumberOfValues);
                    end
                    current_loc=current_loc+IIT.NumberOfValues;
                end
            end
            
        case 5
        otherwise
            %disp(str);
    end
end
fclose(fid);
info.HeaderComment=HeaderComment;
info.SliceInformation=SliceInformation;
info.ImageInformationTags=ImageInformationTags;

% Add Dimensions and Voxel Spacing. Warning, based on only 1 slice!
infof=info.SliceInformation(1);
if(isfield(infof,'ReconResolution'))
    if(isfield(info,'MaxNumberOfSlicesLocations'))
        zs(1)=info.MaxNumberOfSlicesLocations;
        zs(2)=length(SliceInformation)/zs(1);
        if((mod(zs(2),1)>0)||zs(2)==1)
            zs=length(SliceInformation);
        end
    else
        zs=length(SliceInformation);
    end
    
    info.Dimensions=[infof.ReconResolution zs];
else
    info.Dimensions=[info.ScanResolution length(SliceInformation)];
end
if(isfield(infof,'PixelSpacing'))
    if(isfield(infof,'SliceThickness')&&isfield(infof,'SliceGap'))
        zs=infof.SliceThickness+infof.SliceGap;
    else
        zs=0;
    end
    info.Scales=[infof.PixelSpacing zs];
else
    info.Scales=[0 0 0];
end


[folder,filen]=fileparts(info.Filename);
info.FilenameREC=fullfile(folder,[filen '.rec']);
   
    
% Add bith depth
if(infof.ImagePixelSize)
    info.BitDepth=infof.ImagePixelSize;
else
    if(exist(info.FilenameREC,'file'))
        file_info=dir(info.FilenameREC);
        bytes=file_info.bytes;
        info.BitDepth=(bytes/prod(info.Dimensions))*8;
    end
end
end


function [type datatype datalength]=Image_Information_Line(str)
s=find(str=='(',1,'last');
if(isempty(s)), s=length(str); end
type=str(1:s-1);  data=str(s+1:end);
type=regexp(type, '\s+|/|_', 'split');
type_clean='';
for i=1:length(type)
    part=type{i};
    part(part=='#')=[]; part(part==' ')=[]; partu=uint8(part);
    if(~isempty(part))
        check=((partu>=97)&(partu<=122))|((partu>=65)&(partu<=90));
        if(check)
            part=lower(part); part(1)=upper(part(1));
            type_clean=[type_clean part];
        else
            break;
        end
    end
end
type=type_clean;
while(~isempty(data)&&data(1)==' '), data=data(2:end); end
while(~isempty(data)&&data(end)==' '), data=data(1:end-1); end
if(~isempty(data))
    data=data(1:end-1);
    s=find(data=='*',1,'first');
    if(isempty(s)),
        datalength=1; datatype=data;
    else
        datalength=str2double(data(1:s-1));  datatype=data(s+1:end);
    end
else
    datalength=0; datatype=''; type='';
end
end


function [type data]=General_Information_Line(str)
s=find(str==':',1,'first');
if(isempty(s)), s=length(str); end
type=str(1:s-1);  data=str(s+1:end);
type=regexp(type, '\s+|/', 'split');
type_clean='';
for i=1:length(type)
    part=type{i};
    part(part=='.')=[]; part(part==' ')=[]; partu=uint8(part);
    if(~isempty(part))
        check=((partu>=97)&(partu<=122))|((partu>=65)&(partu<=90));
        if(check)
            part=lower(part); part(1)=upper(part(1));
            type_clean=[type_clean part];
        else
            break;
        end
    end
end
type=type_clean;
while(~isempty(data)&&data(1)==' '), data=data(2:end); end
while(~isempty(data)&&data(end)==' '), data=data(1:end-1); end

end