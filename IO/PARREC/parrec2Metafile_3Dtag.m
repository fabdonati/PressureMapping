function parrec2Metafile(rootfilename,gatingType)
%input: rootfilename.PAR and rootfilename.REC (Philips format)
%output: rootfilename.mhd and rootfilename.raw (metafile format in the
%correct dicom space)

%OPTIONAL input: gatingType; values: 'retrospective' (default) or 'prospective'

%check existence of input files and get rid of the extension if it was provided:
if (strcmp(rootfilename(end-3:end),'.PAR') || ...
		strcmp(rootfilename(end-3:end),'.par') || ...
		strcmp(rootfilename(end-3:end),'.rec') || ...
		strcmp(rootfilename(end-3:end),'.REC'))
	rootfilename = rootfilename(1:end-4);
end

if exist([rootfilename '.PAR'],'file')
	inputPAR = [rootfilename '.PAR'];
elseif  exist([rootfilename '.par'],'file')
	inputPAR = [rootfilename '.par'];
else
	error('Input PAR file does not exist!');
end

%load header and in particular the image orientation and position info from PAR file:
hdr = par_read_header(inputPAR);
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
%3 types of [TransformMatrix1, Offset1] for the 3D tag data:
[TransformMatrix1, Offset1] = extractparrecinformationfromfile(hdr);

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
	for j = 1 : hdr.Dimensions(3)
		for k = 1 : hdr.Dimensions(4)
			triggerTime(j,k) = hdr.SliceInformation(j+(k-1)*hdr.Dimensions(3)).TriggerTime;
		end
	end
	triggerTime = triggerTime / 1000.0; %milliseconds --> seconds
	
	%number of data sets in a single PAR-REC file; it is 2 typically if the
	%PAR-REC contains magnitude and phase image.
	noOfDataComp = 1;
	if (mod(hdr.Dimensions(4),2) == 0)
		if (triggerTime(1,1) == triggerTime(1,(hdr.Dimensions(4)/2) + 1))
			warning('Input 3D+t image seems to consist of two type of data, e.g. magnitude and image');
			noOfDataComp = 2;
		end
	end
	
	timeStep = zeros(hdr.Dimensions(3),1);  %time step for each slice (has sense only in retrospective triggering, it should be constant in the prospective one)
	for j = 1 : hdr.Dimensions(3)
		timeStep(j) = (triggerTime(j,end) - triggerTime(j,1)) / ((hdr.Dimensions(4)/noOfDataComp) - 1);
	end
	meanTimeStep = mean(timeStep);
	
	%For the prospective trigger, the mid-point scheme is applied:
	Offset(4) = triggerTime(1,1);  %timeOffset, should be 0 for retrospective trigger
	if strcmp(gatingType,'prospective')
		warning('For the prospective trigger, the mid-point time scheme is applied');
		Offset(4) = Offset(4) + (meanTimeStep / 2.0);
	end
	
	if (noOfDataComp == 2)
		Dimensions = hdr.Dimensions;
		Dimensions(4) = hdr.Dimensions(4) / 2;
		%first component of the data set:
        
        warning('Hack for 3D tag. Only magnitude image is saved.');
        %all three orientation are in one image with a gemoetrical offset:
        dimZ = Dimensions(3)/3;
        
		WriteMhaFile([rootfilename '_A.mhd'], [Dimensions(1) Dimensions(2) dimZ Dimensions(4)], [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
		WriteBinary([rootfilename '_A.raw'],data(:,:,1:dimZ,1:(hdr.Dimensions(4)/2)),dataType);
       
        TransformMatrix2 = TransformMatrix;
        TransformMatrix2(1:3,1:3) = TransformMatrix3D * [0 -1 0; 1 0 0; 0 0 1];  %rotation by 90 degrees in xy: 
        
        WriteMhaFile([rootfilename '_B.mhd'], [Dimensions(1) Dimensions(2) dimZ Dimensions(4)], [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix2,[1 prod(size(TransformMatrix))]));
		WriteBinary([rootfilename '_B.raw'],data(:,:,dimZ+1:2*dimZ,1:(hdr.Dimensions(4)/2)),dataType);
        WriteMhaFile([rootfilename '_C.mhd'], [Dimensions(1) Dimensions(2) dimZ Dimensions(4)], [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
		WriteBinary([rootfilename '_C.raw'],data(:,:,2*dimZ+1:3*dimZ,1:(hdr.Dimensions(4)/2)),dataType);
        
        summedData = data(:,:,1:dimZ,1:(hdr.Dimensions(4)/2)) + data(:,:,dimZ+1:2*dimZ,1:(hdr.Dimensions(4)/2)) + data(:,:,2*dimZ+1:3*dimZ,1:(hdr.Dimensions(4)/2));
        
        WriteMhaFile([rootfilename '_SummedImage.mhd'], [Dimensions(1) Dimensions(2) dimZ Dimensions(4)], [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
		WriteBinary([rootfilename '_SummedImage.raw'],summedData(:,:,:,1:(hdr.Dimensions(4)/2)),dataType);
        
		%second component of the data set:
		%WriteMhaFile([rootfilename '_2.mhd'], Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
		%WriteBinary([rootfilename '_2.raw'],data(:,:,:,(hdr.Dimensions(4)/2)+1:end),dataType);
		%warning(['Assuming that the 2nd component image is the phase image, it will be scaled to cm / sec and written as ' rootfilename '_2scaled.raw']);
		%V = pixel value in REC file, FP = floating point value, DV = displayed value on console
		%RS = rescale slope,           RI = rescale intercept,    SS = scale slope
		%DV = PV * RS + RI
		%data2Scaled = data(:,:,:,(hdr.Dimensions(4)/2)+1:end);
		%data2Scaled = (double(data(:,:,:,(hdr.Dimensions(4)/2)+1:end))*hdr.SliceInformation(end).RescaleSlope + hdr.SliceInformation(end).RescaleIntercept) * sum(hdr.PhaseEncodingVelocity) / (1000.0*pi);		
		%WriteMhaFile([rootfilename '_2scaled.mhd'], Dimensions, [hdr.Scales meanTimeStep], dataType, Offset, reshape(TransformMatrix,[1 prod(size(TransformMatrix))]));
		%WriteBinary([rootfilename '_2scaled.raw'],data2Scaled,dataType);
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