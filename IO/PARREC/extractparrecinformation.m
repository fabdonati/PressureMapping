%> @file extractparrecinfoormation.m
%> 
%> @brief extract direction matrix and origin according to the image angulation given
%> in the par-rec file from philips scanner
%>
%> function [outputdirectionmatrix, outputorigin] = extractdirections (sliceorientation, [ sliceangulations ], [ offcenter ], [ spacing ], [ dimensions ])
%> extract direction matrix and offset according to the 
%> image angulation and other info given in the par-rec file from philips scanner
%> 
%> @param sliceangulations are the slice angulations given in the par-rec file (or
%> the DICOM headers) ap fh rl in degrees
%> @param sliceorientation is an integer that gives the slice orientation
%> allowed values are 1, 2, or 3 respectively correspond to transversal, sagittal and coronal
%> orientations
%> @param spacing is the given spacing in x y and z of the image, to calculate the
%> correct offset
%> @param dimension is the given dimension in x y and z of the image, to calculate the
%> correct offset
%> 
%> @param offcenter is the given offcenter in x y and z, to calculate the
%> correct offset
%> 
%> @param outputdirections is the directions to use (for instance in an mha
%> header
%> 
%> Usage : [M,O] = extractdirections (2, [23.2 34.2 3.1], [12.4 27.0 46.2], [2.1 2.1 4.2], [128 128 30])
%>
%> @author Nicolas Toussaint nicolas.toussaint@kcl.ac.uk


function [outputdirections, outputorigin] = extractparrecinformation (sliceorientation, sliceangulations, offcenter, spacing, dimensions)

%% define the change of system
AFRtoLPS = [
     0     0     1
     1     0     0
     0     1     0
];


%% obscure reference matrix
referencematrix = [
    -1  0  0
     0  0  1 
     0 -1  0
];

%% reference matrices from slice orientations
tra = [
     0  1  0
    -1  0  0
     0  0 -1
    ];

sag = [
    -1  0  0
     0  0  1
     0 -1  0
    ];

cor = [
     0  1  0
     0  0  1
    -1  0  0
    ];

%% choose the orientation matrix from sliceorientation
switch sliceorientation
    case 1,
        Torientation = tra;
    case 2,
        Torientation = sag;
    case 3,
        Torientation = cor;
    otherwise,
        Torientation = eye(3);
end

%% put the angulations in SI units...
sliceangulations = sliceangulations .* pi ./ 180.0;

%% calculate the re-orientation matrices from the angulations
Tap = [
    1  0          0
    0  cos(sliceangulations(1))  -sin(sliceangulations(1))
    0  sin(sliceangulations(1))   cos(sliceangulations(1))
    ];

Tfh = [
    cos(sliceangulations(2))  0  sin(sliceangulations(2))
    0         1  0
    -sin(sliceangulations(2))  0  cos(sliceangulations(2))
    ];

Trl = [
    cos(sliceangulations(3))  -sin(sliceangulations(3))  0
    sin(sliceangulations(3))   cos(sliceangulations(3))  0
    0          0         1
    ];

%% the final reorientation matrix, stunning !
outputdirections = AFRtoLPS * Trl * Tap * Tfh * referencematrix' * Torientation';


%% calculate offset to apply to go from middle of the volume to corner pixel
midoffset = spacing' .* (dimensions' - ones(3,1)) ./ 2.0;
midoffset = outputdirections * midoffset;

%% the final origin of the image.
outputorigin = AFRtoLPS * offcenter' - midoffset;

end