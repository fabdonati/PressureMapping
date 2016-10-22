function PWV = RunAdvectiveCenterline( dataDir, iCase, dCP, werp, pwv ) %#ok<INUSL>

% 'dCP' = distance between centerline points in mm

close all

if nargin < 5;   pwv = [];  end
if nargin < 4;  werp = [];  end
if nargin < 3;   dCP = 10;  end
if nargin < 2; iCase = '2'; end

opts.werp = werp;
opts.pwv  = pwv;

mainDir = cd;
j = 0;
while ~strcmpi(mainDir(end-j+1:end),'PressureMapping')
    mainDir(end-j:end) = [];
    j = 0;
    while ~( mainDir(end-j) == '\' || mainDir(end-j) == '/' );
        j = j + 1;
    end
end

rmpath( genpath( mainDir ) );
addpathi(mainDir,'IO');
addpathi(mainDir,'Visualization');
addpathi(fullfile(mainDir,'ThirdParty'),'Skeleton3D');
addpathi(fullfile(mainDir,'ThirdParty'),'Gerardus');
%addpathi(fullfile(mainDir,'ThirdParty'),'NIfTI_20140122');
addpathi(fullfile(mainDir,'ThirdParty','matlab_bgl-4.0.1'),'matlab_bgl');
addpathi(mainDir,'Tools');

functionIn = ...
    'AdvectiveCenterline( mainDir, dataDir, iCase, ''HT'', ''sys_seg.vtk'', ''sys_seg_ref.vtk'', ''Velocity_rs.mat'', ''HLHS'', 1, dCP, opts );';

if strcmpi(opts.werp,'bWerp') && strcmpi(opts.pwv,'bPwv')
    functionOut = '[ as, ao, ls, dist, werp, uB, PWV ]';
end

if strcmpi(opts.werp,'bWerp') && ~strcmpi(opts.pwv,'bPwv')
    functionOut = '[ as, ao, ls, dist, werp, uB ]';
end

if ~strcmpi(opts.werp,'bWerp') && strcmpi(opts.pwv,'bPwv')
    functionOut = '[ as, ao, ls, dist, ~, ~, PWV ]';
end

if ~strcmpi(opts.werp,'bWerp') && ~strcmpi(opts.pwv,'bPwv')
    functionOut = '[ as, ao, ls, dist ]';
end

eval( [ functionOut ' = ' functionIn ]);

rmpath( genpath( mainDir ) );

display( 'Finished!' );
fprintf( '\n' );

function addpathi(mainDir,folder)
d = dir( mainDir );
i = 1;
while i <= length(d) && ~( strcmpi( d(i).name, folder ) && d(i).isdir )
    i = i + 1;
end
addpath( fullfile(mainDir,d(i).name) );