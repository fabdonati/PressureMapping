function [ as, ao, ls, dist, werp, uB, PWV ] = ...
    AdvectiveCenterline( mainDir, dataDir, iCase, root, VTK, VTKref, velMat, type, out, dCP, opts )
% Input files:
%	bVTK*:		Name of the binary mask in vtk format
%	Velocity: 	Velocity (Ntimeframes, Ncomponents, Nx, Ny, Nz) in mm/s
%   opts:   structure with
%       dt:		Time spacing dt (in s)
%       rho:	Fluid density in kg/mm3
%       mu:		Fluid dynamic viscosity in kg/(mm*s)
%
% Version control
% - v0.1: code clean up by Pablo Lamata, from code by Fabrizio Donati, and
% uploaded to BitBucket on 8th Jan 2015

%%
iCaseName = num2str(str2num(iCase)+1000);
caseName = [ root iCaseName(2:end) ];

loadDir = fullfile( dataDir, 'INPUT', type, caseName );
saveDir = fullfile( mainDir, 'Results', 'Pressure Estimation', type, caseName );
mkdir(saveDir);

bReset = 0;
bRunMain = 1;
if exist(fullfile(saveDir,'dCP.mat'),'file') && ~bReset
    load(fullfile(saveDir,'dCP.mat'));
    if dCP == dCP_old; %#ok<NODEF>
        bRunMain = 0;
        load( fullfile( saveDir, 'adv_lambda.mat' ) );
    end
end

if bRunMain
    delete(fullfile(saveDir,'PWV.mat'));
    %% Loading data
    bVTK    =       fullfile( loadDir, VTK    );
    bVTKref =       fullfile( loadDir, VTKref );
    velMatS = load( fullfile( loadDir, velMat ) );
    
    %% Unwrapping (if necessary) of velocity data
    if ~isfield( velMatS, 'VelocityU' )
        VelocityU = unwrapVelocity( velMatS.Velocity );
        save( fullfile( loadDir, velMat ), 'VelocityU', '-append' );
        velMatS.VelocityU = VelocityU;
    end
    Velocity = velMatS.VelocityU;
   
    %%
    for iFrame = 1 : size(Velocity,1)
        Velocity(iFrame,1,:,:,:) = smooth3(squeeze(Velocity(iFrame,1,:,:,:)),'gaussian',3);
        Velocity(iFrame,2,:,:,:) = smooth3(squeeze(Velocity(iFrame,2,:,:,:)),'gaussian',3);
        Velocity(iFrame,3,:,:,:) = smooth3(squeeze(Velocity(iFrame,3,:,:,:)),'gaussian',3);
    end
    
    %% Parameters
    
    %Time spacing dt (in s)
    in = fopen( fullfile( loadDir, 'time_step.txt' ), 'r' );
    defaultDT = fscanf(in,'%f',[1 1]);
    time_step = defaultDT;
    fclose(in);
    
    % Blood density is 1060 Kg/M3, here converted in Kg/mm3
    defaultRHO =  1060 * 1e-9;
    
    % Blood viscosity is 0.004 Kg/(M*s), here converted in Kg/(mm*s)
    defaultMU  = 0.004 * 1e-3;
    
    %% Option initialisation
    if nargin<10
        opts  = [];
    end
    if ~isfield(opts,'dt'),   opts.dt  = defaultDT;  end
    if ~isfield(opts,'rho'),  opts.rho = defaultRHO; end
    if ~isfield(opts,'mu'),   opts.mu  = defaultMU;  end
    
    %% Further options
    opts.interp2d   = 'linear';   % Interpolation scheme for the plane intersection with binary image ('nn' / 'linear')
    opts.stencil    = 'filtered'; % Finite differences scheme stencil ('standard' / 'filtered')
    opts.timescheme = 'cdt';     % Time derivative scheme ('cdt' for standard cdiff, O((dt)^2) / 'cdt2' for cdiff at mid points, O((dt/2)^2)
    opts.MeshName   = 'AORTA';    % Output name for PPE mesh
    opts.compliance = 0;
    opts.AdveFilter = 0;
    
    [ xx, yy, zz ] = ndgrid( -1 : 1 );
    opts.se = sqrt( xx.^2 + yy.^2 + zz.^2 ) <= 1.0;
    
    %% Import binary mask to calculate the centerline
    % Warning! Look at big and little endian !!!!!
    %bVTK='C:\Users\Alessandro\OneDrive\Matlab\PressureMapping\Data\Nidhin\HT002.vtk';
    [ im.b, hd ]       = io_ReadMedicalImage( bVTK );
    [ imRef.b, hdRef ] = io_ReadMedicalImage( bVTKref );
    
    %%
    %[ im.b, hd ] = allignAnatBinary( im.b, hd, imRef.b, hdRef, [0 -4 -4] );
    
    %% Updating binary mask according to velocity values over time (testing needed)
    %im.b = vel2bin( im.b, im.V, 200, 0 );
    % This might not work because the aorta deform over time, so it will
    % make more sense to only consider a velocity adapted binary mask for each
    % time frame.
    
    %% Smoothing binary mask for better centerline calculation
    im.b = round( smooth3( im.b, 'box', 3) + 0.1);
    
    %% Volume reduction of binary image and velocity matrix to speed up computation
    [ im.b, ro, co, he ] = reduce_volume( im.b, 15 );
    imRef.b =  imRef.b(ro,co,he);
    im.V = Velocity(:,:,ro,co,he);
    
    %%
    nFrames =  size( im.V, 1 );
    switch opts.timescheme
        case 'cdt'
            Vtype = 'V';
        case 'cdt2'
            Vtype = 'Vh';
            im.Vh = ( im.V(1:(nFrames-1),:,:,:,:) + im.V(2:nFrames,:,:,:,:) ) / 2;
            nFrames = size( im.Vh, 1 );
    end
    
    %% Update of relevant header 'hd' fields after volume reduction
    hd.origin    = hd.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
    hd.dim       = size( im.b );
    hd.points    = prod( hd.dim );
    hd.Mv2w(:,4) = hd.origin';
    hd.Mw2v      = GetMw2v(hd.Mv2w);
    
    hdRef.origin    = hdRef.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
    hdRef.dim       = size( imRef.b );
    hdRef.points    = prod( hdRef.dim );
    hdRef.Mv2w(:,4) = hdRef.origin';
    hdRef.Mw2v      = GetMw2v(hdRef.Mv2w);
    
    %% From voxel to world coordinates
    [ im,    hd    ] = transform_v2w( im, hd );
    [ imRef, hdRef ] = transform_v2w( imRef, hdRef );
    
    %% Save velocity field as VTK for each time frame
    fprintf('\n')
    %str = input('Would you like to save the velocity field in VTK?(Y/N) ','s' );
    str = 'N';
    fprintf('\n')
    if strcmpi(str,'y')
        for bWithMask = 0:1;
            B = im.b;
            sWithMask = 'withMask';
            if ~bWithMask
                B = ones(size(B));
                sWithMask = '';
            end
            for iFrame = 1 : nFrames;
                V = squeeze(im.V(iFrame,:,:,:,:));
                V(1,:,:,:) = squeeze( V(1,:,:,:) ).*B;
                V(2,:,:,:) = squeeze( V(2,:,:,:) ).*B;
                V(3,:,:,:) = squeeze( V(3,:,:,:) ).*B;
                filename = [ 'velocity' sWithMask num2str(iFrame) '.vtk'];
                filename = fullfile(loadDir,filename);
                vectorField2vtk( filename, V, hd )
            end
        end
        display('VTK files created!');
        fprintf('\n')
    end
    
    %% B-form of velocity vector field over time
    %{
    t = (0:nFrames-1)*opts.dt;
    x1 = unique(squeeze( im.W(1,:,:,:) ));
    x2 = unique(squeeze( im.W(2,:,:,:) ));
    x3 = unique(squeeze( im.W(3,:,:,:) ));
    V = permute( im.V, [2 3 4 5 1]);
    spV = spaps( {x1,x2,x3,t}, V, 0, [] , 1 );
    %}
    
    %% Plot of segmentation
    h1 = figure(1);
    show_segment_surface( im.b,    hd.origin,    hd.spacing, 0.6250, 0.25 );
    show_segment_surface( imRef.b, hdRef.origin, hdRef.spacing, 0.3750, 0.25 );
    title( caseName )
    
    %% Compute centerline
    [ CenterlineObj, points, normals, nPoints ] = get_centerline( im.b, hd, dCP );
    
    %% Inverting centerline if necessary (from ascending to descending aorta)
    [ f, v ] = isosurface( im.b,  0.5 );
    v = v - ones(size(v));
    v = v(:,[2 1 3]);
    vT = [ v ones(size(v,1),1)]*hd.Mv2w';
    
    bDebug = 0;
    r1 = mesh2crossSection( f, vT, points(1,:),   normals(1,:),   bDebug );
    re = mesh2crossSection( f, vT, points(end,:), normals(end,:), bDebug );
    bInvertPoints = 0;
    if r1 < re
        bInvertPoints = 1;
        fprintf('\n');
        fprintf('Centerline inverted');
        fprintf('\n');
    end
    
    %% Adding extra points at the begiining of centerline
    bExtraPoints = 1;
    [ CenterlineObj, points, normals, nPoints ] = ...
        get_centerline( im.b, hd, dCP, bInvertPoints, bExtraPoints );
    
    %% Plot of centerline
    bDebug = 1;
    r1 = mesh2crossSection( f, vT, points(1,:),   normals(1,:),   bDebug );
    re = mesh2crossSection( f, vT, points(end,:), normals(end,:), bDebug );
    hold on
    plot3( points(:,1), points(:,2), points(:,3), 'r.', 'MarkerSize', 10 )
    %{
    L = ( 0 : size(points,1)-1 )*dCP;
    pointIdx = cell(length(L),1);
    for i = 1 : length(L), pointIdx(i) = { L(i) };end
    text( points(:,1), points(:,2), points(:,3), pointIdx, ...
        'HorizontalAlignment', 'left', 'rotation', 90, 'FontWeight', 'bold', ...
        'FontSize', 8)
    %}
    fnplt( CenterlineObj.F, 0.5 )
    drawnow
    
    %% Main
    t0 = 1;
    tE = nFrames;
    p0 = 1;
    pE = nPoints;
    
    [ as, ao, ai, ls, dist, im.P ] = compute_quantities_1D( ...
        im, Vtype, hd, opts, CenterlineObj, p0, pE, t0, tE );
    as( 1 : t0-1, : ) = [];
    ao( 1 : t0-1, : ) = [];
    ls( 1 : t0-1, : ) = [];
    
    as = as';
    ao = ao';
    ls = ls';
    
    dist = dist';
    
    %% Plot & save
    
    h2 = figure;
    
    subplot(131)
    plot( dist, as )
    hold on
    axis tight
    xlabel( 'Length of aorta (mm)' )
    ylabel( 'mmHg' )
    title( [ caseName ' Advective component' ] )
    
    subplot(132)
    plot( dist, -ao )
    hold on
    axis tight
    xlabel('Length of aorta (mm)')
    ylabel('mmHg')
    title( [ caseName ' Advective Energy Rate' ] )
    
    subplot(133)
    plot( dist, ls/1000*60/1000 )
    hold on
    axis tight
    xlabel( 'Length of aorta (mm)' )
    ylabel( 'l/min' )
    title( [caseName ' Flow'] )
    
    if out == 1
        savefig( h1, fullfile( saveDir, 'centerline.fig' ) );
        savefig( h2, fullfile( saveDir, 'adv_lambda.fig' ) );
        dCP_old = dCP;
        save(    fullfile( saveDir, 'dCP.mat' ), 'dCP_old' );
        save(    fullfile( saveDir, 'adv_lambda.mat' ), ...
            'CenterlineObj', 'points', 'normals', 'ls', 'as', 'ao', 'ai', ...
            'dist', 'opts', 'im', 'hd', 't0', 'tE', 'Vtype', 'out', 'defaultDT' );
    end
end

%% WERP
werp = []; uB = [];
if strcmpi(opts.werp,'bWerp')
    close all
    openfig( fullfile( saveDir, 'centerline.fig' ) );
    openfig( fullfile( saveDir, 'adv_lambda.fig' ) );
    [ werp, uB ] = WERP( ...
        CenterlineObj, points, normals, ls, as, ai, opts, ...
        im, hd, t0, tE, Vtype, out, defaultDT );
end

%% Pulse Wave Velocity
if strcmpi(opts.pwv,'bPwv')
    time_step = defaultDT;
    openfig( fullfile( saveDir, 'centerline.fig' ) );
    load('C:\Users\Alessandro\OneDrive\Matlab\PressureMapping\Data\INPUT\HLHS\aortaPoints.mat');
    for i = 1 : size(aortaPoints,2)-1
        newIni = aortaPoints(str2num(iCase),i);
        newEnd = aortaPoints(str2num(iCase),i+1);
        j = (i-1)*3;
        [ PWV(1+j), PWV(2+j), PWV(3+j) ] = PulseWaveVelocity( ls, dist, t0, tE, time_step, saveDir, newIni, newEnd );
    end
end