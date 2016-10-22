function werp = WerpPressureEstimation_NEW( directory, name, bVTK, bVTKCl, bVTKRef, Velocity, type, region, opts )
%function [werp, uB] = WerpPressureEstimation(directory,name,bVTK,bVTKCl,bVTKRef,Velocity,type,region,opts)
% Input files:
%	bVTK*:		Name of the binary mask in vtk format
%	Velocity: 	Velocity (Ntimeframes, Ncomponents, Nx, Ny, Nz) in mm/s
%   opts:   structure with
%       dt:		Time spacing dt (in s)
%       rho:	Fluid density in kg/mm3
%       mu:		Fluid dynamic viscosity in kg/mm/s
%
% Version control
% - v0.1: code clean up by Pablo Lamata, from code by Fabrizio Donati, and
% uploaded to BitBucket on 8th Jan 2015

%% Output
mkdir( [ directory '/' type '/' name '/PressureEstimation' ] );
folder = [ directory '/' type '/' name '/PressureEstimation/' region '_' ];

%% Input
bVTK    = [ directory '/INPUT/' name '/' bVTK     ];
bVTKCl  = [ directory '/INPUT/' name '/' bVTKCl   ];
bVTKRef = [ directory '/INPUT/' name '/' bVTKRef  ];
load(     [ directory '/INPUT/' name '/' Velocity ] );

in = fopen([directory '/INPUT/' name '/' 'time_step.txt'],'r');
defaultDT = fscanf(in,'%f',[1 1]);
fclose(in);

defaultRHO= 1060  *1e-9;
defaultMU = 0.004 *1e-3;
defaultbManualPlanes = 1;

if nargin < 9, opts = []; end
if ~isfield(opts,'dt' ), opts.dt  = defaultDT;  end
if ~isfield(opts,'rho'), opts.rho = defaultRHO; end
if ~isfield(opts,'mu' ), opts.mu  = defaultMU;  end
if ~isfield(opts,'bManualPlanes'), opts.bManualPlanes = defaultbManualPlanes; end

% Further options:
opts.interp2d   = 'nn';             % Interpolation scheme for the plane intersection with binary image ('nn' / 'linear')
opts.stencil    = 'filtered';       % Finite differences scheme stencil ('standard' / 'filtered')
opts.timescheme = 'cdt2';           % Time derivative scheme ('cdt' for standard cdiff, O((dt)^2) / 'cdt2' for cdiff at mid points, O((dt/2)^2)
opts.MeshName   = 'AORTA';          % Output name for PPE mesh
opts.compliance = 0;
opts.AdveFilter = 0;

[xx,yy,zz] = ndgrid(-1:1); opts.se = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1.0;

% WERP solver
to = 2;
Velocity( isnan(Velocity) ) = 0;
Nt = size(Velocity,1);
Nc = size(Velocity,2);

im.V = Velocity;

switch opts.timescheme
    case 'cdt', Vtype = 'V'; to = 2;
    case 'cdt2', Vtype = 'Vh';
        for t = 1:Nt-1
            for i = 1:Nc
                im.Vh(t,i,:,:,:) = (im.V(t+1,i,:,:,:) + im.V(t,i,:,:,:))/2;
            end
        end
end

[    im.b, hd    ] = io_ReadMedicalImage( bVTK    );
[ imRef.b, hdRef ] = io_ReadMedicalImage( bVTKRef ); % Add reference image with landmarks for visual plane selection and centreline 
[  imCl.b, hdCl  ] = io_ReadMedicalImage( bVTKCl  );

[ ~, ro1, co1, he1 ] = reduce_volume(    im.b );
[ ~, ro2, co2, he2 ] = reduce_volume( imRef.b );
[ ~, ro3, co3, he3 ] = reduce_volume(  imCl.b );

rofirst = min( [ ro1(1) ro2(1) ro3(1) ] );
cofirst = min( [ co1(1) co2(1) co3(1) ] );
hefirst = min( [ he1(1) he2(1) he3(1) ] );
roend = max( [ ro1(end) ro2(end) ro3(end) ] );
coend = max( [ co1(end) co2(end) co3(end) ] );
heend = max( [ he1(end) he2(end) he3(end) ] );
ro = rofirst : roend;
co = cofirst : coend;
he = hefirst : heend;

im.V  = im.V(:,:,ro,co,he);
im.Vh = im.Vh(:,:,ro,co,he);
   im.b =    im.b(ro,co,he);
imRef.b = imRef.b(ro,co,he);
 imCl.b =  imCl.b(ro,co,he);

   hd.origin =    hd.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
hdRef.origin = hdRef.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
 hdCl.origin =  hdCl.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
 
    hd.dim = size(    im.b );
 hdRef.dim = size( imRef.b );
  hdCl.dim = size(  imCl.b );
  
     hd.points = prod(    hd.dim );
  hdRef.points = prod( hdRef.dim );
   hdCl.points = prod(  hdCl.dim );
   
   hd.Mv2w(:,4) = hd.origin';
   hd.Mw2v = GetMw2v(hd.Mv2w);
   
   hdRef.Mv2w(:,4) = hdRef.origin';
   hdRef.Mw2v = GetMw2v(hdRef.Mv2w);
   
   hdCl.Mv2w(:,4) = hdCl.origin';
   hdCl.Mw2v = GetMw2v(hdCl.Mv2w);

[imCl,    hdCl ] = transform_v2w( imCl, hdCl );
[im,      hd   ] = transform_v2w( im,   hd   );

if exist('P','var'),
    imCl.P = P;
else
    if(opts.bManualPlanes)
        figure
        hold on
        Centerline = GetCenterline( bVTKCl, 0, ro, co, he );
        Centerline.nPointsSpline = floor(Centerline.nPointsSpline);
        Pcl = GetAllPoints(Centerline);
        Centerline.PlotSpline();
        for iP = 1:Centerline.nPointsSpline
            Centerline.PlotSplinePoint(iP);
        end
        show_segment_surface(  imCl.b,  hdCl.Mv2w );
        show_segment_surface( imRef.b, hdRef.Mv2w ); %% ADV
    end
    
    [im.P, id1, id2, Pall] = select_planes( bVTKCl, im, hd, imRef, hdRef, opts, ...
                                             ro, co, he ); %% ADV
    ver = show_segment_surface( imCl.b, hdCl.Mv2w );
    vertex = ver.vertices;
    %figure, hold on; show_segment_surface( imCl.b, hdCl.Mv2w );
    for i = 1 : norm(id1-id2)
        Pcl_sel = Pcl( :, id1+i-1 );
        Vec = Pcl( :, id1+i ) - Pcl_sel;
        [ point, Radius(i) ] = cross_sec_new( vertex, Pcl_sel', Vec' );
        if mod(i,2) == 0
            hold on
            plot3(point(:,1),point(:,2),point(:,3),'x');
            plot3(Pcl_sel(1),Pcl_sel(2),Pcl_sel(3),'rs');
        end
    end
    
    figure
    plot(Radius)
    title('Radius')
    xlabel('Central line point')
    ylabel('R (mm)')
    save( [ folder 'radius' num2str(id1) '_' num2str(id2) '.mat' ],'Radius' );
end

plot_image_with_selected_planes( im, hd );

%%
im2 = im;
for i = 1 : size( Pall.point, 1 )
    P(i).point = Pall.point(i,:);
    P(i).slope = Pall.slope(i,:);
end
im2.P = P;
for iFrame = 2 : Nt-1
    iFrame
im2 = intersect_image_with_planes( im2, Vtype, hd, opts, iFrame );
im  = intersect_image_with_planes( im,  Vtype, hd, opts, iFrame );
%% Flux is in mm3/s -> then converted to l/min
    ain  = compute_advective_2d( im.P(  1).V,  -im.P(1).N,   im.P(1).b,   im.P(1).dx, opts.stencil, opts.rho );
    aout = compute_advective_2d( im.P(end).V, im.P(end).N, im.P(end).b, im.P(end).dx, opts.stencil, opts.rho );
for iP = 1 : length(Pall.point)
    Q(iP,iFrame) = compute_lambda_2d( im2.P(iP).V, im2.P(iP).N, im2.P(iP).b, im2.P(iP).dx, opts.stencil );
    DPA(iP) = (ain+aout)/Q(iP);
end

%figure, plot( Q*60/10^6       ); title( 'Flow rate (l/min)'          );
%figure, plot( DPA*1000/133.33 ); title( 'Advective PD in/out (mmHg)' );
end

disp('Preliminary computations of ROI and WALL masks...')
im = define_roi_mask( im, hd );
im = define_wall_mask( im, opts );

Np = length(im.P);
for t = to : Nt-1
    disp( ['Frame ' num2str(t)] );
    im = intersect_image_with_planes( im, Vtype, hd, opts, t );
    
    % Compute lambda
    werp.lambdai(t)  = compute_lambda_2d( im.P(  1).V, -im.P(  1).N, im.P(  1).b, im.P(  1).dx, opts.stencil );
    werp.lambdao(t)  = compute_lambda_2d( im.P(end).V,  im.P(end).N, im.P(end).b, im.P(end).dx, opts.stencil );
    
    % Compute advective energy rate
    for i = 1 : Np/2;
        a(i) = compute_advective_2d( im.P(i).V, -im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho ); 
    end
    for i = Np : -1 : Np/2+1
        a(i) = compute_advective_2d( im.P(i).V,  im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho );
    end
    avar = opts.AdveFilter * max( abs(a(1) - a(2)), abs(a(end) - a(end-1)) );
    werp.a(t)  = ((a(1) + a(end))^2   / ((a(1) + a(end))^2   + avar^2) ) * (a(1) + a(end));
    
    % Compute viscous dissipation rate
    werp.v(t) = compute_viscous_3d( squeeze(eval(['im.' Vtype '(t,:,:,:,:)'])), im.bROI, hd.spacing, opts.stencil, opts.mu);
    
    % Compute kinetic energy rate
    kp = compute_kinetic_energy_3d( squeeze(im.V(t+1,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
    km = compute_kinetic_energy_3d( squeeze(im.V(t-1 * strcmp(opts.timescheme,'cdt'),:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
    
    switch opts.timescheme,
        case 'cdt',  werp.k(t) = (kp - km) / (2 * opts.dt);
        case 'cdt2', werp.k(t) = (kp - km) / (opts.dt);
    end
    
    % Compute pressure drops
    werp.pdo(t)  = - 1 / werp.lambdao(t) * (werp.k(t) + werp.a(t) + werp.v(t));
    werp.pdi(t)  =   1 / werp.lambdai(t) * (werp.k(t) + werp.a(t) + werp.v(t));
    werp.kpdo(t) = - 1 / werp.lambdao(t) * werp.k(t);
    werp.apdo(t) = - 1 / werp.lambdao(t) * werp.a(t);
    werp.vpdo(t) = - 1 / werp.lambdao(t) * werp.v(t);
    werp.kpdi(t) =   1 / werp.lambdai(t) * werp.k(t);
    werp.apdi(t) =   1 / werp.lambdai(t) * werp.a(t);
    werp.vpdi(t) =   1 / werp.lambdai(t) * werp.v(t);
end

%% ADV
plot_results( im, hd, werp, id1, id2, opts.interp2d); %toc

filename = ['werp_' num2str(id1) '_' num2str(id2) '.mat'];
save([folder filename],'werp');
savefig([folder 'ROI.fig'])
