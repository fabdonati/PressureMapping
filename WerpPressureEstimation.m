function [werp] = WerpPressureEstimation(bVTK,Velocity,opts)
% Input files:
%	bVTK:		Name of the binary mask in vtk format
%	Velocity: 	Velocity (Ntimeframes, Ncomponents, Nx, Ny, Nz) in mm/s
%   opts:   structure with
%       dt:		Time spacing dt (in s)
%       rho:	Fluid density in kg/mm3
%       mu:		Fluid dynamic viscosity in kg/mm/s
%
% Version control
% - v0.1: code clean up by Pablo Lamata, from code by Fabrizio Donati, and
% uploaded to BitBucket on 8th Jan 2015

defaultDT = 0.040;
defaultRHO= 1060  *1e-9;
defaultMU = 0.004 *1e-3;
if nargin<3
    opts  = [];
end  
if ~isfield(opts,'dt'),   opts.dt = defaultDT; end
if ~isfield(opts,'rho'),  opts.rho = defaultRHO; end
if ~isfield(opts,'mu'),   opts.mu = defaultMU; end

% Further options:
opts.interp2d   = 'nn';             % Interpolation scheme for the plane intersection with binary image ('nn' / 'linear')
opts.stencil    = 'filtered';       % Finite differences scheme stencil ('standard' / 'filtered')
opts.timescheme = 'cdt2';           % Time derivative scheme ('cdt' for standard cdiff, O((dt)^2) / 'cdt2' for cdiff at mid points, O((dt/2)^2)
opts.MeshName   = 'AORTA';          % Output name for PPE mesh
opts.compliance = 0;
opts.AdveFilter = 0;

[xx,yy,zz] = ndgrid(-1:1); opts.se = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1.0;


% WERP SOLVER
to = 2;  Velocity(isnan(Velocity)) = 0;
Nt = size(Velocity,1);
Nc = size(Velocity,2);

im.V = Velocity;

switch opts.timescheme 
  case 'cdt', Vtype = 'V';  to = 2;
  case 'cdt2', Vtype = 'Vh'; 
    for t = 1:Nt-1, for i = 1:Nc, im.Vh(t,i,:,:,:) = (im.V(t+1,i,:,:,:) + im.V(t,i,:,:,:))/2; end, end
end

[im.b, hd] = io_ReadMedicalImage( bVTK );
[im, hd]   = transform_v2w( im, hd );

if exist('P','var'), 
    im.P = P;
else
    im.P = select_planes( bVTK, im, hd, opts ); 
end, Np = length(im.P);

plot_image_with_selected_planes( im, hd );

disp('Preliminary computations of ROI and WALL masks...')
im = define_roi_mask( im, hd );
im = define_wall_mask( im, opts );

for t = to:Nt-1, disp(['Frame ' num2str(t)]);   
  im = intersect_image_with_planes( im, Vtype, hd, opts, t );
  
  % Compute lambda
  werp.lambdai(t)  = compute_lambda_2d( im.P(1).V, -im.P(1).N, im.P(1).b, im.P(1).dx, opts.stencil );
  werp.lambdao(t)  = compute_lambda_2d( im.P(end).V, im.P(end).N, im.P(end).b, im.P(end).dx, opts.stencil );
  
  % Compute advective energy rate
  for i = 1:Np/2, [a(i), ~] = compute_advective_2d( im.P(i).V, -im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho );  end
  for i = Np:-1:Np/2+1, [a(i), ~] = compute_advective_2d( im.P(i).V, im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho ); end
  avar = opts.AdveFilter * max(abs(a(1) - a(2)),abs(a(end) - a(end-1)));
  werp.a(t)  = ((a(1) + a(end))^2   / ((a(1) + a(end))^2   + avar^2) ) * (a(1) + a(end));
  
  % Compute viscous dissipation rate
  [werp.v(t), ~] = compute_viscous_3d( squeeze(eval(['im.' Vtype '(t,:,:,:,:)'])), im.bROI, hd.spacing, opts.stencil, opts.mu);
  
  % Compute kinetic energy rate
  [kp, ~] = compute_kinetic_energy_3d( squeeze(im.V(t+1,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
  [km, ~] = compute_kinetic_energy_3d( squeeze(im.V(t-1 * strcmp(opts.timescheme,'cdt'),:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
  
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

plot_results( im, hd, werp); toc