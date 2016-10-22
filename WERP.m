function [werp,uB] = WERP(directory,root,VTK,VTKCl,VTKRef,Velmat,type,region,out,opts)
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
%% e.g. >> [werp,ub] = WERP('/home/adv10/DATA/PRESSURE_ESTIMATION/','HT', 'sys_seg.vtk','sys_seg.vtk','sys_seg_ref.vtk','Velocity_rs.mat','HLHSII','Ascending',1);


bLamndaAcrossLength = 0;

%% ADD ALL FOLDERS to path, first PressureMapping, then cardiacmeshing
addpath(genpath('/home/adv10/DATA/PRESSURE_ESTIMATION/Pressuremapping'));
addpath(genpath('/home/adv10/DATA/PRESSURE_ESTIMATION/cardiacmeshing'));

%% Load IDs of planes delimiting the ROI, calculated from Advective_centreline.m 
%% valve plane identified by the minimum in the advecitve comp along the centreline at peak systole
load('ascending_cl.mat');
%hlhs(8,:)=[6 25];
for i = 1:1
    patient = i;
    if i<10
        name = [root '00' num2str(i)];
    else 
        name = [root '0' num2str(i)];
    end    
    
    %% Output
    folder = [directory type '/' name '/PressureEstimation/' region '/' ];
    mkdir(folder);
    
    %% Input
    bVTK = [directory 'INPUT/' name '/' VTK];
    bVTKCl = [directory 'INPUT/' name '/' VTKCl];
    bVTKRef = [directory 'INPUT/' name '/' VTKRef];
    Vmat=load([directory 'INPUT/' name '/' Velmat]);
    Velocity = Vmat.Velocity;
    size(Velocity)

    in=fopen([directory 'INPUT/' name '/time_step.txt'],'r');
    defaultDT = fscanf(in,'%f',[1 1]);
    fclose(in);
    defaultRHO= 1060  *1e-9;
    defaultMU = 0.004 *1e-3;
    defaultbManualPlanes = 1;
    
    if nargin<10
        opts  = [];
    end  
    if ~isfield(opts,'dt'),   opts.dt = defaultDT; end
    if ~isfield(opts,'rho'),  opts.rho = defaultRHO; end
    if ~isfield(opts,'mu'),   opts.mu = defaultMU; end
    if ~isfield(opts,'bManualPlanes'),   opts.bManualPlanes = defaultbManualPlanes; end

    %% Further options:
    opts.interp2d   = 'nn';             % Interpolation scheme for the plane intersection with binary image ('nn' / 'linear')
    opts.stencil    = 'filtered';       % Finite differences scheme stencil ('standard' / 'filtered')
    opts.timescheme = 'cdt2';           % Time derivative scheme ('cdt' for standard cdiff, O((dt)^2) / 'cdt2' for cdiff at mid points, O((dt/2)^2)
    opts.MeshName   = 'AORTA';          % Output name for PPE mesh
    opts.compliance = 0;
    opts.AdveFilter = 0;

    [xx,yy,zz] = ndgrid(-1:1); opts.se = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1.0;


    %% WERP SOLVER SET-up
    to = 2;  
    Velocity(isnan(Velocity)) = 0;
    Nt = size(Velocity,1);
    Nc = size(Velocity,2);

    im.V = Velocity;
    im.Vh = Velocity;

    switch opts.timescheme 
      case 'cdt', Vtype = 'V'; to = 2;
      case 'cdt2', Vtype = 'Vh'; 
        for t = 1:Nt-1, for i = 1:Nc, im.Vh(t,i,:,:,:) = (im.V(t+1,i,:,:,:) + im.V(t,i,:,:,:))/2; end, end
    end

    %% Read image for PD + image for reference landmarks + image for centreline calc
    [im.b, hd] = io_ReadMedicalImage( bVTK );
    [imRef.b, hdRef] = io_ReadMedicalImage( bVTKRef );
    [imCl.b, hdCl] = io_ReadMedicalImage( bVTKCl );
    [imCl, hdCl]   = transform_v2w( imCl, hdCl );
    [imRef, hdRef]   = transform_v2w( imRef, hdRef );
    [im, hd]   = transform_v2w( im, hd );

    %% compute centreline
    [CL, Plane,BinaryOfCenterline] = get_planes_and_centerline(bVTKCl);
    plot_image_with_points( imRef.b,hdRef,Plane ); title([name]);savefig([folder  'CL.fig'])
    Pcl = Plane.point';
    
    %% Compute id for voxels on centreline
    for i = 1:size(Plane.point,1)
      r(i,:) = round(GetMw2v(hd.Mv2w)*[Plane.point(i,:),1]');
    end
    
    %% Extract ROI from centreline planes
    if strcmp(type,'HLHSII')
        ppt = hlhs(patient,:);
    else
        ppt = agematched(patient,:);
    end
    
    if exist('P','var'), 
        imCl.P = P;
    else
        [im.P, id1, id2] = select_planes2(CL, Plane, ppt, opts ); %% ADV
        Pcl_sel = zeros(3,1);
        Pextra = id1-1; %% check orientation of cl-- HACK!
        [ver] = show_segment_surface(imCl.b,hdCl.Mv2w);
        vertex = ver.vertices;
        for i = 1:norm(id1-id2)
            Pcl_sel = Pcl(:,id1+i-1);
            Vec = [Pcl(:,id1+i)-Pcl_sel]'; 
            Pcl_sel = Pcl_sel';
            [point,Radius(i)] = cross_sec_new(vertex,Pcl_sel,Vec);
            if mod(i,2) == 0
                plot3(point(:,1),point(:,2),point(:,3),'x'); 
                plot3(Pcl_sel(:,1),Pcl_sel(:,2),Pcl_sel(:,3),'rs');
            end
        end
        %% Compute cross-sectional radius at each plane in ROI
        Rave=mean(Radius)
        figure, plot(Radius);title('Radius');
        ylabel('R (mm)'); xlabel('length (mm)');
        save([folder name '_radius' num2str(id1) '_' num2str(id2) '.mat'],'Radius'); 
    end, Np = length(im.P);

    disp('Preliminary computations of ROI and WALL masks...')
    im = define_roi_mask( im, hd );
    im = define_wall_mask( im, opts );
    
    if (bLamndaAcrossLength)
        Hlam = figure; hold on;
    end
    
    for t = to:Nt-1, disp(['Frame ' num2str(t)]);   
      im = intersect_image_with_planes( im, Vtype, hd, opts, t );

      %% Compute lambda
      
      werp.lambdai(t)  = compute_lambda_2d( im.P(1).V, -im.P(1).N, im.P(1).b, im.P(1).dx, opts.stencil );
      werp.lambdao(t)  = compute_lambda_2d( im.P(end).V, im.P(end).N, im.P(end).b, im.P(end).dx, opts.stencil );

      if (bLamndaAcrossLength)
          figure(Hlam);
          im2 = im;
          opts.bAllPoints = 1;
          [im2.P, id1, id2] = select_planes(bVTKCl, im, hd, imCl, hdCl, opts);

          im2 = intersect_image_with_planes( im2, Vtype, hd, opts, t );

          nPoints = numel(im2.P);
          for iP = 1:nPoints
              l(iP) = compute_lambda_2d( im2.P(iP).V, -im2.P(iP).N, im2.P(iP).b, im2.P(iP).dx, opts.stencil );
          end
          plot(l); 
      end
      
      %% Compute advective energy rate
      for i = 1:Np/2, [a(i), ~] = compute_advective_2d( im.P(i).V, -im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho );  end
      for i = Np:-1:Np/2+1, [a(i), ~] = compute_advective_2d( im.P(i).V, im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho ); end
      avar = opts.AdveFilter * max(abs(a(1) - a(2)),abs(a(end) - a(end-1)));
      werp.a(t)  = ((a(1) + a(end))^2   / ((a(1) + a(end))^2   + avar^2) ) * (a(1) + a(end));

      %% Compute viscous dissipation rate
      [werp.v(t), ~] = compute_viscous_3d( squeeze(eval(['im.' Vtype '(t,:,:,:,:)'])), im.bROI, hd.spacing, opts.stencil, opts.mu);

      %% Compute kinetic energy rate
      [kp, ~] = compute_kinetic_energy_3d( squeeze(im.V(t+1,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
      [km, ~] = compute_kinetic_energy_3d( squeeze(im.V(t-1 * strcmp(opts.timescheme,'cdt'),:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
      switch opts.timescheme, 
          case 'cdt',  werp.k(t) = (kp - km) / (2 * opts.dt); 
          case 'cdt2', werp.k(t) = (kp - km) / (opts.dt); 
      end

      %% Compute pressure drops WERP
      werp.pdo(t)  = - 1 / werp.lambdao(t) * (werp.k(t) + werp.a(t) + werp.v(t));
      werp.pdi(t)  =   1 / werp.lambdai(t) * (werp.k(t) + werp.a(t) + werp.v(t));
      werp.kpdo(t) = - 1 / werp.lambdao(t) * werp.k(t);
      werp.apdo(t) = - 1 / werp.lambdao(t) * werp.a(t);
      werp.vpdo(t) = - 1 / werp.lambdao(t) * werp.v(t);
      werp.kpdi(t) =   1 / werp.lambdai(t) * werp.k(t);
      werp.apdi(t) =   1 / werp.lambdai(t) * werp.a(t);
      werp.vpdi(t) =   1 / werp.lambdai(t) * werp.v(t);

      %% Compute pressure drops Unsteady Bernoulli
      index = 1;
      [id2_ub] = knnsearch(Plane.point,Pcl(:,id2)');
      [id1_ub] = knnsearch(Plane.point,Pcl(:,id1)');
      blank = zeros(size(BinaryOfCenterline));
      for i = id1_ub:1:id2_ub    
        if blank(r(i,1),r(i,2),r(i,3)) == 0
          uB.V(t,index,:) = MeanV(squeeze(im.V(t,:,:,:,:)),r(i,1),r(i,2),r(i,3));
          uB.point(index,:) = Plane.point(i,:);
          uB.slope(index,:) = Plane.slope(i,:);
          uB.PV(t,index) = dot(squeeze(uB.V(t,index,:)),uB.slope(index,:));
          blank(r(i,1),r(i,2),r(i,3)) = 1;
          index = index + 1;
        end
      end   

    end
    uB.dpP = InertialBernoulli(uB.point,uB.PV,defaultRHO,defaultDT,'centered2',0) + ...
      AdvectiveBernoulli(uB.PV,defaultRHO)';
    uB.dpP_iner =InertialBernoulli(uB.point,uB.PV,defaultRHO,defaultDT,'centered2',0);
    uB.dpP_adve = AdvectiveBernoulli(uB.PV,defaultRHO)';

    
    %% plot
    figure, 
    plot(uB.dpP/0.0075,'k','linewidth',2),hold on
    plot(uB.dpP_iner/0.0075,'b','linewidth',2),hold on
    plot(uB.dpP_adve/0.0075,'g','linewidth',2),hold on
    legend('UB TOTAL','UB INERTIAL','UB ADVECTIVE');
    
    figure,
    show_segment_surface(imRef.b,hdRef.origin,hdRef.spacing); hold on,
    show_segment_surface(im.bROI,hd.origin,hd.spacing,.8);
    surf(im.P(1).gx,im.P(1).gy,im.P(1).gz,im.P(1).Vp); hold on,
    surf(im.P(end).gx,im.P(end).gy,im.P(end).gz,im.P(end).Vp);
    
    plot_results( im, hd, werp, id1, id2, opts.interp2d);

    %% save
    if out == 1
        filename = [name '_werp_' num2str(id1) '_' num2str(id2) '.mat'];
        filenameUB = [name '_ub_' num2str(id1) '_' num2str(id2) '.mat'];
        save([folder filename],'werp');
        save([folder filenameUB],'uB');
        savefig([folder 'ROI.fig'])
    end

    clear Velocity bVTKCl bVTKRef bVTK
    
end
