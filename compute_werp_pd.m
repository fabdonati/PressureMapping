function[ im, hd, werp ]  = compute_werp_pd( filename, bVTK, opts )

load(filename);
opts = get_problem_variables(opts);
if exist('Velocity','var') 
  if isfield(opts,'velocity_convention')
    if strcmp(opts.velocity_convention,'P')
      Velocity(isnan(Velocity)) = 0;
      Nc = size(Velocity,4); Nt = size(Velocity,5); 
      im.V = permute(Velocity,[5 4 1 2 3]);
      if(1), im.V = im.V(:,:,end:-1:1,:,:); end
    elseif strcmp(opts.velocity_convention,'F')
      im.V = Velocity * 1000;
      Nc = size(Velocity,2); Nt = size(Velocity,1);
    end
  else disp('Error! Velocity information is not available'); return
  end
  if exist('bVTK','var')
    [im.b, hd] = io_ReadMedicalImage( bVTK );
    [im, hd]   = transform_v2w( im, hd );
    if exist('P','var'), im.P = P; else im.P = select_planes( bVTK, im, hd, opts ); end
    Np = length(im.P);
    plot_image_with_selected_planes( im, hd );
    im = define_roi_mask( im, hd );
    im = define_wall_mask( im, opts );
    for t = 1:Nt-1, for i = 1:Nc, im.Vh(t,i,:,:,:) = (im.V(t+1,i,:,:,:) + im.V(t,i,:,:,:))/2; end, end
    for t = to:Nt-1
      disp(['Frame ' num2str(t)]);
      im = intersect_image_with_planes( im, 'Vh', hd, opts, t );
      
      %     if t == to
      %       if opts.compliance == 0
      %         disp('Wall variables and compliance model are being generated...');
      %         im = define_compliance_mask( im );
      %       end
      %       generate_cheart_from_image(im.W, im.bROI, opts);
      %     end
      %     save_cheart_velocity(squeeze(im.V(t,:,:,:,:)), opts, t);
      
      % Compute lambda
      werp.lambdai(t)  = compute_lambda_2d( im.P(1).V, -im.P(1).N, im.P(1).b, im.P(1).dx, opts.stencil );
      werp.lambdao(t)  = compute_lambda_2d( im.P(end).V, im.P(end).N, im.P(end).b, im.P(end).dx, opts.stencil );
      
      % Compute advective energy rate
      for i = 1:Np/2
        [a(i), ~] = compute_advective_2d( im.P(i).V, -im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho );
      end
      for i = Np:-1:Np/2+1
        [a(i), ~] = compute_advective_2d( im.P(i).V, im.P(i).N, im.P(i).b, im.P(i).dx, opts.stencil, opts.rho );
      end
      atmp = mean(a);
      aiotmp = a(1) + a(end);
      werp.apre(t) = aiotmp;
      avar = opts.AdveFilter * max(abs(a(1) - a(2)),abs(a(end) - a(end-1)));
      werp.ain(t) = a(1);
      werp.aout(t) = a(6);
      werp.a(t)  = (aiotmp^2   / (aiotmp^2   + avar^2) ) * aiotmp;
      
      % Compute viscous dissipation rate
      [werp.v(t), ~] = compute_viscous_3d( squeeze(im.Vh(t,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.mu);
      
      % Compute kinetic energy rate
      [kp, ~] = compute_kinetic_energy_3d( squeeze(im.V(t+1,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
      [km, ~] = compute_kinetic_energy_3d( squeeze(im.V(t,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
      werp.k(t) = (kp - km) / (opts.dt);
      
      % Compute pressure drops
      werp.pd(t)  = - 1 / werp.lambdao(t) * (werp.k(t) + werp.a(t) + werp.v(t));
      werp.kpd(t) = - 1 / werp.lambdao(t) * werp.k(t);
      werp.apd(t) = - 1 / werp.lambdao(t) * werp.a(t);
      werp.vpd(t) = - 1 / werp.lambdao(t) * werp.v(t);
      
    end
  else
    disp('Error! Binary mask is not defined'); return
  end
else
  disp('Error! Velocity is not defined'); return
end
