function [ werp, uB ] = WERP( ... 
    CenterlineObj, points, normals, ls, as, ai, opts, ... 
    im, hd, t0, tE, Vtype, out, defaultDT )

%%
% Blood density is 1060 Kg/M3, here converted in Kg/mm3
defaultRHO =  1060 * 1e-9;

%% Selection of input and output points of centerline
display('Calculating WERP')
display('----------------');
display('Select ROI');

id1mm = input('First point on centerline (in mm): ');
id2mm = input(' Last point on centerline (in mm): ');

P1 = ( fnval(CenterlineObj.F,id1mm) )';
P2 = ( fnval(CenterlineObj.F,id2mm) )';
id1 = knnsearch(points,P1);
id2 = knnsearch(points,P2);

%%
nPoints = size(points,1);

%% Compute lambda
werp.lambdai = -ls( id1 , : );
werp.lambdao =  ls( id2, : );

%% Compute advective energy rate
ao = -as.*ls/(1e3/133.33) - ( ai'*ones(1,nPoints) )';
a1 = -ao(id1,  :);
a2 = -ao(id1+1,:);
aF =  ao(id2-1,:);
aE =  ao(id2,  :);
avar = opts.AdveFilter * max( abs(a1-a2), abs(aE-aF) );
werp.a = ( ((a1 + aE).^2) ./ ((a1 + aE).^2 + avar.^2) ) .* (a1 + aE);

disp('Computing ROI mask...')
im = define_roi_mask( im, hd, id1, id2 );
%disp('Computing the WALL mask...')
%im = define_wall_mask( im, opts, id1, id2 );

%% Compute id for voxels on centerline
r = round( [ points ones(nPoints,1) ]*hd.Mw2v' );

%%
for t = t0 : tE
    disp(['Frame ' num2str(t)]);
    %% Compute viscous dissipation rate
    werp.v(t) = compute_viscous_3d( squeeze(eval(['im.' Vtype '(t,:,:,:,:)'])), ...
        im.bROI, hd.spacing, opts.stencil, opts.mu);
    
    %% Compute kinetic energy rate
    kp = compute_kinetic_energy_3d( squeeze(im.V(t+1,:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
    km = compute_kinetic_energy_3d( squeeze(im.V(t-1 * strcmp(opts.timescheme,'cdt'),:,:,:,:)), im.bROI, hd.spacing, opts.stencil, opts.rho);
    switch opts.timescheme,
        case 'cdt'
            werp.k(t) = (kp - km) / (2 * opts.dt);
        case 'cdt2'
            werp.k(t) = (kp - km) / (opts.dt);
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
end

%% Compute pressure drops Unsteady Bernoulli
uB.point =  points( id1 : id2, : );
uB.slope = normals( id1 : id2, : );

for t = t0 : tE
    disp(['Frame ' num2str(t)]);
    index = 1;
    blank = im.b * 0;
    for i = id1 : id2
        if blank( r(i,1), r(i,2), r(i,3) ) == 0
            uB.V(t,index,:) = mean3d( squeeze(im.V(t,:,:,:,:)),r(i,1),r(i,2),r(i,3) );
            uB.PV(t,index) = dot(squeeze(uB.V(t,index,:)),uB.slope(index,:));
            blank( r(i,1), r(i,2), r(i,3) ) = 1;
            index = index + 1;
        end
    end
end

uB.dpP_iner = InertialBernoulli(uB.point,uB.PV,defaultRHO,defaultDT,'centered2',0);
uB.dpP_adve = AdvectiveBernoulli(uB.PV,defaultRHO)';
uB.dpP      = uB.dpP_iner + uB.dpP_adve;

%% Plot
figure
hold on
plot( uB.dpP/0.0075,      '-ko', 'linewidth', 2 )
plot( uB.dpP_iner/0.0075, '-bo', 'linewidth', 2 )
plot( uB.dpP_adve/0.0075, '-go', 'linewidth', 2 )
title('Unsteady Bernoulli')
legend('UB Total','UB Inertial','UB Advective', 'Location', 'NorthEast');

figure(1)
show_segment_surface( im.bROI, hd.origin, hd.spacing, .8 ); hold on

surf( im.P( id1 ).gx, im.P( id1 ).gy, im.P( id1 ).gz, squeeze(max( im.P( id1 ).Vp,[],1)), ... 
    'FaceColor', 'interp', 'FaceLighting', 'none', 'LineStyle', 'none' );
colormap(jet)
surf( im.P( id2 ).gx, im.P( id2 ).gy, im.P( id2 ).gz, squeeze(max( im.P( id2 ).Vp,[],1)), ... 
    'FaceColor', 'interp', 'FaceLighting', 'none', 'LineStyle', 'none' );
colormap(jet)
drawnow
plot_results( werp, id1, id2, opts.interp2d );

%{
%% Save
if out == 1
    filename = [name '_werp_' num2str(id1) '_' num2str(id2) '.mat'];
    filenameUB = [name '_ub_' num2str(id1) '_' num2str(id2) '.mat'];
    save([folder filename],'werp');
    save([folder filenameUB],'uB');
    savefig([folder 'ROI.fig'])
end
%}