function [ Q, Qi, PWV_TTF, PWV_TTP, PWV_XCor ] = PulseWaveVelocity( im, hd, nCLpts, caseId )
%% Function to estimate pulse wave velocity given an image 'im' with header 'hd'

in = fopen( fullfile('Data','Adelaide','INPUT',caseId,'time_step.txt'), 'r' );
time_step = fscanf(in,'%f',[1 1]);
fclose(in);

%% Reduction of the size of image 'im' in order to speed up computation
[ im.b, ro, co, he ] = reduce_volume( im.b, 5 );
Mvenc = max( max(im.V(:)), abs(min(im.V(:)) ) );
im.V = im.V(:,:,ro,co,he);

% Unwrapping
im.V(:,1,:,:,:) = unwrap( im.V(:,1,:,:,:)/Mvenc*pi )/pi*Mvenc;
im.V(:,2,:,:,:) = unwrap( im.V(:,2,:,:,:)/Mvenc*pi )/pi*Mvenc;
im.V(:,3,:,:,:) = unwrap( im.V(:,3,:,:,:)/Mvenc*pi )/pi*Mvenc;

%% Update of relevant 'hd' fields after volume reduction of 'im'
hd.origin = hd.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
hd.dim = size( im.b );
hd.points = prod( hd.dim );
hd.Mv2w(:,4) = hd.origin';
hd.Mw2v = GetMw2v(hd.Mv2w);

%%
[ im, hd ] = transform_v2w( im, hd );
figure
subplot(1,2,1)
show_segment_surface( im.b, hd.Mv2w, [], 0.6250, 0.25 );
hold on

%% Calculation of points on centerline and corresponding directional vectors
BinaryOfCenterline = Skeleton3D( im.b );
nrrdIM = Build_nrrd( BinaryOfCenterline, hd );
[ pp, Length ] = SmoothSkeleton( nrrdIM, 1 );
evPoints0 = Length*0;
Length = Length*1;

nPoints = nCLpts;
evPoints = evPoints0 : ( Length - evPoints0 )/( nPoints - 1 ) : Length;
Point = fnval( pp, evPoints )';
ppder = fnder( pp );
Slope = fnval( ppder, evPoints )';
norm_Slope = sqrt( sum( Slope.*Slope, 2 ) );
Slope = diag(1./norm_Slope)*Slope;

%%
[ f, v ] = isosurface( im.b,  0.5 );
v = v - ones(size(v));
v = v(:,[2 1 3]);
vT = [ v ones(size(v,1),1)]*hd.Mv2w';

%%
bRecalibrateCenterline = 0;
if bRecalibrateCenterline
    sm = 20; % smoothing factor
    fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');
    [ pp, Length ] = recalibrateCenterline( f, vT, Point, Slope, sm );
    nPoints = nCLpts;
    evPoints = evPoints0 : ( Length - evPoints0 )/( nPoints - 1 ) : Length;
    Point = fnval( pp, evPoints )';
    ppder = fnder( pp );
    Slope = fnval( ppder, evPoints )';
    norm_Slope = sqrt( sum( Slope.*Slope, 2 ) );
    Slope = diag(1./norm_Slope)*Slope;
    hold on
    plot3( Point(:,1), Point(:,2), Point(:,3), '-r', 'LineWidth', 4 )
end

bDebug = 1;
r1 = mesh2crossSection( f, vT, Point(1,:),   Slope(1,:),   bDebug );
re = mesh2crossSection( f, vT, Point(end,:), Slope(end,:), bDebug );

if r1 < re
    sm = 20; % smoothing factor
    Point = Point(end:-1:1,:);
    pp = spaps( evPoints, Point', sm );
    ppder = fnder( pp );
    Slope = fnval( ppder, evPoints )';
    norm_Slope = sqrt( sum( Slope.*Slope, 2 ) );
    Slope = diag(1./norm_Slope)*Slope;
    fprintf('\n');
    fprintf('Centerline inverted');
    fprintf('\n');
end

for iP = 1 : nPoints
    im.P(iP).point = Point(iP,:);
    im.P(iP).slope = Slope(iP,:);
end

%% We permute im.V dimensions to improve performance in the next section
im.V = permute( im.V, [2 3 4 5 1] );
nTimeFrames = size( im.V, 5 );

%% Calculating flow at each point of the center line
% Flow is in mm3/s -> then converted to l/min
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');
Q = zeros( nTimeFrames, nPoints );
for iP = 1 : nPoints
    progress = floor(iP/nPoints*100);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('Progress:%3d%% of centerline', progress );
    [ im, ~, ~, bmask ] = intersect_image_with_plane( im, hd, iP );
    point = im.P(iP).point;
    slope = im.P(iP).slope;
    b = im.P(iP).b;
    sx = double(im.P(iP).dx(1) * im.P(iP).dx(2));
    for t = 1 : nTimeFrames;
        V = im.V(:,:,:,:,t);
        Vel = reshape( V ,size(V,1), size(V,2)*size(V,3)*size(V,4) );
        Vp = slope*Vel; % Component of velocity along normal
        Vp = reshape( Vp, size(V,2), size(V,3), size(V,4) );
        vPnrrd = Build_nrrd( Vp, hd );
        Vp2d = scinrrd_intersect_plane( vPnrrd, point, slope );
        Vp2d( isnan(Vp2d) ) = 0;
        Vp2d = Vp2d .* b;
        Vp = sum( Vp2d(:) );
        Q(t,iP) = Vp * sx; % Lambda
    end
end
fprintf('\n\n');

%% Displaying final data
Q = Q/1000; % from mm3/s to ml/s

subplot(1,2,2)
surf( evPoints, ([1:nTimeFrames]-1)*time_step*1000, Q );

axis tight
axis square
colormap( jet(10000) )
%camlight
xlabel('Centerline length (mm)')
ylabel('Time (ms)')
title( 'Flow rate (ml/s)' )
drawnow

savefig(['Results\Pulse Wave Velocity\shape&flow' caseId ]);

%% Smoothing spline surface fit to original data
t = 1 : nTimeFrames;
s = evPoints;
z = double(Q);
sp = spaps( {t,s}, z, 0);
%tstep = (t(end)-t(1))/1000;
tstep = 1;
%sstep = (s(end)-s(1))/200;
tishift = 12;
ti = t(1)-tishift : tstep : t(end);
%si = s(1)   : sstep : s(end);
si = s;
sp = fnxtr(sp,3);
Qi = fnval( sp, {ti, si} );

figure
surf( si, ti*time_step*1000, Qi,'LineStyle','none')
axis tight
axis square
colormap( jet(10000) )
%camlight
xlabel('Centerline length (mm)')
ylabel('Time (ms)')
title( 'Flow rate (ml/s)' )
drawnow

%%
%  Estimation of PWV with Time-To-Foot, Time-To-Peak, and XCor as in
%  "Estimation of Global Aortic Pulse Wave Velocity by Flow-Sensitive 4D MRI"
%  by Markl et al., Magnetic Resonance in Medicine 63:1575–1582 (2010)
for j = 1:size(Qi,2);
    spi = spaps( ti, Qi(:,j), 0 );
    %spid = fnder(spi);
    %TTF(j) = ( -fnval( spi, t(1) ) / fnval( spid, t(1) ) + t(1) )*time_step*1000; % Time-To-Foot
    z = fnzeros(spi);
    TTF(j) = z(1)*time_step*1000;
    [~, TTP(j)] = fnmin(fncmb(spi,-1)); % Time-To-Peak
    TTP(j) = TTP(j) * time_step * 1000;
end

ti0 = tishift/tstep + 1;
for j = 1:size(Qi,2);
    % XCor (cross correlation)
    [ acor, lag ] = xcorr( Qi(ti0:end,j), Qi(ti0:end,1) );
    [ ~, I ] = max( abs(acor) );
    XCor(j) = lag(I)*tstep * time_step * 1000;
end

%% Plotting PWV estimation results

 siz = round(length(si)/2);
 si(1:siz) = [];
 TTF(1:siz) = [];
 TTP(1:siz) = [];
 XCor(1:siz) = [];

figure

subplot(1,3,1)
f_TTF = fit( si', TTF', 'poly1' );
plot( f_TTF, si, TTF )
view(-90, 90) %# Swap the axes
set(gca, 'ydir', 'reverse'); %# Reverse the y-axis
%legend( 'Data', 'Fitted line', 'Location', 'northwest' )
legend off
axis square
xlabel('Centerline (mm)')
ylabel('Time (ms)')
PWV_TTF = 1 / f_TTF.p1;
title(sprintf('Estimation with TTF\n Pulse Wave Velocity = %0.2f m/s\n',PWV_TTF))

subplot(1,3,2)
f_TTP = fit( si', TTP', 'poly1' );
plot( f_TTP, si, TTP )
view(-90, 90) %# Swap the axes
set(gca, 'ydir', 'reverse'); %# Reverse the y-axis
%legend( 'Data', 'Fitted line', 'Location', 'northwest' )
legend off
axis square
xlabel('Centerline (mm)')
ylabel('Time (ms)')
PWV_TTP = 1 / f_TTP.p1;
title(sprintf('Estimation with TTP\n Pulse Wave Velocity = %0.2f m/s\n',PWV_TTP))

subplot(1,3,3)
f_XCor = fit( si', XCor', 'poly1' );
plot( f_XCor, si, XCor )
view(-90, 90) %# Swap the axes
set(gca, 'ydir', 'reverse'); %# Reverse the y-axis
%legend( 'Data', 'Fitted line', 'Location', 'northwest' )
legend off
axis square
xlabel('Centerline (mm)')
ylabel('Time (ms)')
PWV_XCor = 1 / f_XCor.p1;
title(sprintf('Estimation with XCor\n Pulse Wave Velocity = %0.2f m/s\n',PWV_XCor))
savefig(['Results\Pulse Wave Velocity\PWVEstimation' caseId ]);

end