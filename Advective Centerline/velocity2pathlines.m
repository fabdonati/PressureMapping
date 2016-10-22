function P = velocity2pathlines( im, opts, P, xvec, yvec, zvec )

% Function to calculate the trajectory of points P given the velocity
% vector field over time
dt = opts.dt;
rho = opts.rho;
mu = opts.mu;
x1 = unique(squeeze( im.W(1,:,:,:) ));
x2 = unique(squeeze( im.W(2,:,:,:) ));
x3 = unique(squeeze( im.W(3,:,:,:) ));
nFrames = size(im.V,1);
t = ( 0 : .25 : nFrames-1 )*dt;
for it = 1 : length(t)-1;
    t0 = floor(t(it)/dt+1);
    t1 = ceil(t(it)/dt+1);
    m = t(it)/dt+1 - t0;
    Vx = im.V(t0,1,:,:,:)*(1-m) + im.V(t1,1,:,:,:)*m;
    Vy = im.V(t0,2,:,:,:)*(1-m) + im.V(t1,2,:,:,:)*m;
    Vz = im.V(t0,3,:,:,:)*(1-m) + im.V(t1,3,:,:,:)*m;
    Vx = squeeze(Vx); Vx = permute( Vx, [2 1 3]);
    Vy = squeeze(Vy); Vy = permute( Vy, [2 1 3]);
    Vz = squeeze(Vz); Vz = permute( Vz, [2 1 3]);
    vxP = interp3( xvec, yvec, zvec, Vx, P(:,1,it), P(:,2,it), P(:,3,it) );
    vyP = interp3( xvec, yvec, zvec, Vy, P(:,1,it), P(:,2,it), P(:,3,it) );
    vzP = interp3( xvec, yvec, zvec, Vz, P(:,1,it), P(:,2,it), P(:,3,it) );
    vP = [vxP vyP vzP];
    P(:,:,it+1) = P(:,:,it) + vP*( t(it+1) - t(it) );
end