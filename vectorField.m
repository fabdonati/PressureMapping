function vectorField(caseId, tFrame)

close all

load( [ 'TestData\Adelaide\' caseId '\Velocity_rs.mat' ] );
im.V = Velocity;
nT = size( im.V, 1 );
for t = 1:nT-1
    im.Vh(t,:,:,:,:) = ( im.V(t+1,:,:,:,:) + im.V(t,:,:,:,:) ) / 2;
end
[ im.b, hd ] = io_ReadMedicalImage( [ 'TestData\Adelaide\' caseId '\sys_seg.vtk' ] );
% Some images have landmark in im.b, the following statement transform them
% to 1
im.b  = ceil( im.b/max(im.b(:)) );

%% Reduction of the size of image 'im' in order to speed up computation
[ im.b, ro, co, he ] = reduce_volume( im.b );
im.V  = im.V(:,:,ro,co,he);
im.Vh = im.Vh(:,:,ro,co,he);

%% Update of relevant 'hd' fields after volume reduction of 'im'
hd.origin = hd.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
hd.dim = size( im.b );
hd.points = prod( hd.dim );
hd.Mv2w(:,4) = hd.origin';
hd.Mw2v = GetMw2v(hd.Mv2w);

%%
for t = 1 : size(im.V,1);
    for iC = 1 : 3
        im.V(t,iC,:,:,:) = squeeze(im.V(t,iC,:,:,:)).*im.b;
    end
end

%%
[ im, hd ] = transform_v2w( im, hd );
u = squeeze( im.V(tFrame,1,:,:,:) );
v = squeeze( im.V(tFrame,2,:,:,:) );
w = squeeze( im.V(tFrame,3,:,:,:) );
x = squeeze( im.W(1,:,:,:) );
y = squeeze( im.W(2,:,:,:) );
z = squeeze( im.W(3,:,:,:) );
makevtk_struc_grid_FEX