function b = vel2bin( b, V, ths, bDebug )

if nargin < 3
    ths = 100;
end

if nargin < 4
    bDebug = 0;
end

b0 = b;

nFrames = size(V,1);
bin = zeros(nFrames,size(b,1),size(b,2),size(b,3));
cube = strel( ones(5,5,5) );
for iFrame = 1 : nFrames
    Vx = squeeze( V( iFrame, 1, :, :, : ) );
    Vy = squeeze( V( iFrame, 2, :, :, : ) );
    Vz = squeeze( V( iFrame, 3, :, :, : ) );
    Vm = sqrt( Vx.^2 + Vy.^2 + Vz.^2 );
    bd = imdilate( b, cube );
    bin(iFrame,:,:,:)  = b + ( (bd - b).*Vm > ths );
end

b = ceil( ( ceil(squeeze(mean(bin,1)) + b )/2 ) );

b = round( smooth3(b,'box',5) + 0.1 );

if bDebug
    [ f, v ] = isosurface( b, 0.5 );
    patch('Faces',f,'Vertices',v,'facecolor', [1 1 0], 'facealpha', 0.33,'edgecolor','none');
    axis equal
    camlight
    hold on
    
    [ f0, v0 ] = isosurface( b0,  0.5 );
    patch('Faces',f0,'Vertices',v0,'facecolor', [1 0 0], 'facealpha', 1,'edgecolor','none');
end