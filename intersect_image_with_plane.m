function [ im, v0, v1, bmask ] = intersect_image_with_plane( im, hd, iP )

P = im.P(iP);

Pnrrd = Build_nrrd(im.b,hd);
options.interp = 'nearest';
[ btmp, gx, gy, gz, midx, v0, v1 ] = scinrrd_intersect_plane( Pnrrd, P.point, P.slope, options );
btmp(isnan(btmp)) = 0; % Sobstitute NaNs with zeros
btmp = ceil(btmp); % Transforming all nonzero entries to 1
% Select the component of the cross section that contains the center point
btmp = bwselect( btmp, midx(2), midx(1), 8 );

im.P(iP).b = btmp;

im.P(iP).gx = gx;
im.P(iP).gy = gy;
im.P(iP).gz = gz;

im.P(iP).dx = evaluate_pixeldim( gx, gy, gz );

points =  [ gx( btmp ), gy(btmp), gz( btmp ) ];
pts = [ points ones(size(points,1),1) ] * hd.Mw2v';
pts = pts + ones(size(pts));
pts = floor( pts );
%pts = pts*[0 1 0; 1 0 0; 0 0 1];
    
B = im.b;
bmask = zeros(size(B));
for i = 1 : size(pts,1);
    bmask( pts(i,1), pts(i,2), pts(i,3) )=1;
end

%{
hold on
plot3( gx( btmp ), ...
    gy( btmp ), ...
    gz( btmp ), 'r*' );
plot3( gx(midx(1), midx(2)), ...
       gy(midx(1), midx(2)), ...
       gz(midx(1), midx(2)), 'bo' )
axis equal
axis tight
drawnow
%}