function [ r, cs ] = calculate_cross_section( vol, P, N )

% P is a point on the smoothed central line
% vec is the direction vector of the central line at point P

R = P + N;

% Compute rotation angles
angleh = 0; Rot = eye(3);
u(1) =    R(2) - P(2);
u(2) = -( R(1) - P(1) );
u(3) = 0;
uN = u / norm(u);
v = [ 1 0 0 ];
if u(1)~= 0 || u(2) ~= 0
    ruv = vrrotvec(u,v);
    angleh = ruv(4)/(2*pi)*360;
    Rot = vrrotvec2mat(ruv);
end

if uN(2) < 0
    angleh = acos(uN(1)) / (2*pi)*360;
end
if uN(2) > 0
    angleh = 360 - acos(uN(1)) / (2*pi)*360;
end
if uN(2) == 0
    if uN(1) == 1
        angleh = 0;
    end
    if uN(1) == -1
        angleh = 180;
    end
end

% Rotations around point P
%[vol_r , P_r ] = volrotateP_c( vol  , P  , anglex );
[vol_r , P_r ] = volrotateP_h( vol  , round(P)  , angleh );

R_r = (Rot*u')' + [ P_r(2) -P_r(1) 0 ];
R_r = [ -R_r(2) R_r(1) R(3) ];

angley = 0;
u(1) = R_r(2) - P_r(2);
u(2) = R_r(3) - P_r(3);
u(3) = 0;
uN = u / norm(u);
v = [ 0 1 0 ];
if u(1)~= 0 || u(2) ~= 0
    ruv = vrrotvec(u,v);
    angley = ruv(4)/(2*pi)*360;
end

if uN(1) < 0
    angley = 360 - ( 90 - asin(uN(2)) / (2*pi)*360 );
end
if uN(1) > 0
    angley = 90 - asin(uN(2)) / (2*pi)*360;
end
if uN(1) == 0
    if uN(2) == 1
        angley = 0;
    end
    if uN(2) == -1
        angley = 180;
    end
end

[vol_rr, P_rr] = volrotateP_r( vol_r, round(P_r), angley );

% Cross section
cs  = vol_rr( :, :, P_rr(3) );
% Select the component of the cross section that contains the rotation point
cs2 = double( bwselect( cs, P_rr(2), P_rr(1), 8 ) );
cs = cs.*cs2;
%figure, imshow(cs)
%cs = volsqueeze(cs);
%cs(:,:,end) = [];
%cs(:,:,1) = [];

bou = bwboundaries(cs);
b = bou{1,1};
CP = round(mean(b,1));
%{
u = [ 0 0 1 ];
v = N;

ruv = vrrotvec(u,v);
R = vrrotvec2mat(ruv);
%CP = ( R*[ CP - P_rr(1:2), 0]' )' + P;
w = size(cs,1);
l = size(cs,2);
h = size(cs,3);
[ x, y, z ] = ind2sub( [w,l,h], find( cs(:) ) );
%{
display_settings;
[ f, v ] = isosurface( vol_rr, 0.5 );
patch('Faces',f, 'Vertices', v, 'facecolor', [1 0 0], 'facealpha', 0.25,'edgecolor','none');
hold on
plot3(y,x,z*P_rr(3),'bo')
plot3( b(:,2), b(:,1), ones(size(b,1))*P_rr(3),'bo')
%}
%{
cs = [ x y z*P_rr(3) ];
cs = cs - ones(size(cs))*diag(P_rr);
cs = ( R*cs' )' + ones(size(cs))*diag(P);
% Calculate length of cross section
L = 0;
for j = 1 : size(b,1)-1
    L = L + norm( b(j+1,:) - b(j,:) );
end
%}

%}

nump = size(b,1);

dist = zeros(nump,1);
for i = 1 : nump
    dist(i) = norm(b(i,:) - CP);
end

%r = mean(dist);
r = max(dist);
%r = mean( [ mean(dist) max(dist) ] );