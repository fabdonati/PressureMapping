function [ A, hdA ] = allignAnatBinary( A, hdA, R, hdR, shift )

if nargin < 5
    shift = [ 0 0 0 ];
end

% function allignAnatBinary( aVTK, rVTK, shift )
% 
% if nargin < 3
%     shift = [ 0 0 0 ];
% end

% rVTK = 'C:\Users\Alessandro\OneDrive\Matlab\PressureMapping\Data\INPUT\HLHS\HT002\sys_seg_ref.vtk';
% aVTK = 'C:\Users\Alessandro\OneDrive\Matlab\PressureMapping\Data\Nidhin\HT002.vtk';
%[ imR.b, hdR ] = io_ReadMedicalImage( rVTK );
%[ imA.b, hdA ] = io_ReadMedicalImage( aVTK );

A( A~=1 ) = 0;

w = size(A,1);
l = size(A,2);
h = size(A,3);
[ x, y, z ] = ind2sub( [w,l,h], find( A(:) ) );

C = [ x y z ];
C = C - ones( size(C) );
C = C(:,[2 1 3]);
C = [ C ones(size(C,1),1) ] * hdA.Mv2w';

C = [ C ones( size(C,1),1) ] * hdR.Mw2v';
C = C(:,[2 1 3]);
C = C + ones( size(C) );

C = round( C );
C(:,2) = C(:,2) + shift(2);
C(:,3) = C(:,3) + shift(3);

B = zeros( size(R) );
for i = 1 : size(C,1); B( C(i,1), C(i,2), C(i,3) ) = 1; end

A = B;
hdA.origin  = hdR.origin;
hdA.spacing = hdR.spacing;
hdA.Mv2w    = hdR.Mv2w;
hdA.Mw2v    = hdR.Mw2v;