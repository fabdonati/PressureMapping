function [ Centerline, BinaryOfCenterline, pp, Length, nPoints ] = ...
                                            GetCenterline( binarymaskfile, bDebug, ro, co, he )

if nargin<2
    bDebug = 0;
end

[im,hd] = io_ReadMedicalImage( binarymaskfile );

im           = im(ro,co,he);
hd.origin    = hd.origin - [ ro(1)-1  co(1)-1 he(1)-1 ];
hd.dim       = size( im );
hd.points    = prod( hd.dim );
hd.Mv2w(:,4) =  hd.origin';
hd.Mw2v = GetMw2v(hd.Mv2w);

% Old version depending on Gerardus:
% nrrdIM = Build_nrrd(im,hd);
% BinaryOfCenterline = itk_imfilter('skel',nrrdIM);
% nrrdIM.data = BinaryOfCenterline; 
% New version depending on Matlab Central code:
% http://uk.mathworks.com/matlabcentral/fileexchange/43400-skeleton3d
BinaryOfCenterline = Skeleton3D(im);

%% ADV
% io = input('Do you wish to add a correction point to the centreline? [0=no;1=yes]  ');
% if io==1
%     cp = input('Correction point coordinates (e.g. [x y z]:  ');%[60 55 16]
%     BinaryOfCenterline2 = BinaryOfCenterline;
%     BinaryOfCenterline(cp(1),cp(2),cp(3)) = 1;
% end
% if(1)
%     show_segment_surface(BinaryOfCenterline);
%     hold on
%     show_segment_surface(BinaryOfCenterline2);
% end
nrrdIM = Build_nrrd(BinaryOfCenterline,hd);
[ pp, Length, nPoints ] = SmoothSkeleton( nrrdIM, bDebug) ;
nPoints = round( nPoints * 1 );
Centerline = CenterlineSpline();
Centerline = Centerline.SetSpline( pp, nPoints, Length );
