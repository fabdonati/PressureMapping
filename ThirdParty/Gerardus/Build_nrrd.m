function nrrd = Build_nrrd(im,imh)
% Function to build the nrrd image format from the matrix and header format
% Key references:
% - read_image_vtk2.m and parseheader.m (adaptation of D.C.Barber's)
% - scinrrd_load.m (from Gerardus)
%
% The main issue is to change from the x,y,z indexing in images, to the
% row, column, slice indexing in matlab:

bGerardus = 1;
PermuteOrder = [2 1 3];
if(bGerardus)
    % Call Ramon's code
    im   = permute(im,PermuteOrder);
    res(PermuteOrder) = imh.spacing;
    offset(PermuteOrder) = imh.origin ;
    %offset = offset - res;
    RotationMatrix = imh.TransformMatrix;
    nrrd = scinrrd_im2nrrd(im, res, offset, RotationMatrix);    
else    
    
    nrrd.data   = permute(im,PermuteOrder);
    imspacing(PermuteOrder) = imh.spacing;
    imsize(PermuteOrder)    = size(im);
    AxisMin(PermuteOrder)   = imh.origin;
    AxisMax(PermuteOrder)   = reshape(imh.origin,1,3) + reshape((size(im) - [1 1 1]).*reshape(imspacing(PermuteOrder),1,3),1,3);
    center      = (AxisMax-AxisMin)/2;
    for I = 1:3
        nrrd.axis(I,1).size = imsize(I);
        nrrd.axis(I,1).spacing = imspacing(I);
        % The min and max are the left corner of the first and last voxels!:
        nrrd.axis(I,1).min = AxisMin(I) - imspacing(I)/2;
        nrrd.axis(I,1).max = AxisMax(I) - imspacing(I)/2;
        nrrd.axis(I,1).min = AxisMin(I);
        nrrd.axis(I,1).max = AxisMax(I);
        nrrd.axis(I,1).center = center(I);
        nrrd.axis(I,1).unit = 'no unit';
    end
    nrrd.axis(1).label = 'axis 2';
    nrrd.axis(2).label = 'axis 1';
    nrrd.axis(3).label = 'axis 3';
end
