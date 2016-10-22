% Script to test the code to find an intersection plane

TestImage = 'TestData\aorta2.vtk';

% Read the image:
[im, hd] = io_ReadMedicalImage(TestImage);
% Build the NRRD format image:
Pnrrd = Build_nrrd(im,hd);
% Get the intersection:
point = hd.origin + (hd.dim .* hd.spacing)/2;
slope = [1 1 1];
options.interp = 'linear';

iP = 0;
for InPlaneAngle = 0:15:75
    options.InPlaneAngle = pi*InPlaneAngle/180;
    [tmp, gx, gy, gz, midx] = scinrrd_intersect_plane(Pnrrd,point,slope,options);
    iP = iP+1;
    subplot(2,3,iP)
    imagesc(tmp); axis equal;
    %surf(gx, gy, gz, tmp, 'EdgeColor', 'none'); 
end