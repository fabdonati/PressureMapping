function [ Centerline, Plane, BinaryOfCenterline ] = ...
    get_planes_and_centerline( binarymaskfile )
% Function to get the equation of the entry and exit planes of a vessel
% segment, providing:
% - Image of the mask of the vessel segment
% - (optional) 3D location of the points indicating location of the entry and exit
% - a is the most external layer and b is 1-layer inside the volume

bDebug = 0;

[ Centerline, BinaryOfCenterline ] = GetCenterline( binarymaskfile, bDebug );

Point = ( Centerline.GetAllPoints )';
Slope = ( Centerline.GetAllSlopes )';
Slope = ((Slope')*diag( sqrt(sum(Slope.^2,2)) ))';

for i = 1 : Centerline.nPointsSpline
    Plane.slope(i,:) = Slope(i,:);
    Plane.point(i,:) = Point(i,:);
end

if(bDebug)
    [im,hd] = io_ReadMedicalImage(binarymaskfile);
    
    figure
    hold on
    subplot(211)
    title('Skeleton output')
    show_segment_surface(im,hd.Mv2w);
    show_segment_surface(BinaryOfCenterline,hd.Mv2w,[],0.5);
    
    subplot(212)
    hold on
    title('Centerline spline')
    Centerline.PlotSpline();
    show_segment_surface(im,hd.Mv2w);
    PlotPlane(EntryPlaneA.point,EntryPlaneA.slope,10);
    PlotPlane(ExitPlaneA.point,ExitPlaneA.slope,10);
end

    function [ Centerline, BinaryOfCenterline ] = GetCenterline( binarymaskfile, bDebug )
        [im,hd] = io_ReadMedicalImage( binarymaskfile );
        BinaryOfCenterline = Skeleton3D( im );
        nrrdIM = Build_nrrd( BinaryOfCenterline, hd );
        [ F, Length ] = SmoothSkeleton( nrrdIM, bDebug );
        nPoints = 100;
        Centerline = CenterlineSpline();
        Centerline = Centerline.SetSpline( F, nPoints, Length );
    end

    function PlotPlane( point, slope, scale )
        if nargin < 3
            scale = 1;
        end
        hold on
        plot3(point(1), point(2), point(3), '*');
        vector = scale*slope;
        quiver3( point(1), point(2), point(3), vector(1), vector(2), vector(3) );
    end
end