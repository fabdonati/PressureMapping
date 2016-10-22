function [ Centerline, Point, Slope, nPoints ] = ... 
    get_centerline( im, hd, dCP, bInvertPoints, bExtraPoints, bDebug )
% FUnction to get the equation of the entry and exit planes of a vessel
% segment, providing:
% - Image of the mask of the vessel segment
% - (optional) 3D location of the points indicating location of the entry and exit
% - a is the most external layer and b is 1-layer inside the volume

if nargin < 4; bInvertPoints = 0; end
if nargin < 5; bExtraPoints  = 0; end
if nargin < 6; bDebug        = 0; end

[ Centerline, BinaryOfCenterline, nPoints ] = ... 
    GetCenterline( im, hd, dCP, bDebug, bInvertPoints, bExtraPoints );

Point = ( Centerline.GetAllPoints )';
Slope = ( Centerline.GetAllSlopes )';
Slope = normr(Slope);

if(bDebug)  
    figure
    
    subplot(211)
    hold on
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

    function [ Centerline, BinaryOfCenterline, nPoints ] = ... 
            GetCenterline( im, hd, dCP, bDebug, bInvertPoints, bExtraPoints )
        BinaryOfCenterline = Skeleton3D( im );
        nrrdIM = Build_nrrd( BinaryOfCenterline, hd );
        [ F, Length, nPoints, points ] = SmoothSkeleton( nrrdIM, bDebug, 0, bInvertPoints );
        nExtraPoints = 0; if bExtraPoints; nExtraPoints = 2; end
        nPoints = floor( Length / dCP ) + 1;
        Length = (nPoints-1)*dCP;
        param = -dCP*nExtraPoints : dCP : Length; % Initial point is usually 0
        points = fnval(F,param);
        nPoints = nPoints + nExtraPoints;
        Length = Length + dCP*nExtraPoints;
        param = 0 : dCP : Length;
        F = csaps( param, points, 1 );
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