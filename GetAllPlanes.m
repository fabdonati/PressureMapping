function [Plane,BinaryOfCenterline] = GetAllPlanes(binarymaskfile,Points)
% FUnction to get the equation of the entry and exit planes of a vessel
% segment, providing:
% - Image of the mask of the vessel segment
% - (optional) 3D location of the points indicating location of the entry and exit
% - a is the most external layer and b is 1-layer inside the volume
bDebug = 0;

[Centerline, BinaryOfCenterline] = GetCenterline(binarymaskfile,bDebug);

if nargin>=2
    % We do not want the first and last point, but the ones closest to
    % Points:
    if numel(Points)~=6
        fprintf('ERROR! Points introduced in GetPlanes are not as expected (6 values = 2 points * 3 coordinates)\n');
        return;
    else
        [a b] = size(Points);
        if a>b, Points = Points'; end
        Point1 = Points(1,:);
        Point2 = Points(2,:);
        % TO BE IMPLEMENTED IN CenterlineSpline class:
        %P0 = Centerline.FindClosestPoint(Point1);
        %P1 = Centerline.FindClosestPoint(Point2);
    end
else
    for i = 1:Centerline.nPointsSpline
      P(i) = i;
      Point(i,:) = Centerline.GetPoint(P(i));
    end
end

% Estimate the planes:
h = exp(-10);
for i = 2:Centerline.nPointsSpline - 1
  Plane.slope(i,:) = normalise( Centerline.GetPoint(P(i)+h) - Centerline.GetPoint(P(i))) ;
  Plane.point(i,:) = Point(i,:);
end


if(bDebug)
  figure; hold on;
  [im,hd] = io_ReadMedicalImage(binarymaskfile);
%   subplot(211)
%   title('Skeleton output')
%   
%   show_segment_surface(im,hd.Mv2w);
%   show_segment_surface(BinaryOfCenterline,hd.Mv2w,[],0.5);

%  subplot(212); hold on;
  title('Centerline spline')
  Centerline.PlotSpline();
  show_segment_surface(im,hd.Mv2w);
  PlotPlane(EntryPlaneA.point,EntryPlaneA.slope,10); 
  PlotPlane(ExitPlaneA.point,ExitPlaneA.slope,10);
end

function [a] = normalise(b)
a = b./sqrt(sum(b.^2));

function [Centerline,BinaryOfCenterline] = GetCenterline(binarymaskfile,bDebug)
[im,hd] = io_ReadMedicalImage(binarymaskfile);

nrrdIM = Build_nrrd(im,hd);
BinaryOfCenterline = itk_imfilter('skel',nrrdIM);
%% ADV comment out/in
nrrdIM.data = BinaryOfCenterline; 
[F Lenght nPoints] = SmoothSkeleton(nrrdIM,bDebug);
nPoints = nPoints * 1;
Centerline = CenterlineSpline();
Centerline =  Centerline.SetSpline(F,nPoints,Lenght);

function [] = PlotPlane(point,slope,scale)

if nargin < 3
  scale = 1;
end

hold on;
plot3(point(1), point(2), point(3), '*');

vector = scale*slope;

quiver3(point(1), point(2), point(3), vector(1), vector(2), vector(3));

% Add a squared plane to be drawn as well..