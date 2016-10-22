function [ F, Length, nPoints, points ] = ... 
    SmoothSkeleton( sknrrd, bDebugSplineFit, bCorrectScaleSkeleton, bInvertPoints )
% Function to calculate the spline coefficients of the line fitting the
% skeleton
% INPUT:
% - sknrrd: skeleton, a nrrd image with the voxels corresponding to the
% centerline of the skeleton.
%
% OUTPUT:
% - F: spline coefficients, see csaps.m
% - nPoints: number of voxels originally in the skeleton
%
% By Pablo Lamata, Oxford June 2011.
%
% Control version:
% 12/07/2011: control of the definition of the initial and final point of
% the line (final one: the closest to the bounding box).

if nargin<4
    bInvertPoints=0;
end

if nargin<3
    % Initial runs needed to change from metres to milimeters. Once the
    % step3 was corrected, this was not any longer needed.
    bCorrectScaleSkeleton=0;
end


if(nargin<2)
    bDebugSplineFit = 0;
end

if(bCorrectScaleSkeleton)
    Scale  =1000;
    fprintf('   Skeleton image is scaled by factor %i!\n',Scale);    
    for Iax = 1:3
        sknrrd.axis(Iax).spacing  = sknrrd.axis(Iax).spacing* Scale;
        sknrrd.axis(Iax).min      = sknrrd.axis(Iax).min    * Scale;
        sknrrd.axis(Iax).max      = sknrrd.axis(Iax).max    * Scale;
    end
end

%--------------------------------------------------------------------------
% Calculate the approximated spline through the Skeleton of the
% SkeletonImage (in voxel coordinates):
%--------------------------------------------------------------------------
% 1. Order the voxels of the skeleton:
% Robust Gerardus version (clpc446):
%[SkIm, cc] = skeleton_label(sknrrd.data,[],[sknrrd.axis.spacing],10/180*pi);
% Updated Gerardus version (laptop):
[SkIm, cc] = skeleton_label(sknrrd.data,sknrrd.data,[sknrrd.axis.spacing],10/180*pi);
%--------------------------------------------------------------------------
% 2. Generate the 3D coordinates of the ordered points:
IndexBranchOfInterest = find(cc.BranchLength==max(cc.BranchLength));
LengthNotConsidered = (sum(cc.BranchLength) - cc.BranchLength(IndexBranchOfInterest))/sum(cc.BranchLength);
if LengthNotConsidered>0.05
    %fprintf('WARNING! there was a %1.1f %% of skeleton not included!\n',100*LengthNotConsidered);
    bJoinSecondBranch = 1;
else
    bJoinSecondBranch = 0;
end
[r, c, s] = ind2sub(size(SkIm), cc.PixelIdxList{IndexBranchOfInterest});
RowColumSliceCoordinats = [reshape(r,numel(r),1),reshape(c,numel(c),1),reshape(s,numel(s),1)];
points = scinrrd_index2world(RowColumSliceCoordinats,sknrrd.axis);
param  = cc.PixelParam{IndexBranchOfInterest};
if(bJoinSecondBranch)
    [ ~, I ] = sort(cc.BranchLength);
    IndexSecondBranch = I(end-1);
    [r, c, s] = ind2sub(size(SkIm), cc.PixelIdxList{IndexSecondBranch});
    RowColumSliceCoordinats = [reshape(r,numel(r),1),reshape(c,numel(c),1),reshape(s,numel(s),1)];
    points2 = scinrrd_index2world(RowColumSliceCoordinats,sknrrd.axis);
    P11 = points(1,:);    P12 = points(end,:);
    P21 = points2(1,:);   P22 = points2(end,:);
    % Study the concatenation order:
    D(1) = sum(P11-P21);
    D(2) = sum(P11-P22);
    D(3) = sum(P12-P21);
    D(4) = sum(P12-P22);
    D=abs(D);
    I = find(D==min(D));
    switch I
        case 1
            % Need to invert the order of points 1
            points = invertPoints(points);            
        case 2
            % Need to invert the order of points 1 and 2
            points = invertPoints(points); 
            points2 = invertPoints(points2); 
        case 3
            % Do nothing, this is the rigth order
        case 4
            points2 = invertPoints(points2); 
    end
    param2  = cc.PixelParam{IndexSecondBranch};
    points = [points;points2];
    param = [param param2+max(param)];
    LengthNotConsidered = (sum(cc.BranchLength) - max(param)/sum(cc.BranchLength));
    if LengthNotConsidered>0.05
        fprintf('WARNING! there was a %1.1f %% of skeleton not included! (despite joining two branches) \n',100*LengthNotConsidered);
    end
end    
%--------------------------------------------------------------------------
% 2.B: reorder the points, defining the last to be the closest to the
% bounding box:
P1 = points(1,:);
P2 = points(end,:);
C1 = [sknrrd.axis(1).min sknrrd.axis(2).min sknrrd.axis(3).min];
C2 = [sknrrd.axis(1).max sknrrd.axis(2).max sknrrd.axis(3).max];
%Corner1 = scinrrd_index2world(RowColumSliceCoordinats,sknrrd.axis);
permuteOrder = [2 1 3];
BB = [C1(permuteOrder(1)) C2(permuteOrder(1)) C1(permuteOrder(2)) C2(permuteOrder(2)) C1(permuteOrder(3)) C2(permuteOrder(3))];
DistanceOptions.bYZ=1;
D1 = GetMinimumDistance2BoundingBox(P1,BB,DistanceOptions);
D2 = GetMinimumDistance2BoundingBox(P2,BB,DistanceOptions);
if D2>D1
    fprintf('Inversion of points! in order to have point 1 in ascending aorta\n');
    points = invertPoints(points); 
end

if(bInvertPoints)
    points = invertPoints(points); 
end

%--------------------------------------------------------------------------
% 3. Fit the spline (different options were analyised, 2 and 3 are best)
splineoption = 3;
nPoints = numel(points)/3;
[a,b]=size(points);
if a>b, points = points'; end
m0= 1; m = numel(points)/3;     
switch splineoption
    case 1                
        F = spline((1:m),points);
    case 2
        % An empirical rule for the number of spline segments: one
        % every 6 points
        nPointsForEachSegment = 6;
        pieces = ceil(nPoints/nPointsForEachSegment);
        F = splinefit((1:m),points,pieces);
    case 3
        m0 = min(param);
         m = max(param);
        % An empirical rule for the smoothness, the choice of h:
        h = 10*mean([sknrrd.axis.spacing]);
        P = 1/(1+(h^3)/6);
        F = csaps( param, points, P*10 );
    otherwise
        fprintf('ERROR! Not valid option for spline approximation\n');
end

% Get Length:
Length = max(param);

% Trajectory
if(bDebugSplineFit)
    step = (m-m0)/(3*nPoints);
    t = m0 : step : m;
    Ft = ppval(F,t);
    hold on
    plot3( Ft(1,:), Ft(2,:), Ft(3,:), '-g', 'LineWidth', 3 );        
    plot3( points(1,:), points(2,:), points(3,:), 'b.', 'MarkerSize', 20 );
    axis equal
    drawnow
end

    function points = invertPoints( points )
        points = points(end:-1:1,:);
    end
end