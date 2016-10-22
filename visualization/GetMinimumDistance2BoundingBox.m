function [Distance] = GetMinimumDistance2BoundingBox(point,BoundingBox,options)
% Function to calculate the minimum distance of a point inside an image to
% the bounding box.
%
% INPUT:
% - point: (x,y,z) coordinates
% - BoundingBox: [Xmin Xmax Ymin Ymax Zmin Zmax]
% - options: choice of which plane to compare
%     options.bYZ: only the YZ planes
%
% OUTPUT:
% - Distance
%
% By Pablo Lamata, Oxford 12 July 2011

if nargin==3
    if isfield(options,'bYZ')
        bYZ = options.bYZ;
    else
        bYZ = 0;
    end
end

% The indexes of the BoundingBox elements that define each of the six 
% planes, i.e. plane XY
%PlaneCoordinate1 = [1 2 1 2 3 4];
%PlaneCoordinate2 = [3 4 5 6 5 6];
OutOfPlaneCoord  = [5 6 3 4 1 2];
WhichPlanes      = [1 1 1 1 1 1];
if bYZ
    WhichPlanes      = [0 0 0 0 1 1];
end
Distances = zeros(1,6);
for iPlane = 1:6
    if(WhichPlanes(iPlane)) % If this is a plane selected for comparison
        iCoor = OutOfPlaneCoord(iPlane);
        B = BoundingBox(iCoor);
        switch iCoor
            case {1,2}, A = point(1); % This is an X coordinate
            case {3,4}, A = point(2);
            case {5,6}, A = point(3);
        end
        Distances(iPlane) = abs(A-B);
%         iCoor(1) = PlaneCoordinate1(iPlane);
%         iCoor(2) = PlaneCoordinate2(iPlane);
%         for iC=1:2
%             B = BoundingBox(iCoor(iC));
%             switch iCoor(iC)
%                 case {1,2}, A = point(1); % This is an X coordinate
%                 case {3,4}, A = point(2);
%                 case {5,6}, A = point(3);
%             end
%             Distances(iPlane) = Distances(iPlane) + (A-B)^2;
%         end    
    else
        Distances(iPlane) = Inf;
    end
end
Distance = min(Distances);