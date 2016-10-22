function [ pp, Length ] = recalibrateCenterline( f, vT, Point, Slope, sm )

nPoints = size(Point,1);
for iP = 1 : nPoints
    progress = floor(iP/nPoints*100);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('%3d%% of centerline recalibrated', progress );
    [ ~, Point(iP,:) ] = mesh2crossSection( f, vT, Point(iP,:), Slope(iP,:), 0 );
end
PointDist = [0; cumsum( sqrt( sum( ( Point(2:end,:) - Point(1:end-1,:) ).^2, 2) ) ) ];
Length = PointDist(end);
pp = spaps( PointDist, Point', sm );