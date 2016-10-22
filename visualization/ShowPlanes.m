function [] = ShowPlanes(BinaryImageName,planes)

    %figure('Color',[1 1 1]);    
    hold on;
    if nargin>1
        [a nPlanes]= size(planes);
        for iPlane = 1:nPlanes
            points = planes(iPlane).plane.Plane(:,1:3);
            %pointB = planes(iPlane).plane.Plane(2,1:3);
            %pointC = planes(iPlane).plane.Plane(3,1:3);
            points(4,:) = points(1,:) + (points(3,:)-points(2,:));
            Centre = mean(points,1);
            points(5,:) = points(1,:);
            plot3(points(:,1), points(:,2), points(:,3),'ro-');
            %plot3(Origin1(1),Origin1(2),Origin1(3),'ko','LineWidth',5)
            %plot3(Origin2(1),Origin2(2),Origin2(3),'go','LineWidth',5)
            %plot3(Origin2(1)-CROP(1),Origin2(2)-CROP(2),Origin2(3)-CROP(3),'bo','LineWidth',5)
            %plot3(SC(1)*points(:,PO(1))+OF(1), SC(2)*points(:,PO(2))+OF(2), SC(3)*points(:,PO(3))+OF(3),'ro-');
            %plot3(pointB(2),pointB(1),pointB(3),'go');
            %plot3(pointC(2),pointC(1),pointC(3),'co');
            text(Centre(1),Centre(2),Centre(3),planes(iPlane).plane.Description);
        end
    end
    if exist(BinaryImageName,'file')
        ShowBinaryMask(BinaryImageName);
    else
        fprintf('warning: binary image does not exist (%s)\n',BinaryImageName);
    end
