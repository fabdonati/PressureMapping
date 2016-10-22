function [Angle,H] = GetInflow2SeptumAngle(RootDir,iCase)

% Script to compute the angle between a flow jet through a plane, and
% another plane.

bViewSeptumVelocities = 1;

Pname1 = 'Plane_entry';
Pname2 = 'Plane_IAS';

DirName = [RootDir sprintf('Case%i/',iCase)];
CaseDir = [DirName 'case/'];
ContextDir = [DirName 'context/'];

RootNam{1} = sprintf('PFO%03i_case',iCase);
RootNam{2} = sprintf('PFO%03i',iCase);
RootNam{3} = 'case';

% Find out which is the correct naming:
for ii = 1:3
    RN = RootNam{ii};
    [Imageheader] = ReadGeoFile(CaseDir,RN,1,0);
    if isstruct(Imageheader)
        RootName = RootNam{ii};
        break;
    end
end

contextfile = sprintf('PFO%03i.ctx',iCase);
VelocityFile = [RootName '.fVelocity0004'];
MaskFile = [RootName '.mask0004'];

bRenderAngle = 1;

%% Read the data
% Read the planes:
planes = ReadPlanesFromEnsightContextFile([ContextDir contextfile]);
EntryPlane = GetPlaneFromName(planes,Pname1);
SeptumPlane = GetPlaneFromName(planes,Pname2);

% Read the velocity data, need to read slowly due to need of rotation
% matrix:
[Imageheader] = ReadGeoFile(CaseDir,RootName,1,0);
options.dimensionOfValues = 3;
[imVel] = ReadEnsightValues(CaseDir,VelocityFile,Imageheader,options);
%VelMag = sqrt(sum(im.^2,4));
% Mask = VelMag;
% Threshold = 0.3;
% Mask(Mask>=Threshold) = 1;
% Mask(Mask<Threshold) = 0;
Mask = ReadEnsightValues(CaseDir,MaskFile,Imageheader,1);

%% Compute the angle
% Compute the vector of the flow jet through plane 1
MaskNrrd = Build_nrrd(Mask,Imageheader);
V1nrrd =  Build_nrrd(squeeze(imVel(:,:,:,1)),Imageheader);
V2nrrd =  Build_nrrd(squeeze(imVel(:,:,:,2)),Imageheader);
V3nrrd =  Build_nrrd(squeeze(imVel(:,:,:,3)),Imageheader);

Centre = GetPlaneCentrePoint(EntryPlane);
normal = GetPlaneNormal(EntryPlane); 
[Mask2D gx gy gz midx] = scinrrd_intersect_plane(MaskNrrd,Centre,normal,'linear');   

[V1 gx gy gz midx] = scinrrd_intersect_plane(V1nrrd,Centre,normal,'linear');   
[V2 gx gy gz midx] = scinrrd_intersect_plane(V2nrrd,Centre,normal,'linear');   
[V3 gx gy gz midx] = scinrrd_intersect_plane(V3nrrd,Centre,normal,'linear'); 

% Visualise velocity vectors at the septum plane:
CentreSP = GetPlaneCentrePoint(SeptumPlane);
normalSP = GetPlaneNormal(SeptumPlane);
[V1SP gxSP gySP gzSP midxSP] = scinrrd_intersect_plane(V1nrrd,CentreSP,normalSP,'linear');   
[V2SP gxSP gySP gzSP midxSP] = scinrrd_intersect_plane(V2nrrd,CentreSP,normalSP,'linear');   
[V3SP gxSP gySP gzSP midxSP] = scinrrd_intersect_plane(V3nrrd,CentreSP,normalSP,'linear'); 


V1 = V1.*Mask2D;
V2 = V2.*Mask2D;
V3 = V3.*Mask2D;
% Valid points those inside the mask:
Ivalid1 = find(~isnan(Mask2D));
% And also those pionts within a certain distance of the mid point (avoid
% other vessels)
DistThreshold = 15; % in mm!
Distances = sqrt( (gx - gx(midx(1),midx(2))).^2 + (gy - gy(midx(1),midx(2))).^2  + (gz - gz(midx(1),midx(2))).^2 );
Ivalid2 = find(Distances<DistThreshold);
Ivalid = intersect(Ivalid1,Ivalid2);
Vjet = [mean(V1(Ivalid)) mean(V2(Ivalid)) mean(V3(Ivalid))];
Vjet = Vjet / sqrt(sum(Vjet.^2));

% Compute the angle between the two directions
CentreS = GetPlaneCentrePoint(SeptumPlane);
normalS = GetPlaneNormal(SeptumPlane);
Angle = 180 * acos( dot(Vjet,normalS) ) / pi;

if(bViewSeptumVelocities)
    figure('color',[1 1 1])
    ShowPlanes('',planes);
    quiver3(gxSP(Ivalid), gySP(Ivalid), gzSP(Ivalid), V1SP(Ivalid), V2SP(Ivalid), V3SP(Ivalid), 10);
    
    
end


%% Rendering
if (bRenderAngle)
    P2 = Centre + 15*normal;
    PJ = Centre + 15*Vjet;
    PS = CentreS + 15*normalS;
    Hfig = figure('color',[1 1 1],'OuterPosition',[100 100 1000 1000]);
    subplot(221)
    ShowPlanes('',planes);    
    show_segment_surface(Mask,Imageheader.Mv2w);
    title('Mask and planes'); view(90,0)
%     subplot(222)
%     ShowPlanes('',planes); hold on;
%     plot3(Centre(1),Centre(2),Centre(3),'g*');
%     plot3([Centre(1) P2(1)],[Centre(2) P2(2)],[Centre(3) P2(3)],'g');
%     plot3([CentreS(1) PS(1)],[CentreS(2) PS(2)],[CentreS(3) PS(3)],'b');
%     title('Planes, entry in green, septum in blue'); view(90,0); axis equal;
    subplot(222) 
    ShowPlanes('',planes);
    surf(gx, gy, gz, Mask2D, 'EdgeColor', 'none'); hold on
    plot3([Centre(1) P2(1)],[Centre(2) P2(2)],[Centre(3) P2(3)],'g');
    show_segment_surface(Mask,Imageheader.Mv2w);
    view(90,0); axis equal;
    title('Isosurface with 2D plane drawing mask');
    
    subplot(223)
    ShowPlanes('',planes);
    quiver3(gx(Ivalid), gy(Ivalid), gz(Ivalid), V1(Ivalid), V2(Ivalid), V3(Ivalid), 10);
    view(90,0); axis equal;
    title('Entry yet');
    
    subplot(224)
    ShowPlanes('',planes);
    plot3([CentreS(1) PS(1)],[CentreS(2) PS(2)],[CentreS(3) PS(3)],'b','LineWidth',3);
    plot3([Centre(1) PJ(1)],[Centre(2) PJ(2)],[Centre(3) PJ(3)],'b','LineWidth',3);
    view(90,0); axis equal; 
    title(sprintf('Case %i: angle is %1.2f',iCase,Angle));
    FigName = sprintf('%sCase%i.png',DirName,iCase);
    export_fig(FigName,'-png','-m2',Hfig);
else
    H = NaN;
end

function Plane = GetPlaneFromName(planes,Pname1)
Plane = [];
for iP = 1:numel(planes)
    p = planes(iP).plane;
    name = p.Description;
    names(iP,1:numel(name)) = name;
    if strcmp(name,Pname1)
        Plane = p;
        return;
    end
end
fprintf('error! Plane named %s not found in the list provided:\n',Pname1);
fprintf('    %s\n',names);

function [Centre,normal] = GetPlaneCentrePoint(plane)
    points = plane.Plane(:,1:3);
    points(4,:) = points(1,:) + (points(3,:)-points(2,:));
    Centre = mean(points,1);
function [normal] = GetPlaneNormal(plane)    
    points = plane.Plane(:,1:3);
    V1 = points(2,:) - points(1,:);
    V2 = points(3,:) - points(1,:);
    normal = cross(V1 , V2);
    normal = normal ./ sqrt(sum(normal.^2));