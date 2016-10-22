function [EnsightData] = io_ReadEnsightDirectory(dirname,rootname)
% Funciton to load the case, and associated data, in an Ensight folder
% INPUT:
% - dirname: directory with the data
% - rootname: the common root name in all files

bVel = 1;
bMag = 1;

if nargin<2
    % In Pablo's machine, this is the most common root name:
    rootname = 'OpenCMISS';
end


GeoFile = [rootname '.geo'];
[Imageheader,points] = ReadGeoFile(dirname,GeoFile,1);
CaseInfo = ReadCaseFile(dirname,rootname);
if(bVel)
    Velocity = NaN*ones([Imageheader.dim(1:3), CaseInfo.nSteps,3]);
end
if(bMag)
    Magnitude = NaN*ones([Imageheader.dim(1:3), CaseInfo.nSteps]);
end
for iFrame = 1:CaseInfo.nSteps
    if(bVel)
        fileName = sprintf('%s.Velocity%04i',rootname,iFrame-1);
        opt.dimensionOfValues = 3;
        [values,ReadResult,Mag] = ReadEnsightValues(dirname,fileName,Imageheader,opt);
        for ic = 1:3
            Velocity(:,:,:,iFrame,ic) = squeeze(values(:,:,:,ic));
        end
    end
    if(bMag)
        fileName = sprintf('%s.Magnitude%04i',rootname,iFrame-1);
        opt.dimensionOfValues = 1;
        [values,ReadResult,Mag] = ReadEnsightValues(dirname,fileName,Imageheader,opt);
        for ic = 1:3
            Magnitude(:,:,:,iFrame) = values(:,:,:);
        end
    end
    clear 'values'
end

EnsightData.CaseInfo = CaseInfo;
if(bVel)
    EnsightData.Velocity = Velocity;
end
if(bMag)
    EnsightData.Magnitude = Magnitude;    
end
EnsightData.Imageheader = Imageheader;