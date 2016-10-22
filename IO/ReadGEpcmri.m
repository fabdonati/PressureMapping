function [hd,Mag,Phase] = ReadGEpcmri(DicomDir,MagFolder,PhaseFolders,nSlices,nT)
    % Function to read PCMRI data, velocity values from GE format

    % Input: where the data is located. A dicom directory where individual
    % directories with the magnitude and velocity data are:
%     DicomDir = 'F:\PPE\PAH-15T-4Dflow\PAH-15T-4Dflow\FU_YUN_PAH-320\22439\';
%     MagFolder = '434';
%     PhaseFolders(1,:) = '531';
%     PhaseFolders(2,:) = '532';
%     PhaseFolders(3,:) = '533';
%     nSlices = 16;
%     nT = 20;

    %1. read the data:
    MagDir = fullfile(DicomDir,MagFolder);
    [hd,Mag] = GetImage(MagDir,nSlices,nT);
    
    Phase = zeros(size(Mag,1),size(Mag,2),size(Mag,3),size(Mag,4),3);
    for iC = 1:3
        PhaDir = fullfile(DicomDir,PhaseFolders(iC,:));
        [hd,ph] = GetImage(PhaDir,nSlices,nT);
        Phase(:,:,:,:,iC) = ph;
    end
    
%     % Compute the phase angiography, averaged in time:
%     PhaseMag = sqrt(sum(Phase.^2,5));
%     MeanPhaseMag = mean(PhaseMag,4);
% 
% 
%     % TODO: get this information from the GE header
%     opts.dt = 40; %CaseInfo.deltaT/1000;
%     cd(CaseDir);
%     [werp] = WerpPressureEstimation(BinaryMask,Velocity,opts);

function [hd,Image] = GetImage(directory,nSlices,nT)
    orderconvention = 'TimeSlice';
    ListOfFiles = dir(directory);
    for iFile = 3:numel(ListOfFiles)
        dicomfile = fullfile(directory,ListOfFiles(iFile).name);
        [im,hd] = dicomread(dicomfile);
        if iFile == 3
            Zmax = 300;
            Image = zeros(size(im,1),size(im,2),Zmax);
        end
        iSlice = iFile - 2;
        Image(:,:,iSlice) = im;
    end
    Image(:,:,iFile+1:end) = [];
    dicomH = dicominfo(dicomfile);
    % Reorder the data accordingly to the temporal frames
    temp = Image;
    Image = zeros(size(temp,1),size(temp,2),nSlices,nT);
    switch orderconvention
        case 'TimeSlice',
            for iT = 1:nT
                i1 = iT;
                D = nT;
                i2 = nT*nSlices;
                Image(:,:,:,iT) = temp(:,:,i1:D:i2);
            end
        case 'SliceTime', 
                    % TODO!
                    i1 = iT;
                        
        
    end
    % The header only has information of one slice:
    hd.Dimension = size(Image);
    hd = ParseHeader(dicomH);
