% Script to prepare the data for Fabrizio:

clear;
close all;
fclose all;

bSaveVelocity = 1;
bSaveMagnitude = 0;
bSaveImAverageVel= 0;
bSaveAnatomicalLocations = 0;
bSaveResults = 0;
bQualityAssessment = 0;
bSaveHR = 0;

bSavePCMRexplorer = 0;
bDebug = 1;

cohortID = 4;
DestinyRootDir = 'F:\PPE\FabrizioCases';
%DestinyRootDir = 'C:\Data\FabrizioCases';
cohort =  'Alex';
bOriginalData = 0;
switch cohortID
    case 1
        % Healthy cases: 
        iStudy = 12;
        [cases,RootName,bNotAlex,nDigitsName,optionsCalculateComponents,subfolder] = PPEcaseMapping(iStudy);
        AdditionalName = '';
    case 2
        % Disease cases:
        cases = [1 2 79];
        AdditionalName = '';
    case 3
        % BAV cases
        cases = 1:278;
        cohort = 'BAV';     
        AdditionalName = 'BAV';
    case 4
        % BAV cases
%        cases = [017, 213, 106, 238, 056, 009, 051, 181, 176, 130, 22, 71, 87, 102, 8, 138, 139, 143, 154, 159, 161, 174, 192, 204, 220, 239, 244, 246, 263, 290];
%        cases = [];
        cases = [213];
        %cases = [71];
        %cases = 238;
        alreadyhad = [17, 56, 71, 106, 154, 238, 246,];
        %cases = [21, 37, 83, 116, 120, 247, 
        %
            cases = [ 172, 252,  256];

        cohort = 'BAV';     
        AdditionalName = 'BAV20';
        bOriginalData = 1;
end
if(bSaveHR)
    HR = [];
end
for iCase = 1:numel(cases)
    CaseID = cases(iCase);
    Directory = RetrieveCaseFolder(CaseID,cohort);
    paths = GetPathsPPEdata();
    Destiny = fullfile(DestinyRootDir,sprintf('Case%s%03i',AdditionalName,CaseID));
    
    % Read the original velocity data before the cropping:
    if(bOriginalData)
        VelDir = fullfile(Directory,paths.Step02out);
    else
        % Get the cropped domain obtained through Krittian's workflow:
        VelDir = fullfile(Directory,paths.Step04out);           
    end
    MaskFile = fullfile(Directory,paths.BinaryImWithPath); 
    nameData = 'Data.mat';
    nameRegional = 'RegionalAnalysis.mat';
    if(bSaveVelocity)||bSaveMagnitude
        if exist(fullfile(VelDir,'OpenCMISS.geo'),'file');
            if ~exist(Destiny,'dir'), mkdir(Destiny); end            
            [hd,ImagePoints,bError]=ReadGeoFile(VelDir,'OpenCMISS.geo',1);
            if (bError)
                % Try to get the correct GEO from the original binary
                % format
                 BinaryDir = ls(fullfile(Directory,'..','Split files','EnSight*'));
                 BinaryDir = fullfile(Directory,'..','Split files',BinaryDir);
                 GeoBinary = ls(fullfile(BinaryDir,'*.geo'));
                 [hd,ImagePoints,bError2] = ReadGeoFile(BinaryDir,GeoBinary,1);
            end
            CaseInfo = ReadCaseFile(VelDir,'OpenCMISS');
            if(bSaveMagnitude)
                dimensionOfValues = 1;
                Magnitude = NaN * ones(hd.dim(1),hd.dim(2),hd.dim(3),dimensionOfValues,CaseInfo.nSteps);
                BinaryName = 'mag';
                AsciiName = 'Magnitude';
            end
            if(bSaveVelocity)
                dimensionOfValues = 3;
                Velocity = NaN * ones(hd.dim(1),hd.dim(2),hd.dim(3),dimensionOfValues,CaseInfo.nSteps);
                BinaryName = 'vel';
                AsciiName = 'Velocity';
            end
            for iF = 1:CaseInfo.nSteps
                iFrameNumber = iF-1;
                if(bError)
                    VelDir = BinaryDir;
                    fileName = ls(fullfile(VelDir,sprintf('EnSight*%02i.%s',iFrameNumber,BinaryName)));
                    optionsReadEnsightValues.bBinaryFormat = 1;
                else
                    fileName = sprintf('OpenCMISS.%s%04i',AsciiName,iFrameNumber);
                end
                optionsReadEnsightValues.dimensionOfValues = dimensionOfValues;
                data = ReadEnsightValues(VelDir,fileName,hd,optionsReadEnsightValues);
                Velocity(:,:,:,:,iF) = data;
                %if iF==5
                    if (bSaveMagnitude) && bSavePCMRexplorer
                        namefile = fullfile(Destiny,sprintf('Mag%03i.vtk',iFrameNumber));
                        write_image_vtk2(namefile,data,hd,'ascii');
                    end          
                    fileNamePeakVel = sprintf('OpenCMISS.Velocity%04i',iFrameNumber);
                    opt2.dimensionOfValues = 3;
                    data = ReadEnsightValues(VelDir,fileNamePeakVel,hd,opt2);
                    VelMag = squeeze(sum(data.^2,4));
                    if(bSavePCMRexplorer)
                        namefile = fullfile(Destiny,sprintf('VelMag%03i.vtk',iFrameNumber));
                        write_image_vtk2(namefile,VelMag,hd,'ascii');
                    end
                %end
            end
            if(bSaveVelocity)
                save(fullfile(Destiny,nameData),'Velocity','hd','CaseInfo','ImagePoints');
            end
        end
    end    
    if (bSaveImAverageVel)
        file2save = fullfile(Directory,paths.ImageAverageVelocity);
        if exist(file2save,'file')
            copyfile(file2save,Destiny);
        else
            fprintf('WARNING! No image with average velocity available in case %s\n',Directory);
        end
        file2save = fullfile(Directory,paths.Path2mask,paths.BinaryImName);
        if exist(file2save,'file')
            copyfile(file2save,Destiny);
        else
            fprintf('WARNING! No segmentation available in case %s\n',Directory);
        end
    end
       
    if (bSaveAnatomicalLocations)
        % Anatomical locations:
        optionsSpline.bBetweenAnatomicalPlanes = 1;
        optionsSpline.bCropAtFixedDistance =  1;    
        [iRef,iInitAscending,iEndDescending,distanceBtw,iAnatomical] = GetCharacteristicSplinePoints(CaseID,optionsSpline);
        [~,~,~,~,~,Centerline] = LoadSplineFile(CaseID);
        save(fullfile(Destiny,nameData),'iAnatomical','Centerline','-append');
        if(bDebug)
            H = figure('color',[1 1 1]);
            [im,hd] = io_ReadMedicalImage(MaskFile);
            show_segment_surface(im,hd.Mv2w);
            hold on
            Centerline.PlotSpline();
            for iP = 1:numel(iAnatomical)
                Centerline.PlotSplinePoint(iAnatomical(iP));
            end
            view(90,90)
            title(sprintf('Case %i',CaseID));
            export_fig(fullfile(Destiny,'Anatomy.png'),'-png',H);
            close(H);
        end
        copyfile(MaskFile, fullfile(Destiny,'mask.vtk'));
    end
    if(bSaveResults)
        % Save the pressure drops in each segment, total and the components
        OutDir = fullfile(Directory,paths.Path2output);        
        resultsFile = fullfile(OutDir,nameRegional);
        if exist(resultsFile,'file')
            copyfile(resultsFile, fullfile(Destiny,nameRegional));
        else
            fprintf('Error, no regional analysis in %s \n',OutDir);
        end
    end
    if(bSaveHR)
        % Save the Heart Rythm
        % point to one directory with original dicoms:
        DicomDir = fullfile(Directory,'..','Split files','Mag');
        % Get one of the files:
        if ~exist(DicomDir,'dir')
            fprintf('Error, no dicom directory in %s\n',DicomDir);
        else
            listfiles = dir(DicomDir);
            iFile = 3; % first two files are'.' and '..'
            bHRfound = 0;
            %iValidDicom = 0;
            while iFile < numel(listfiles) && ~bHRfound            
                filename = listfiles(iFile).name;
                %check if this is a dcm:
                if strcmp(filename(end-3:end),'.dcm')
                    hd = dicominfo(fullfile(DicomDir,filename));
                    %iValidDicom =  iValidDicom + 1;
                    %HRvalues(iValidDicom) = hd.NominalInterval;
                    HR(iCase).value = hd.NominalInterval;
                    HR(iCase).CaseID = CaseID;
                    bHRfound = 1;
                end
                iFile = iFile + 1;
            end
        end
    end
    if(bQualityAssessment)
        load(fullfile(Destiny,nameRegional));
        load(fullfile(Destiny,nameData));
        % check that the anatomical locations are the same:
        iAnat1 = iAnatomical;
        iAnat2 = RegionalAnalysis.iAnatomical;
        if sum(iAnat1 - iAnat2) ~= 0
            fprintf('ERROR! different anatomical locations in case %i\n',CaseID);
        else
            fprintf('Consistent anatomical locations in case %i\n',CaseID);
        end
        %% check the correctness of the centerline:
        if(0)
            H = figure('color',[1 1 1]);
            [im,hd] = io_ReadMedicalImage(MaskFile);
            show_segment_surface(im,hd.Mv2w);
            hold on
            Centerline.PlotSpline();
            for iP = 1:numel(iAnatomical)
                Centerline.PlotSplinePoint(iAnatomical(iP));
            end
            view(90,90)
        end
    end
end

if(bSaveHR)   
    save(fullfile(DestinyRootDir,'HeartRate.mat'),'HR');
end