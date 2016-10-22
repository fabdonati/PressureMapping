function [AverageSquaredValue,Imageheader] = CalculateAverageSumSquared(directoryIN,RootName,CaseInfo,NameConvention,directoryOUT,NameMagnitudeIN,NameOut)
% Function to calculate the average image of the Sum of Squared values of
% velocity in a phase contrast MRI sequence, using the sequence of frames
% in ensight format:

if nargin<4
    NameConvention = 2;
end
if nargin<5
    directoryOUT = directoryIN;
end
if nargin<6
    NameMagnitudeIN = 'Speed_SumSquares';
end
if nargin<7
    NameOut = 'Im_SSVaverage';
end
GEOfileName = [RootName '.geo'];
bReadPoints=1;
nFrames = CaseInfo.nSteps;
T0 = CaseInfo.T0filename;
[Imageheader] = ReadGeoFile(directoryIN,GEOfileName,bReadPoints);
ImSize = Imageheader.dim;

values = zeros(ImSize(1),ImSize(2),ImSize(3),nFrames);
for iFrame=1:nFrames
    iFr2read = iFrame-1+T0;
    switch NameConvention
        case 1            
            fileNameIn   = [RootName '.' NameMagnitudeIN '_' sprintf('%02i',iFr2read)];
        case 2
            fileNameIn   = [RootName '.' NameMagnitudeIN sprintf('%04i',iFr2read)];
    end
    values(:,:,:,iFrame) = ReadEnsightValues(directoryIN,fileNameIn,Imageheader);
end

AverageSquaredValue = mean(values,4);

imageOut = [directoryOUT NameOut '.vtk'];
write_image_vtk2(imageOut,AverageSquaredValue,Imageheader,'ascii');