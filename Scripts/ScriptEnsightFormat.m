% Script to compute pressure differences in a vessel


CaseDir = 'F:\PPE\BAV\OXBAV001\pressure_mapping\04_EnsightMasking\OUTPUT\';
BinaryMask = 'aorta2.vtk';

[Imageheader,points] = ReadGeoFile(CaseDir,'OpenCMISS.geo',1);
CaseInfo = ReadCaseFile(CaseDir,'OpenCMISS');
Velocity = NaN*ones([CaseInfo.nSteps,3,Imageheader.dim(1:3)]);
for iFrame = 1:CaseInfo.nSteps
    fileName = sprintf('OpenCMISS.Velocity%04i',iFrame-1);
    [values,ReadResult,Magnitude] = ReadEnsightValues(CaseDir,fileName,Imageheader,3);
    for ic = 1:3
        Velocity(iFrame,ic,:,:,:) = values(:,:,:,ic);
    end
end

opts.dt = CaseInfo.deltaT/1000;
cd(CaseDir);
[werp] = WerpPressureEstimation(BinaryMask,Velocity,opts);
