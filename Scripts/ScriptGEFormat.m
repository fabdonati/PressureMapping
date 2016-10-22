% Script to read velocity values from GE data

% Script to extract GE data

    DicomDir = 'F:\PPE\PAH-15T-4Dflow\PAH-15T-4Dflow\FU_YUN_PAH-320\22439\';
    MagFolder = '434';
    PhaseFolders(1,:) = '531';
    PhaseFolders(2,:) = '532';
    PhaseFolders(3,:) = '533';
    nSlices = 16;
    nT = 20;
    
    OutDirectory = 'F:\PPE\PAH-15T-4Dflow\PCMRexplorer\';


[hd,Mag,Phase] = ReadGEpcmri(DicomDir,MagFolder,PhaseFolders,nSlices,nT);
WritePCMRexplorerFormat(OutDirectory,Mag,Phase,hd);

a=1;
opts.dt = CaseInfo.deltaT/1000;
cd(CaseDir);
[werp] = WerpPressureEstimation(BinaryMask,Velocity,opts);
