% Script to test the WritePCMRexplorerFormat

DataDir = 'F:\PPE\BAV\OXBAV006\pressure_mapping\02_EnsightRegistration\OUTPUT\';
OutDir = 'F:\euHeart\code\PressureMapping\TestPCMR';
[EnsightData] = io_ReadEnsightDirectory(DataDir);
WritePCMRexplorerFormat(OutDir,EnsightData.Magnitude,EnsightData.Velocity,EnsightData.Imageheader);