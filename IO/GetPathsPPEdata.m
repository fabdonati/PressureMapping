function [paths] = GetPathsPPEdata()

paths.DirectoryTemplate = 'E:\data\PressureCases\Code and data structures\GenericCase\';

paths.RootName    = 'OpenCMISS';
paths.Path2mask   = '04_EnsightMasking/OUTPUT/';
paths.Path2step7  = '07_EnsightVisualisation/';
paths.Path2output = '08_PostProcessing/';
paths.Path2extra  = '09_AdditionalData/';
paths.Path2PresComps = 'PressureComponents/';
paths.Path2originalEnsight = '02_EnsightRegistration/OUTPUT/';
paths.Path2ConversionSettings = '03_OpenCMISSDomain/';
paths.Step02out = '02_EnsightRegistration/OUTPUT/';
paths.Step03outEnsight = '03_OpenCMISSDomain/OUTPUT/Ensight/';
paths.Step05outOpenCMISS = '05_OpenCMISSAreas/OUTPUT/Ensight/';
paths.Step04out = '04_EnsightMasking/OUTPUT/';
paths.Step06out = '06_OpenCMISSSimulation/output/';
paths.Magnitudes = [paths.Path2output 'Magnitudes/'];
paths.Speed_Sum = [paths.Path2output 'Speed_Sum/'];

%This data is deleted for space reasons:
%paths.Path2cropedEnsight = '03_OpenCMISSDomain/OUTPUT/Ensight/';
paths.Path2cropedEnsight = '04_EnsightMasking/OUTPUT/';

% And also some naming:
%dirOutput = [Directory Path2output];

% Specific files
paths.BinaryImName          = 'aorta2.vtk';
paths.BinaryImWithPath      = '04_EnsightMasking/OUTPUT/aorta2.vtk';
paths.ImageAverageVelocity  = '04_EnsightMasking/OUTPUT/Im_SSVaverage.vtk';
paths.nameSplineFile        = [paths.Path2output 'spline.mat'];
paths.Skeleton              = '08_PostProcessing/aorta-sk-pruned2.mat';
% Some cases might have manual hacks of the original skeleton:
paths.AutomSkeleton         = '08_PostProcessing/aorta-sk-pruned2Original.mat';
% Beats per minute
paths.bpmFile               = [paths.Path2extra 'bpm.txt'];
% Conversion settings:
paths.ConvSet               = [paths.Path2ConversionSettings 'ConversionSettings'];
% Original Case file:
paths.OriginalCaseFile      = [paths.Path2originalEnsight 'OpenCMISS.case'];
% Flag to indicate if conversion settings are locked:
paths.LockConvSetts         = [paths.Path2ConversionSettings 'bLock.txt'];