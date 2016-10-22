clear all

%% To write the right format of Velocity.mat (5D single) run WTHOUT ADDING CardiacMeshing to the path. 
%% To write the binary image test*.vtk run WITH CardiacMeshing in the path.
%% To run Werp ADD ALL FOLDERS to path including CardiacMeshing

VTK = 0;
MAT = 0;
FLOWVIS = 1;

folder = '/home/adv10/DATA/PRESSURE_ESTIMATION/HLHSII/HT001/4dflow/';
%folder = '/Arch_Project/HLHSII/HT013/4Dflow/';
name = 'kt8_pca_sparse_';
fld = input('Folder name (initials + /):  ','s');
mkdir(fld);

% COMPUTE SEPARATE MATRIX FOR SYS AND DIAS
ES = input('End systolic frame number:  ');

if MAT == 1
    rmpath(genpath('/home/adv10/DATA/PRESSURE_ESTIMATION/cardiacmeshing'));
end

% Read par-rec and extract velocities in cm/s
for i = 1:3
    parrec2Metafile([folder name num2str(i) '.par'],'prospective');
end

[img1 hdr1] = ReadData3D([folder name '1_2scaled.mhd']);
[img2 hdr2] = ReadData3D([folder name '2_2scaled.mhd']);
[img3 hdr3] = ReadData3D([folder name '3_2scaled.mhd']);

% Get encoding velocity
Venc = max(max(max(max(img1))));

%% Build Velocity matrix for pressure estimation in mm/s 
% [time_i  [u_i  v_i  w_i]  [x_i  y_i  z_i]]
time = importdata([folder name '1_timing.txt']);
for i = 1:ES
    VelocityS(i,1,:,:,:) = img1(:,:,:,i)*10;
    VelocityS(i,2,:,:,:) = img2(:,:,:,i)*10;
    VelocityS(i,3,:,:,:) = img3(:,:,:,i)*10;
end
for i = ES+1:length(time)
    VelocityD(i-ES,1,:,:,:) = img1(:,:,:,i)*10;
    VelocityD(i-ES,2,:,:,:) = img2(:,:,:,i)*10;
    VelocityD(i-ES,3,:,:,:) = img3(:,:,:,i)*10;
end

if FLOWVIS == 1
    for i =1:length(time)    
        phase(:,:,:,i,1) =  img1(:,:,:,i);
        phase(:,:,:,i,2) =  img2(:,:,:,i);
        phase(:,:,:,i,3) =  img3(:,:,:,i);
        vmag(:,:,:,i)=sqrt(sum( img1(:,:,:,i).^2, 4));
    end
    mag = [fld 'vmag.mat'];
    flow = [fld 'flow.mat'];
    save(mag, 'vmag');
    save(flow,'phase');
    hdr1.origin=hdr1.Offset(1:3);
    hdr1.spacing=hdr1.PixelDimensions(1:3);
    hdr1.PatientPosition = 'xxx';
    hdr1.InstitutionName = 'sth';
    hdr1.SeriesDescription = 'xxx';
    hdr1.StudyDate = 'xxx';
    hdr1.StudyTime = 'xxx';
    hdr1.ProtocolName = 'xxx';
    hdr1.PatientSex = 'xxx';
    hdr1.PatientWeight = 'xxx';
    hdr1.PatientID= 'xxx';
    mkdir /home/adv10/DATA/PRESSURE_ESTIMATION/INPUT/' fld 'FLOWVIS/
    WritePCMRexplorerFormat(['/home/adv10/DATA/PRESSURE_ESTIMATION/INPUT/' fld 'FLOWVIS/'],vmag,phase,hdr1);
    copyfile(['/home/adv10/DATA/PRESSURE_ESTIMATION/INPUT/' fld 'FLOWVIS/'],['/home/adv10/vm_share/' fld]);
end
% Average velocity map at PEAK SYSTOLE
%for t = 1:size(Velocity,1)
%t=8;
ts=round(ES/2)+1;
vs(:,:,:,1)=squeeze(VelocityS(ts,1,:,:,:));
vs(:,:,:,2)=squeeze(VelocityS(ts,2,:,:,:));
vs(:,:,:,3)=squeeze(VelocityS(ts,3,:,:,:));
vmags = sqrt(sum( vs.^2 , 4 ));
% Average velocity map at PEAK DIASTOLE
%for t = 1:size(Velocity,1)
%t=8;
td=round((length(time)-ES)/2);
vd(:,:,:,1)=squeeze(VelocityD(td,1,:,:,:));
vd(:,:,:,2)=squeeze(VelocityD(td,2,:,:,:));
vd(:,:,:,3)=squeeze(VelocityD(td,3,:,:,:));
vmagd = sqrt(sum( vd.^2 , 4 ));
if VTK == 1
    addpath(genpath('/auto/nas/heart-nas/adv/PRESSURE_ESTIMATION/cardiacmeshing'));
    write_image_vtk2([fld 'test_sys.vtk'],vmags,struct('origin',hdr3.Offset(1:3),'spacing',hdr1.PixelDimensions(1:3)),'ascii');
    write_image_vtk2([fld 'test_dias.vtk'],vmagd,struct('origin',hdr3.Offset(1:3),'spacing',hdr1.PixelDimensions(1:3)),'ascii');
end

% Export matrix 
matrixS = [fld 'Velocity_rs.mat'];
matrixD = [fld 'Velocity_rd.mat'];

Velocity = single(VelocityS);
save(matrixS,'Velocity');

clear Velocity
Velocity = single(VelocityD);
save(matrixD,'Velocity');

fid = fopen([fld 'time_step.txt'],'w');
fprintf(fid,'%f ',hdr3.PixelDimensions(4));
fclose(fid)
