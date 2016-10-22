clear all

%% ADD ALL FOLDERS to path including cardiacmeshing
folder = '/nas/adv/PRESSURE_ESTIMATION/HLHSII/William_Odowd/4dflow/';
name = 'kt8_pca_sparse_';
fld = input('Folder name (initials + /):  ','s');
mkdir(fld);

% COMPUTE SEPARATE MATRIX FOR SYS AND DIAS
ES = input('End systolic frame number:  ');

% Read par-rec and extract velocities in cm/s
for i = 1:3
    parrec2Metafile([folder name num2str(i) '.par'],'prospective');
end

[img1 hdr1] = ReadData3D([folder name '1_2scaled.mhd']);
[img2 hdr2] = ReadData3D([folder name '2_2scaled.mhd']);
[img3 hdr3] = ReadData3D([folder name '3_2scaled.mhd']);

% Get encoding velocity
Venc = max(max(max(max(img1))));

% Build Velocity matrix for pressure estimation in mm/s 
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

% Average velocity map at PEAK SYSTOLE
%for t = 1:size(Velocity,1)
%t=8;
t=round(ES/2)+1;
vs(:,:,:,1)=squeeze(VelocityS(t,1,:,:,:));
vs(:,:,:,2)=squeeze(VelocityS(t,2,:,:,:));
vs(:,:,:,3)=squeeze(VelocityS(t,3,:,:,:));
vmags = sqrt(sum( vs.^2 , 4 ));
write_image_vtk2([fld 'test_sys.vtk'],vmags,struct('origin',hdr3.Offset(1:3),'spacing',hdr1.PixelDimensions(1:3)),'ascii');

% Average velocity map at PEAK DIASTOLE
%for t = 1:size(Velocity,1)
%t=8;
t=round((length(time)-ES)/2);
vd(:,:,:,1)=squeeze(VelocityD(t,1,:,:,:));
vd(:,:,:,2)=squeeze(VelocityD(t,2,:,:,:));
vd(:,:,:,3)=squeeze(VelocityD(t,3,:,:,:));
vmagd = sqrt(sum( vd.^2 , 4 ));
write_image_vtk2([fld 'test_dias.vtk'],vmagd,struct('origin',hdr3.Offset(1:3),'spacing',hdr1.PixelDimensions(1:3)),'ascii');

% Export matrix 
matrixS = [fld 'Velocity_rs.mat'];
matrixD = [fld 'Velocity_rd.mat'];

Velocity = VelocityS;
save(matrixS,'Velocity');

clear Velocity
Velocity = VelocityD;
save(matrixD,'Velocity');

fid = fopen([fld 'time_step.txt'],'w');
fprintf(fid,'%f ',hdr3.PixelDimensions(4));
fclose(fid)
