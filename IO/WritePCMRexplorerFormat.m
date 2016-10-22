function [] = WritePCMRexplorerFormat(OutDirectory,Mag,Flow,hd)
% Function to write PCMRI data compatible with the PCMR explorer software


if ~exist(OutDirectory,'dir')
    mkdir(OutDirectory);
end

VelocityScaleFactor = 1;
if isfield(hd,'VelocityUnits')
    switch hd.VelocityUnits
        case {'m/s'}
            VelocityScaleFactor = 100;
        case {'dm/s'}
            VelocityScaleFactor = 10;
        case {'cm/s'}
            VelocityScaleFactor = 1;
        case {'mm/s'}
            VelocityScaleFactor = 0.1;
        otherwise
            fprintf('Velocity Scale Factor not recognised in image header: %s\n',hd.VelocityUnits);
    end
end
% The software understands velocity in cm/s, and Ensight data as used in
% the FEM-PPE workflow works in m/s. Therefore, the correction factor will
% here be 100
Flow = VelocityScaleFactor * Flow;

[nx,ny,nz,nT] = size(Mag);
% PCMR explorer seems to like the big endian, as in linux:
endian  = 'b';
type = 'float';


% for iT = 1:nT
%     sT = sprintf('%i',iT-1);
%     write_image_vtk2(fullfile(OutDirectory,['mag_' sT '.vtk']),Mag(:,:,:,iT),hd,[],endian,type);
%     for iC = 1:3
%         sC = sprintf('%i',iC-1);
%         write_image_vtk2(fullfile(OutDirectory,['vct_' sT '_' sC '.vtk']),Flow(:,:,:,iT,iC),hd,[],endian,type);
%     end
% end

emptyimage = zeros(nx,ny,nz);
FlowMag = sqrt(sum(Flow,5));
VelMagAverage = mean(FlowMag,4);
write_image_vtk2(fullfile(OutDirectory,'premask.vtk'),emptyimage,hd,[],endian,type);


%WritePCMRheaders(OutDirectory,hd);
headerfile = fullfile(OutDirectory,'header.sth');
fid = fopen(headerfile,'w');

if fid==-1
    fprintf('WARNING! Not possible to open %s\n',headerfile)
else
    if isfield(hd,'Venc'), 
        Venc = hd.Venc; 
    else
        Venc = 150;
    end

    if isfield(hd,'dT'),
        dT = hd.dT;
    else
        dT = 40;
    end

    fprintf(fid,'time_step_count %i\n',nT);
    fprintf(fid,'venc %i\n',Venc);
    if isfield(hd,'PatientPosition'), 
        PatPos = hd.PatientPosition;
    else
        PatPos = 'unknown';
    end
    fprintf(fid,'patient_orientation %s\n',PatPos);
    fprintf(fid,'image_orientation ASL\n');
    fprintf(fid,'vector_mag_min %1.2f\n',min(Mag(:)));
    fprintf(fid,'vector_mag_max %1.2f\n',max(Mag(:)));
    fprintf(fid,'length_time_step %1.2f\n',dT);
    fprintf(fid,'extended_fields {\n');
    if isfield(hd,'PatientPosition'), 
        InsName = hd.InstitutionName;
    else
        InsName = 'unknown';
    end
    fprintf(fid,'{"Institution Name" "%s"}\n',InsName);
    %fprintf(fid,'{"Institutional Department Name" "OCMR2 "}\n',);
    %fprintf(fid,'{"Requesting Physician" ""}\n',nT);
    %fprintf(fid,'{"Manufacturer" "SIEMENS "}\n',nT);
    %fprintf(fid,'{"Manufacturer's Model Name" "TrioTim "}\n',nT);
    %fprintf(fid,'{"Station Name" "MRC35031"}\n',nT);
    %fprintf(fid,'{"Study Description" "Clinical protocols^OXBAV"}\n',nT);
    if isfield(hd,'SeriesDescription'), 
        SerDes = hd.SeriesDescription;
    else
        SerDes = 'unknown';
    end
    fprintf(fid,'{"Series Description" "%s"}\n',SerDes);
    if isfield(hd,'StudyDate'), 
        StuDat = hd.StudyDate;
    else
        StuDat = 'unknown';
    end
    fprintf(fid,'{"Study Date" "%s"}\n',StuDat);
    if isfield(hd,'StudyTime'), 
        StuTim = hd.StudyTime;
    else
        StuTim = 'unknown';
    end
    fprintf(fid,'{"Study Time" "%s"}\n',StuTim);
    if isfield(hd,'ProtocolName'), 
        ProNam = hd.ProtocolName;
    else
        ProNam = 'unknown';
    end
    fprintf(fid,'{"Protocol Name" "%s"}\n',ProNam);
    %fprintf(fid,'{"Patient's Name" "OXBAV238^AS "}\n',nT);
    if isfield(hd,'PatientID'), 
        PatID = hd.PatientID;
    else
        PatID = 'unknown';
    end
    fprintf(fid,'{"Patient ID" "%s"}\n',PatID);
    %fprintf(fid,'{"Patient's Birth Date" "19961104"}\n',nT);
    if isfield(hd,'PatientSex'), 
        PatSex = hd.PatientSex;
    else
        PatSex = 'unknown';
    end
    fprintf(fid,'{"Patients Sex" "%s"}\n', PatSex);
    if isfield(hd,'PatientWeight'), 
        PatWei = hd.PatientWeight;
    else
        PatWei = 'unknown';
    end
    fprintf(fid,'{"Patients Weight" "%i"}\n',PatWei);
    %fprintf(fid,'{"Medical Alerts" ""}\n',nT);
    fprintf(fid,'}\n');

    fclose(fid);
end

