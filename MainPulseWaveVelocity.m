function [ Q, Qi, PWV_TTF, PWV_TTP, PWV_XCor ] = ... 
    MainPulseWaveVelocity( caseId, nCLpts )

close all

load( [ 'Data\Adelaide\INPUT\' caseId '\Velocity_rs.mat' ] );
im.V = Velocity;
[ im.b, hd ] = io_ReadMedicalImage( [ 'Data\Adelaide\INPUT\' caseId '\sys_seg.vtk' ] );
im.b = vel2bin( im.b, Velocity, 200, 1 );
im.b = round( smooth3( im.b, 'box', 5 ) + 0.2 );
[ Q, Qi, PWV_TTF, PWV_TTP, PWV_XCor ] = PulseWaveVelocity( im, hd, nCLpts, caseId );
