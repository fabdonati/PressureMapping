clc
if ( ~exist('im') || ~exist('hd') )
    clear
    load('TestData\Cases1\20051\Data.mat');
end
[ Q, PWV_TTF, PWV_TTP, PWV_XCor ] = PulseWaveVelocity( im, hd );