function [ PWV_TTF, PWV_TTP, PWV_XCor ] = ...
    PulseWaveVelocity( ls, dist, t0, tE, time_step, saveDir, newIni, newEnd )

%%
bRunMain = 1;
if exist(fullfile(saveDir,'PWV.mat'),'file')
    load(fullfile(saveDir,'PWV.mat'))
    bRunMain = 0;
end

if bRunMain
    %% Pulse Wave Velocity
    
    Q = transpose(ls)/1000; % from mm3/s to ml/s
    
    figure
    subplot(1,2,1)
    surf( dist, ((t0:tE)-1)*time_step*1000, Q );
    
    axis tight
    axis square
    colormap( jet(10000) )
    xlabel('Centerline length (mm)')
    ylabel('Time (ms)')
    title( 'Flow rate (ml/s)' )
    drawnow
    
    %% Smoothing spline surface fit to original data
    t = t0 : tE;
    s = dist;
    z = double(Q);
    TOL = 100;
    sp = spaps( {t,s}, z, { TOL, TOL/length(t)*length(s)*10 } );
    tstep = 0.005;
    tishift = 4;
    ti = t(1)-tishift : tstep : t(end);
    %si = [ s(1) : (s(end)-s(1))/300 : s(end) ]';
    si = s;
    sp = fnxtr(sp,3);
    Qi = fnval( sp, {ti, si} );
    
    subplot(1,2,2)
    surf( si, ti*time_step*1000, Qi,'LineStyle','none')
    axis tight
    axis square
    colormap( jet(10000) )
    xlabel('Centerline length (mm)')
    ylabel('Time (ms)')
    title( 'Flow rate (ml/s)' )
    drawnow
    
    figure
    contour( si, ti*time_step*1000, Qi,[0 0],'LineWidth',4,'LineColor',[0 0 1])
    axis tight
    xlabel('Centerline length (mm)')
    ylabel('Time (ms)')
    title( 'Flow rate (ml/s) = 0' )
    drawnow
    %%
    %  Estimation of PWV with Time-To-Foot, Time-To-Peak, and XCor as in
    %  "Estimation of Global Aortic Pulse Wave Velocity by Flow-Sensitive 4D MRI"
    %  by Markl et al., Magnetic Resonance in Medicine 63:1575–1582 (2010)
    for j = 1:size(Qi,2);
        spi = spaps( ti, Qi(:,j), 0 );
        spid = fnder(spi);
        TTF(j) = ( -fnval( spi, t(1) ) / fnval( spid, t(1) ) + t(1) )*time_step*1000; % Time-To-Foot
        %{
        %ztmp = ( -fnval( spi, t(1) ) / fnval( spid, t(1) ) + t(1) );
        z = fnzeros(spi);
        z = z(1,:);
        z = z(~(z>tE/2));
        TTF(j) = z(end)*time_step*1000;
        %}
        [~, TTP(j)] = fnmin(fncmb(spi,-1)); % Time-To-Peak
        TTP(j) = TTP(j) * time_step * 1000;
    end
    
    ti0 = round(tishift/tstep) + 1;
    for j = 1:size(Qi,2);
        % XCor (cross correlation)
        [ acor, lag ] = xcorr( Qi(ti0:end,j), Qi(ti0:end,1) );
        [ ~, I ] = max( abs(acor) );
        XCor(j) = lag(I)*tstep * time_step * 1000;
    end
    %%
    newIni_old = newIni;
    newEnd_old = newEnd;
    clear newIni newEnd
    save(fullfile(saveDir,'PWV.mat'));
    newIni = newIni_old;
    newEnd = newEnd_old;
end

%% Plotting PWV estimation results
%openfig( fullfile( saveDir, 'centerline.fig' ) );
%newIni = input('First point on centerline (in mm): ');
%newEnd = input(' Last point on centerline (in mm): ');
load(fullfile(saveDir,'PWV.mat'))
[~,idxI] = min(abs(si - newIni));
[~,idxE] = min(abs(si - newEnd));

siO = si;
TTFO = TTF;
TTPO = TTP;
XCorO = XCor;

si   =   si(idxI:idxE);
TTF  =  TTF(idxI:idxE);
TTP  =  TTP(idxI:idxE);
XCor = XCor(idxI:idxE);

spTTFO = spaps(siO,TTFO,0);
spTTPO = spaps(siO,TTPO,0);
spXCorO = spaps(siO,XCorO,0);
spTTFOder = fnder(spTTFO);
spTTPOder = fnder(spTTPO);
spXCorOder = fnder(spXCorO);
TTFder = fnval(spTTFOder,si);
TTPder = fnval(spTTPOder,si);
XCorder = fnval(spXCorOder,si);

figure

subplot(1,3,1)
f_TTF = fit( si, TTF', 'poly1' );
plot( f_TTF, si, TTF )
hold on
plot( siO, TTFO )
view(-90, 90) %# Swap the axes
set(gca, 'ydir', 'reverse'); %# Reverse the y-axis
legend off
axis equal
axis tight
xlabel('Centerline (mm)')
ylabel('Time (ms)')
PWV_TTF = 1 / f_TTF.p1;
%PWV_TTF = 1/median(TTFder);
%PWV_TTF = 1/mean(TTFder);
title(sprintf('Estimation with TTF\n Pulse Wave Velocity = %0.2f m/s\n',PWV_TTF))

subplot(1,3,2)
f_TTP = fit( si, TTP', 'poly1' );
plot( f_TTP, si, TTP )
hold on
plot( siO, TTPO )
view(-90, 90) %# Swap the axes
set(gca, 'ydir', 'reverse'); %# Reverse the y-axis
legend off
axis equal
axis tight
xlabel('Centerline (mm)')
ylabel('Time (ms)')
PWV_TTP = 1 / f_TTP.p1;
%PWV_TTP = 1/median(TTPder);
%PWV_TTP = 1/mean(TTPder);
title(sprintf('Estimation with TTP\n Pulse Wave Velocity = %0.2f m/s\n',PWV_TTP))

subplot(1,3,3)
f_XCor = fit( si, XCor', 'poly1' );
plot( f_XCor, si, XCor )
hold on
plot( siO, XCorO )
view(-90, 90) %# Swap the axes
set(gca, 'ydir', 'reverse'); %# Reverse the y-axis
legend off
axis equal
axis tight
xlabel('Centerline (mm)')
ylabel('Time (ms)')
PWV_XCor = 1 / f_XCor.p1;
%PWV_XCor = 1/median(XCorder);
%PWV_XCor = 1/mean(XCorder);
title(sprintf('Estimation with XCor\n Pulse Wave Velocity = %0.2f m/s\n',PWV_XCor))