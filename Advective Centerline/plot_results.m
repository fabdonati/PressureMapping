function plot_results( werp, p1, p2, scheme)

nFrames = size(werp.pdo,2);

%%
figure

ax1 = subplot(121);
plot(werp.pdo  * 1e3/133.33,'k','linewidth',2),hold on
plot(werp.kpdo * 1e3/133.33,'b','linewidth',2),hold on
plot(werp.apdo * 1e3/133.33,'g','linewidth',2),hold on
plot(werp.vpdo * 1e3/133.33,'r','linewidth',2),hold on
axis tight
legend('WERP PD','KINETIC','ADVECTIVE','VISCOUS','Location','SouthEast')
s = sprintf(' OUTFLUX RESULTS BETWEEN POINT %i AND %i WITH %s',p1,p2, scheme);
title(s);
%title(['OUTFLUX RESULTS' num2str(p1)-num2str(p1)])

ax2 = subplot(122);
plot(werp.pdi  * 1e3/133.33,'k','linewidth',2),hold on
plot(werp.kpdi * 1e3/133.33,'b','linewidth',2),hold on
plot(werp.apdi * 1e3/133.33,'g','linewidth',2),hold on
plot(werp.vpdi * 1e3/133.33,'r','linewidth',2),hold on
axis tight
legend('WERP PD','KINETIC','ADVECTIVE','VISCOUS','Location','SouthEast')
s = sprintf(' INFLUX RESULTS BETWEEN POINT %i AND %i WITH %s',p1,p2,scheme);
title(s)
%title('INFLUX RESULTS')

ax1.YLim = [ min( [ ax1.YLim(1); ax2.YLim(1)],[],1) max( [ ax1.YLim(2); ax2.YLim(2)],[],1) ];
ax2.YLim = [ min( [ ax1.YLim(1); ax2.YLim(1)],[],1) max( [ ax1.YLim(2); ax2.YLim(2)],[],1) ];

%%
figure

subplot(321), plot(werp.lambdai/1e9), hold on,plot(linspace(0,length(werp.pdo),10),0*linspace(0,length(werp.pdo),10),'k'); title('INFLUX');
xlim([1 nFrames])
subplot(322), plot(werp.lambdao/1e9), hold on,plot(linspace(0,length(werp.pdo),10),0*linspace(0,length(werp.pdo),10),'k'); title('OUTFLUX');
xlim([1 nFrames])
subplot(323), plot(werp.a/1e6); title('ADVECTIVE ENERGY RATE');
xlim([1 nFrames])
subplot(324), plot(werp.k/1e6); title('KINETIC ENERGY RATE');
xlim([1 nFrames])
subplot(325), plot(werp.v/1e6); title('VISCOUS DISSIPATION RATE');
xlim([1 nFrames])
subplot(326), plot(werp.pdi*1000); hold on, plot(werp.pdo*1000,'r'); legend('INFLUX','OUTFLUX'), title('PRESSURE DROP'); 
xlim([1 nFrames])
s = sprintf( ' RESULTS BETWEEN POINT %i AND %i W %s', p1, p2, scheme );
title(s);

%%
%{
for i = 1:size(im.P(1).gx,1)
  for j = 1:size(im.P(1).gx,2)
    quiver3(im.P(1).gx(i,j),im.P(1).gy(i,j),im.P(1).gz(i,j),...
      im.P(1).N(1,i,j)*im.P(1).b(i,j),im.P(1).N(2,i,j)*im.P(1).b(i,j),im.P(1).N(3,i,j)*im.P(1).b(i,j),10,'k'); hold on,
  end
end
for i = 1:size(im.P(end).gx,1)
  for j = 1:size(im.P(end).gx,2)
    quiver3(im.P(end).gx(i,j),im.P(end).gy(i,j),im.P(end).gz(i,j),...
      im.P(end).V(1,i,j)*im.P(end).b(i,j),im.P(end).V(2,i,j)*im.P(end).b(i,j),im.P(end).V(3,i,j)*im.P(end).b(i,j),1,'k'); hold on,
  end
end
%}