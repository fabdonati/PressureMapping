function plot_results_all( werp, Na, bScaled )

lims.lambdai = [min(min(werp.lambdai/1e9)); max(max(werp.lambdai/1e9))];
lims.lambdao = [min(min(werp.lambdao/1e9)); max(max(werp.lambdao/1e9))];
lims.lambda  = [min(-lims.lambdai(1),lims.lambdao(1)); max(-lims.lambdai(2),lims.lambdao(2))];
lims.a       = [min(min(werp.a/1e6)); max(max(werp.a/1e6))];
lims.k       = [min(min(werp.k/1e6)); max(max(werp.k/1e6))];
lims.v       = [min(min(werp.v/1e6)); max(max(werp.v/1e6))];
lims.pd      = [min(min(min(werp.pdi*1000/133.33)), min(min(werp.pdo*1000/133.33))); ...
  max(max(max(werp.pdi*1000/133.33)), max(max(werp.pdo*1000/133.33)));];

colors = {'k','b','r','g','m','c','y'};
titles = {'SV','AA1','AA2','Arch','DA1','DA2','DA3'};

% % % figure,
% % % for iA = 1:size(werp.lambdai,1)
% % %   subplot(1,7,iA), 
% % %   plot(-werp.lambdai(iA,:) / 1e9, colors{iA}, 'linewidth', 2), hold on
% % %   plot(werp.lambdao(iA,:) / 1e9, [colors{iA} '--*'], 'linewidth', 2), hold on, 
% % %   plot(linspace(0,size(werp.lambdai,2),10), linspace(0,0,10),'k');
% % %   title(titles{iA}); legend('-INFLUX','OUTFLUX')
% % %   ylim([lims.lambda])
% % % end

% % % % INFLUX PRESSURE DROPS
% % % figure
% % % subplot(221)
% % % for iA = 1:size(werp.kpdi,1), plot(werp.kpdi(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('INFLUX KINETIC PRESSURE DROP')
% % % 
% % % subplot(222)
% % % for iA = 1:size(werp.apdi,1), plot(werp.apdi(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('INFLUX ADVECTIVE PRESSURE DROP')
% % % 
% % % subplot(223)
% % % for iA = 1:size(werp.vpdi,1), plot(werp.vpdi(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('INFLUX VISCOUS PRESSURE DROP')
% % % 
% % % subplot(224)
% % % for iA = 1:size(werp.pdi,1), plot(werp.pdi(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('INFLUX PRESSURE DROP')

% % % % OUTFLUX PRESSURE DROPS
% % % figure
% % % subplot(221)
% % % for iA = 1:size(werp.kpdo,1), plot(werp.kpdo(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('OUTFLUX KINETIC PRESSURE DROP')
% % % 
% % % subplot(222)
% % % for iA = 1:size(werp.apdo,1), plot(werp.apdo(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('OUTFLUX ADVECTIVE PRESSURE DROP')
% % % 
% % % subplot(223)
% % % for iA = 1:size(werp.vpdo,1), plot(werp.vpdo(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('OUTFLUX VISCOUS PRESSURE DROP')
% % % 
% % % subplot(224)
% % % for iA = 1:size(werp.pdo,1), plot(werp.pdo(iA,:)  * 1e3/133.33, colors{iA}, 'linewidth', 2), hold on, end
% % % legend('SV','AA1','AA2','Arch','DA1','DA2','DA3'), title('OUTFLUX PRESSURE DROP')

% PLOT DROPS AT EACH LOCATION
if nargin < 3, bScaled = 0; end

             
figure() ;  clf ;
set( gcf, 'Color', [0.8,0.8,0.8], 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] ) ;

for iA = 1:size(werp.pdo,1)
  subplot(5,7,iA), grid on, hold on,
  plot(-werp.lambdai(iA,:)/1e9,'b'), hold on,
  plot(werp.lambdao(iA,:)/1e9,'r'), hold on,
  plot(linspace(0,length(werp.pdo(iA,:)),10),0*linspace(0,length(werp.pdo(iA,:)),10),'k'); hold on
  xlim([0 size(werp.pdo,2)])
  if bScaled, ylim([lims.lambda]), end
  title('FLUXES');
  legend('INFLUX','OUTFLUX')
 
  subplot(5,7,iA+Na), grid on, hold on,
  plot(werp.a(iA,:)/1e6); 
  xlim([0 size(werp.pdo,2)])
  if bScaled, ylim([lims.a]), end
  title('ADVECTIVE');
  
  subplot(5,7,iA+2*Na), grid on, hold on,
  plot(werp.k(iA,:)/1e6); 
  xlim([0 size(werp.pdo,2)])
  if bScaled, ylim([lims.k]), end
  title('KINETIC');
  
  subplot(5,7,iA+3*Na), grid on, hold on,
  plot(werp.v(iA,:)/1e6); 
  xlim([0 size(werp.pdo,2)])
  if bScaled, ylim([lims.v]), end
  title('VISCOUS');
  
  subplot(5,7,iA+4*Na), grid on, hold on,
  plot(werp.pdi(iA,:)*1000/133.33); hold on, 
  plot(werp.pdo(iA,:)*1000/133.33,'r'); 
  xlim([0 size(werp.pdo,2)])
  if bScaled, ylim([lims.pd]), end
  legend('IN','OUT','location','SouthEast'), title('PRESSURE DROP');
end

axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', [0.8,0.8,0.8], 'YColor', [0.8,0.8,0.8] ) ;
text( 0.175, 0, 'SV',  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
text( 0.290, 0, 'AA1', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
text( 0.405, 0, 'AA2', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
text( 0.520, 0, 'Arch','FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
text( 0.630, 0, 'DA1', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
text( 0.750, 0, 'DA2', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
text( 0.865, 0, 'DA3', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;