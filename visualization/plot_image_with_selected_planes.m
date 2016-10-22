function [] = plot_image_with_selected_planes( im, hd ) 

if isfield(im,'b') && isfield(im,'P')
  if isfield(im.P,'point') && isfield(im.P,'slope'), Np = length(im.P);
    figure, show_segment_surface(im.b,hd);
    for i = 1:Np, hold on
      P = im.P(i);
      quiver3(P.point(1),P.point(2),P.point(3),P.slope(1),P.slope(2),P.slope(3),10,'Linewidth',2), hold on,
      text(P.point(1),P.point(2),P.point(3),['\bullet ' num2str(i)],'FontSize',14)
    end
  else
    disp('Error! No planes assigned!')
    return
  end
else
  disp('here')
  disp('Error! No binary image assigned!')
  return
end
