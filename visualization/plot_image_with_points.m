function plot_image_with_points( im, hd, point, slope )

figure
show_segment_surface( im, hd );
hold on
for i = 1:size(point,1)
    hold on
    quiver3( point(i,1), point(i,2), point(i,3), slope(i,1), slope(i,2), slope(i,3),20,'b','Linewidth',2)
    hold on
    text( point(i,1), point(i,2), point(i,3), ['\bullet ' num2str(i)],'FontSize',14);
end
