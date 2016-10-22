function [P, p1, p2, tmp] = select_planes( bVTK, im, hd, im2, hd2, opts, ro, co, he )
% Function to indicate the planes between which pressure is computed.
% This can be done by:
% - Automatic default option: extremes of the centerline of the binary mask
% - Manual definition

[ tmp, iExtremePlanes ] = get_planes( bVTK, ro, co, he );
plot_image_with_points( im2.b, hd2, tmp );
show_segment_surface(im2.b,hd2.Mv2w);

bManualPlanes = 0;
if isfield(opts,'bManualPlanes')
    if opts.bManualPlanes == 1
        bManualPlanes = 1;
    end
end
if(~bManualPlanes)
    % Automatic default option:
    pi = iExtremePlanes(1);
    po = iExtremePlanes(2)-1;
    if isfield(opts,'iPoint'), pi = opts.iPoint; end
    if isfield(opts,'oPoint'), po = opts.oPoint; end
    for i = 1:3
        P(i).point = tmp.point(pi+i-1,:);
        P(i).slope = tmp.slope(pi+i-1,:);
    end
    for i = 6:-1:4
        P(i).point = tmp.point(po-6+i,:);
        P(i).slope = tmp.slope(po-6+i,:);
    end
else
    o = input('How many points? ');
    
    p1 = input('Select first point: ');
    for i = 1:o
        P(i).point = tmp.point(p1+(i-1),:);
        P(i).slope = tmp.slope(p1+(i-1),:);
    end
    
    p2 = input('Select last point: ');
    for i = o*2:-1:o+1
        P(i).point = tmp.point(p2+(i-o*2),:);
        P(i).slope = tmp.slope(p2+(i-o*2),:);
    end
end
