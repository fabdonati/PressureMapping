function im = intersect_binarymask_with_planes( im, hd, opts )

if isfield(im,'P')
    b = im.b;
    bPnrrd = Build_nrrd( b, hd );
    for iP = 1 : length( im.P )
        point = im.P(iP).point;
        slope = im.P(iP).slope;
        [ btmp, gx, gy, gz, midx ] = scinrrd_intersect_plane( bPnrrd, point, slope, opts.interp2d );
        %   [ ~, gx, gy, gz ]       = scinrrd_intersect_plane( bPnrrd, point, slope, 'linear' );
        im.P(iP).gx = gx;
        im.P(iP).gy = gy;
        im.P(iP).gz = gz;
        im.P(iP).dx = evaluate_pixeldim( im.P(iP).gx, im.P(iP).gy, im.P(iP).gz );
        % Detect the "blob" that is connected to voxel midx in btmp:
        btmp( isnan(btmp) ) = 0;
        CC = bwconncomp( btmp );
        for i = 1 : length( CC.PixelIdxList )
            if ~isempty( find( cell2mat( CC.PixelIdxList(i) ) == sub2ind( size(btmp), midx(1), midx(2) ) ) )
                btmp(cell2mat(CC.PixelIdxList(i))) = 1;
            else btmp(cell2mat(CC.PixelIdxList(i))) = 0;
            end
        end
        %figure
        %imshow(btmp);
        %drawnow;
        im.P(iP).b = btmp;
        
        for i = 1 : size( im.P(iP).b, 2 )
            for j = 1 : size( im.P(iP).b, 3 )
                im.P(iP).N(:,i,j) = slope;
            end
        end
        
    end
else
    disp('Error! No plane or velocity image assigned!')
    return
end