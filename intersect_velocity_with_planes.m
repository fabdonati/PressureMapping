function im = intersect_velocity_with_planes( im, Vtype, hd, opts, frame )

if isfield(im,'P') && isfield(im,Vtype)
    V = eval(['im.' Vtype]);
    V = squeeze( V(frame,:,:,:,:) );
    Nc = size( V, 1 );
    for iP = 1 : length( im.P )
        point = im.P(iP).point;
        slope = im.P(iP).slope;
        
        for c = 1 : Nc
            vPnrrd = Build_nrrd( squeeze( V(c,:,:,:) ), hd );
            vtmp = scinrrd_intersect_plane( vPnrrd, point, slope, opts.interp2d );
            vtmp( isnan(vtmp) ) = 0;
            im.P(iP).V(frame,c,:,:) = im.P(iP).b .* vtmp;
        end
        
        for i = 1 : size( im.P(iP).V, 3 )
            for j = 1 : size( im.P(iP).V, 4 )
                im.P(iP).N(:,i,j) = slope;
                im.P(iP).Vp(i,j)  = dot( squeeze( im.P(iP).V(frame,:,i,j) ), slope );
            end
        end
        
    end
else
    disp('Error! No plane or velocity image assigned!')
    return
end