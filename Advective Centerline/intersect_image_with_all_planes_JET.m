function Plane = intersect_image_with_all_planes_JET( ...
    im, Vtype, hd, opts, iFrame, first_iA, last_iA, Plane )

Pnrrd = Build_nrrd( im.b, hd );

if isfield( im, Vtype )
    V = eval( ['im.' Vtype] );
    V = permute(V,[3 4 5 2 1]);
    V = V(:,:,:,:,iFrame);
    nC = size( V, 4 ); % Number of components of velocity field
    fprintf( '\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n' );
    for iP = first_iA : last_iA
        ip = iP - first_iA + 1;
        fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
        fprintf( 'Plane %3d of %3d', iP, last_iA );
        Plane(ip).V = [];
        for iC = 1 : nC
            Vel = V( :, :, :, iC );
            Vnrrd = Build_nrrd( Vel, hd );
            % Note: try to use interp3 instead of "scinrrd_intersect_plane"
            tmp = scinrrd_intersect_plane( Vnrrd, Plane(ip).P, Plane(ip).NN, opts.interp2d );
            tmp( isnan(tmp) ) = 0;
            Plane(ip).V(iC,:,:) = tmp;
        end
        %{
        Vx = Plane(ip).V(1,:,:);
        Vy = Plane(ip).V(2,:,:);
        Vz = Plane(ip).V(3,:,:);
        Vp = sqrt( Vx.^2 + Vy.^2 + Vz.^2 );
        Plane(ip).Vp(iFrame,:,:) = squeeze(Vp);
        %}
        btmp = Plane(ip).b;
        %{
        bord = imdilate( btmp, strel(ones(3) ) ) - btmp;
        errM = -ones( size(btmp) );
        Vp0 = Vp( btmp == 1 );
        [ f0, xi0 ] = ksdensity( Vp0, 'support', 'positive' );
        pp0 = csaps( xi0, f0, -1 );
        for i = 1 : 8;
            %             errM = -ones( size(btmp) );
            %             [ f0, xi0 ] = ksdensity( Vp( btmp == 1 ), 'support', 'positive' );
            %             pp0 = csaps( xi0, f0, -1 );
            [ rBord, cBord ] = find( bord == 1);
            for iBord = 1 : length(rBord)
                %Vp0 = Vp( btmp == 1 );
                Vp2 = [ Vp0; Vp( rBord(iBord), cBord(iBord) ) ];
                [ f, xi ] = ksdensity( Vp2, 'support', 'positive' );
                pp = csaps( xi, f, -1 );
                L = max(min(xi0),min(xi));
                U = min(max(xi0),max(xi));
                x = L : (U-L)/99 : U;
                pdf0 = fnval( pp0,x );
                pdf  = fnval( pp, x );
                err = abs( 1 - dot( pdf, pdf0 ) / dot( pdf0, pdf0 ) );
                if  err <= 0.01 && ...
                        dot( Plane(ip).V(:,rBord(iBord), cBord(iBord)), Plane(ip).N ) > 0
                    btmp( rBord(iBord), cBord(iBord) ) = 1;
                end
                errM( rBord(iBord), cBord(iBord) ) = err;
            end
            
            if min(min(errM(errM >=0))) > 0.01
                break
            end
            bord = imdilate( btmp, strel(ones(3) ) ) - btmp;
            bord = bord.*(errM == -1);
        end
        %[ min(min(errM(errM >=0))) i ]
        %}
        Plane(ip).bb = btmp;
    end
    fprintf('\n\n');
else
    disp('Error! No plane or velocity image assigned.')
    return
end