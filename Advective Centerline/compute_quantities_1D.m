function [ as, ao, ai, ls, dist, Plane ] = ...
    compute_quantities_1D( im, Vtype, hd, opts, CenterlineObj, p0, pE, t0, tE, bDebug )

%%
% 'Velocity' is in mm/s
% 'as' is in mmHg
% 'ls' is in l/min

if nargin < 10
    bDebug = 0;
end

if(bDebug)
    figure
    hold on
end

ls = zeros( tE, pE );

%% WARNING: it is not accurate
dist = zeros(1, pE - p0 + 1);
for iP = 2 : pE - p0 + 1;
    dist( iP ) = CenterlineObj.GetDistance( iP ) - CenterlineObj.GetDistance( 1 );
end

%%
nPoints = CenterlineObj.nPointsSpline;
Plane(:).N = zeros( nPoints, 3);
Pnrrd = Build_nrrd( im.b, hd );
for iP = 1 : nPoints;
    Plane(iP).P  =        CenterlineObj.GetPoint( iP );
    Plane(iP).N  = normr( CenterlineObj.GetSlope( iP ) );
    Plane(iP).NN = Plane(iP).N;
    [ btmp, gx, gy, gz, midx ] = ...
        scinrrd_intersect_plane( Pnrrd, Plane(iP).P, Plane(iP).N, opts.interp2d );
    % Detect the component that is connected to pixel midx in btmp:
    btmp( isnan(btmp) ) = 0;
    btmp = round(btmp);
    btmp = double( bwselect( btmp, midx(2), midx(1), 8 ) );
    Plane(iP).b  = btmp;
    Plane(iP).bb = Plane(iP).b;
    Plane(iP).gx = gx;
    Plane(iP).gy = gy;
    Plane(iP).gz = gz;
    Plane(iP).dx = evaluate_pixeldim( gx, gy, gz );
end

ai = zeros(  1, tE );
as = zeros( tE, pE - p0 + 1 );
ao = zeros( tE, pE - p0 + 1 );

for iFrame = t0 : tE
    disp(' ');
    disp( [ 'Time frame ' num2str(iFrame) ' of ' num2str(tE) ] );
    disp(' ');
    Plane = intersect_image_with_all_planes_JET( im, Vtype, hd, opts, iFrame, p0, pE, Plane );
    ai(iFrame) = compute_advective_2d( Plane(1).V, -Plane(1).N, Plane(1).bb, Plane(1).dx, opts.stencil, opts.rho );
    
    for iP = 1 : pE - p0 + 1;
        ls( iFrame, iP ) = compute_lambda_2d(    Plane(iP).V, Plane(iP).N, Plane(iP).bb, Plane(iP).dx );
        ao( iFrame, iP ) = compute_advective_2d( Plane(iP).V, Plane(iP).N, Plane(iP).bb, Plane(iP).dx, opts.stencil, opts.rho );
        as( iFrame, iP ) = -( (ai(iFrame) + ao( iFrame, iP ) ) / ls( iFrame, iP ) )*(1e3/133.33);
        
        bDebug = 1;
        bPlaneVectorField = 0;
        bStreamlines = 0;
        bPathlines = 1;
        if bDebug
            if bPlaneVectorField
                vx = squeeze( Plane(iP).V(1,(Plane(iP).b>0)) );
                vy = squeeze( Plane(iP).V(2,(Plane(iP).b>0)) );
                vz = squeeze( Plane(iP).V(3,(Plane(iP).b>0)) );
                hold on
                quiver3(Plane(iP).gx(Plane(iP).b>0), ...
                    Plane(iP).gy(Plane(iP).b>0), ...
                    Plane(iP).gz(Plane(iP).b>0), vx', vy', vz', 3);
                axis equal
                drawnow
            end
            if bStreamlines || bPathlines
                xvec = unique(squeeze( im.W(1,:,:,:) ));
                yvec = unique(squeeze( im.W(2,:,:,:) ));
                zvec = unique(squeeze( im.W(3,:,:,:) ));
                ix = dsearchn(xvec,Plane(iP).gx(Plane(iP).b>0));
                iy = dsearchn(yvec,Plane(iP).gy(Plane(iP).b>0));
                iz = dsearchn(zvec,Plane(iP).gz(Plane(iP).b>0));
                [X,Y,Z] = meshgrid(xvec,yvec,zvec);
                U = squeeze(im.V(iFrame,1,:,:,:)).*im.b;
                V = squeeze(im.V(iFrame,2,:,:,:)).*im.b;
                W = squeeze(im.V(iFrame,3,:,:,:)).*im.b;
                U = permute(U,[2 1 3]);
                V = permute(V,[2 1 3]);
                W = permute(W,[2 1 3]);
                nFrames = size(im.V,1);
                x1 = unique(squeeze( im.W(1,:,:,:) ));
                x2 = unique(squeeze( im.W(2,:,:,:) ));
                x3 = unique(squeeze( im.W(3,:,:,:) ));
                if bStreamlines
                    for i = 1 : length(ix)
                        [ sx, sy, sz ] = meshgrid(xvec(ix(i)),yvec(iy(i)),zvec(iz(i)));
                        XYZ = stream3(X,Y,Z,U,V,W,sx,sy,sz);
                        hold on
                        streamline(XYZ)
                        drawnow
                    end
                end
                if bPathlines
                    stPoints = [ xvec(ix) yvec(iy) zvec(iz) ];
                    P = velocity2pathlines( im, opts, stPoints, xvec, yvec, zvec );
                    P = permute(P,[2 3 1]);
                    hold on
                    for istP = 1 : size(P,3)
                        plot3(P(1,:,istP),P(2,:,istP),P(3,:,istP));
                        drawnow
                    end
                end
            end
        end
    end
end