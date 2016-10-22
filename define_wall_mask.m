function im = define_wall_mask( im, opts, id1, id2 )

%%
if nargin < 4
    id2 = length( im.P );
end
if nargin < 3
    id1 = 1;
end

%%
if isfield(opts,'se'), se = opts.se;
    if isfield(im,'W') && isfield(im,'bROI') && isfield(im,'P')
        
          W = im.W;
        ROI = im.bROI;
         pi = im.P(id1);
         po = im.P(id2);
         
        Eroded        = imerode(ROI,se);
        Boundary      = ROI - Eroded;
        noInletOutlet = ... 
            DefineSubBinaryMask( W, ROI, pi.P + 2*pi.N, pi.N, po.P - 2*po.N, po.N );
        WallMask      = noInletOutlet .* Boundary;
        im.bWall      = WallMask;
    else
        disp('Error! No ROI defined');
        return;
    end
else
    disp('Error! No se defined');
    return;
end