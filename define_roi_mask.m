function im = define_roi_mask( im, hd, id1, id2 )

%%
if nargin < 4
    id2 = length( im.P );
end
if nargin < 3
    id1 = 1;
end

%%
if isfield(im,'W') && isfield(im,'b')
    W = im.W;
    B = im.b;
    Nx = size( W, 2 );
    Ny = size( W, 3 );
    Nz = size( W, 4 );
    if isfield(hd,'Mw2v')
        Mw2v = hd.Mw2v;
    else
        Mw2v = GetMw2v(hd.Mv2w);
    end
    if isfield(im,'P')
        pi = im.P(id1);
        po = im.P(id2);
        pp.P = ( pi.P + po.P )/2;
        for i = 1 : Nx
            for j = 1 : Ny
                for k = 1 : Nz
                    btmp(i,j,k) = ... 
                        dot( pi.N', W(:,i,j,k) - pi.P' ) >= 0 ... 
                                         & ...
                        dot( po.N', W(:,i,j,k) - po.P' ) <= 0 ... 
                                         & ...
                                      B(i,j,k);
                end
            end
        end
        CC = bwconncomp(btmp);
        vox  = round( Mw2v * [ pp.P, 1]');
        voxi = round( Mw2v * [ pi.P, 1]');
        voxo = round( Mw2v * [ po.P, 1]');
        for i = 1 : length( CC.PixelIdxList )
            if length( cell2mat( CC.PixelIdxList(i) ) ) > 1 && ...
                    isempty(find(cell2mat(CC.PixelIdxList(i))==sub2ind(size(btmp),voxi(1),voxi(2),voxi(3)))) && ...
                    isempty(find(cell2mat(CC.PixelIdxList(i))==sub2ind(size(btmp),voxo(1),voxo(2),voxo(3)))) && ...
                    isempty(find(cell2mat(CC.PixelIdxList(i))==sub2ind(size(btmp),vox(1),vox(2),vox(3))))
                btmp(cell2mat(CC.PixelIdxList(i))) = 0;
            else btmp(cell2mat(CC.PixelIdxList(i))) = 1;
            end
        end
        im.bROI = btmp;
    else
        disp('Error! Inlet/outlet planes are not defined!')
        return;
    end
else
    disp('Error! World coordinates are not defined!')
    return;
end