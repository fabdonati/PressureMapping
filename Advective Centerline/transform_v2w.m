function [ im, hd ] = transform_v2w( im, hd )
% 'im'      = binary image
% 'hd'      = header of image
% 'spacing' = voxel dim
% 'or'      = origin
%%
if isfield(im,'b')
    Nx = size( im.b, 1 );
    Ny = size( im.b, 2 );
    Nz = size( im.b, 3 );
    im.W  = zeros( 3, Nx, Ny, Nz);
    if isfield( hd, 'spacing' )   &&   isfield( hd, 'origin' )
        hd.Mv2w = [ hd.spacing(1)         0               0         hd.origin(1); ...
            0         hd.spacing(2)         0         hd.origin(2); ...
            0               0         hd.spacing(3)   hd.origin(3) ];
        for i = 1:Nx
            for j = 1:Ny
                im.W(:,i,j,:) = hd.Mv2w*[ (i-1)*ones(1,Nz); (j-1)*ones(1,Nz); (1:Nz)-1; ones(1,Nz) ];
            end
        end
        
    else
        disp('Error! No voxeldim or origin assigned!')
        return
    end
else
    disp('Error! No binary image assigned!')
    return
end