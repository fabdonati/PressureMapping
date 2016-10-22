function[ SubBinaryMask ] = DefineSubBinaryMask( imW, Mask, InletPoint, InletVersor, OutletPoint, OutletVersor)

Nx = size(imW,ndims(imW)-2); Ny = size(imW,ndims(imW)-1); Nz = size(imW,ndims(imW));
SubBinaryMask = zeros(Nx,Ny,Nz);

if ndims(imW) == 4
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
                SubBinaryMask(i,j,k) = dot( InletVersor', imW(:,i,j,k) - InletPoint' )  >= 0 & ...
                    dot( OutletVersor', imW(:,i,j,k) - OutletPoint' ) <= 0 & Mask(i,j,k);
            end
        end
    end
    
elseif ndims(imW) == 3
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
                SubBinaryMask(i,j,k) = dot( InletVersor', imW(i,j,k) - InletPoint' )  >= 0 & ...
                    dot( OutletVersor', imW(i,j,k) - OutletPoint' ) <= 0 & Mask(i,j,k);
            end
        end
    end
end