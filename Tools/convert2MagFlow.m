function [ Mag, Flow ] = convert2MagFlow( V, B, bB )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B = ceil(B/max(B(:)));
[ nT, nC, nx, ny, nz ] = size(V);
Flow = permute(V,[3 4 5 1 2]);
if(bB)
    for iT = 1 : nT
        for iC = 1 : nC
            Flow(:,:,:,iT,iC) = Flow(:,:,:,iT,iC).*B;
        end
    end
end
Mag = zeros(nx,ny,nz,nT);
for iT = 1 : nT
    Mag(:,:,:,iT) = sqrt( Flow(:,:,:,iT,1).^2 + ...
                          Flow(:,:,:,iT,2).^2 + ...
                          Flow(:,:,:,iT,3).^2 );
end

