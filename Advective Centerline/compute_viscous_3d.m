function [ve, veIm] = compute_viscous_3d( V, B, dx, stencil, mu)

Nx = size(V,2); Ny = size(V,3); Nz = size(V,4);

% Pixel volume (isotropic) ...
dV = dx(1)*dx(2)*dx(3); % Voxel volume

% zeroing returned image and value ...
ve = 0; veIm = zeros(Nx,Ny,Nz);

% Looping over the image and computing viscous components ...
for k = 2:Nz-1
    if(max(max(B(:,:,k))) < 0.5) % checking for any contributions to the sum
        continue
    end
    
    for j = 2:Ny-1
        if(max(B(:,j,k)) < 0.5) % checking for any contributions to the sum
            continue
        end
        
        for i = 2:Nx-1
            % Computing the viscous energy components ...
            veIm(i,j,k) = mu * dV * B(i,j,k) * norm(gradV(V,i,j,k,dx,stencil),'fro')^2;
            ve = ve + veIm(i,j,k);
        end
    end
end