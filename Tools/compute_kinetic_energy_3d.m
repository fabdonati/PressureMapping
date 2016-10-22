function [ke keV] = compute_kinetic_energy_3d(V, B, dx, stencil, rho)


Nx = size(V,2); Ny = size(V,3); Nz = size(V,4);

% Pixel volume (isotropic) ...
dV = dx(1)*dx(2)*dx(3); % Voxel volume

% zeroing returned image and value ...
ke = 0; keV = zeros(Nx,Ny,Nz);

% Looping over the image and computing kinetic components ...
for k = 2:Nz-1
  if(max(max(B(:,:,k))) < 0.5) % checking for any contributions to the sum
    continue
  end
  
  for j = 2:Ny-1
    if(max(B(:,j,k)) < 0.5) % checking for any contributions to the sum
      continue
    end
    
    for i = 2:Nx-1
      % Computing the kinetic energy components ...
      keV(i,j,k) = 0.5 * rho * dV * B(i,j,k) * dot(mean3d(V,i,j,k,stencil),mean3d(V,i,j,k,stencil));
      
      % Adding the kinetic energy components to the sum ...
      ke = ke + keV(i,j,k);
    end
  end
end
