function val  = generate_hexa_scalar ( im )

Nx = size(im,1); Ny = size(im,2); Nz = size(im,3);
c = 1;

if rem(Nx,2) == 0 
  Nx = Nx - 1;
end

if rem(Ny,2) == 0 
  Ny = Ny - 1;
end

if rem(Nz,2) == 0
  Nz = Nz - 1;
end

val = zeros(Nx * Ny * Nz, 1);

for k = 1:Nz
  for j = 1:Ny 
    for i = 1:Nx
      val(c) = im(i,j,k);
      c = c+1;
    end
  end
end