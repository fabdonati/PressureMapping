function n  = generate_hexa_nodes( im )

if size(im,1) > 1
  Nx = size(im,2); Ny = size(im,3); Nz = size(im,4);
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
  
  n = zeros(Nx * Ny * Nz, 3);
  
  for k = 1:Nz
    for j = 1:Ny
      for i = 1:Nx
        n(c,:) = im(:,i,j,k);
        c = c+1;
      end
    end
  end
else
  Disp('Error! Missing coordinate field!');
end