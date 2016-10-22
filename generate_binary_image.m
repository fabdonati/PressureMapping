function im = generate_binary_image( im )

if isfield(im,'I')
  Nx = size(im.I,1); Ny = size(im.I,2); Nz = size(im.I,3);
  im.B = zeros(Nx,Ny,Nz);
  for i = 1:Nx
    for j = 1:Ny
      for k = 1:Nz
        if Im(i,j,k) > 0
          im.B(i,j,k) = 1;
        end
      end
    end
  end
end