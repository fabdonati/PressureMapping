function im = define_compliance_mask( im ) 
tol = 1;
pp = [];

if isfield(im,'bWall') && isfield(im,'P') && isfield(im,'W')
  B = im.bWall; W = im.W; Np = length(im.P);
  Nx = size(B,1); Ny = size(B,2); Nz = size(B,3);
  N = zeros(size(W)); C = zeros(size(B));
  for i = 1:Nx
    for j = 1:Ny
      for k = 1:Nz
        if B(i,j,k) == 1
          for p = 1:Np
            point = im.P(p).point;
            slope = im.P(p).slope;
            if abs(dot(normaliseV(W(:,i,j,k)' - point), slope)) <= tol
              tol = abs(dot(normaliseV(W(:,i,j,k)' - point), slope));
              N(:,i,j,k) = normaliseV(W(:,i,j,k)' - point);
              pp = p; xpoint = point;
            end
          end
          if sqrt(sum(W(:,i,j,k) + 2* N(:,i,j,k) - xpoint').^2) < sqrt(sum(W(:,i,j,k) - xpoint').^2)
            N(:,i,j,k) = - N(:,i,j,k);
          end
          C(i,j,k) = (pp - 1) / (Np - 1);
          tol = 1;
          pp = [];
        end
      end
    end
  end
end

im.NWall = N;
im.C = C;