function plane = intersect_image_with_plane( V, B, hd, P, opts )

for c = 1:3
  Pnrrd = Build_nrrd(squeeze(V(c,:,:,:)),hd);
  [tmp, ~, ~, ~, ~] = scinrrd_intersect_plane(Pnrrd,P.point,P.slope,opts.interp2d);
  Pnrrd = Build_nrrd(B,hd); [btmp, gx, gy, gz, midx] = scinrrd_intersect_plane(Pnrrd,P.point,P.slope,opts.interp2d);
  % Detect the "blob" that is connected to voxel midx in btmp
  btmp(isnan(btmp)) = 0;
  tmp(isnan(tmp)) = 0;
  CC = bwconncomp(btmp);
  for i = 1:length(CC.PixelIdxList)
    if ~isempty(find(cell2mat(CC.PixelIdxList(i))==sub2ind(size(btmp),midx(1),midx(2)),1)), btmp(cell2mat(CC.PixelIdxList(i))) = 1;
    else btmp(cell2mat(CC.PixelIdxList(i))) = 0;
    end
  end
  plane.b = btmp;
  plane.V(c,:,:) = btmp .* tmp;
  plane.gx = gx; plane.gy = gy; plane.gz = gz;
end
for i = 1:size(plane.V,2)
  for j = 1:size(plane.V,3)
    plane.M(i,j)   = sqrt( plane.V(1,i,j)^2 + plane.V(2,i,j)^2 + plane.V(3,i,j)^2);
    plane.Vp(i,j)  = dot(plane.V(:,i,j),P.slope);
    plane.N(:,i,j) = P.slope;
  end
end
plane.dx = evaluate_pixeldim(gx,gy,gz);
plane.bJ = zeros(size(plane.b));

if isfield(opts,'threshold')
  thr = opts.threshold;
else 
  thr = 0.1;
end
for i = 1:size(plane.M,1)
  for j = 1:size(plane.M,2)
    if plane.M(i,j) > thr * max(plane.M(:))
      plane.bJ(i,j) = 1;
    end
  end
end

[plane.i, plane.j] = ind2sub(size(plane.M),find(plane.M(:) == max(plane.M(:)),1));
plane.midx = midx;