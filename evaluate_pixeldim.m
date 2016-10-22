function dx = evaluate_pixeldim( gx, gy, gz )

list = find(~isnan(gx) & ~isnan(gy) & ~isnan(gz));
for l = 1:length(list)
  [i, j] = ind2sub(size(gx),list(l));
  if ~isnan(gx(i+1,j)) && ~isnan(gy(i+1,j)) && ~isnan(gz(i+1,j)) && ...
      ~isnan(gx(i,j+1)) && ~isnan(gy(i,j+1)) && ~isnan(gz(i,j+1))
    dx(1) = pdist([ gx(i,j), gy(i,j), gz(i,j) ; gx(i+1,1), gy(i+1,1), gz(i+1,1) ]);
    dx(2) = pdist([ gx(i,j), gy(i,j), gz(i,j) ; gx(i,j+1), gy(i,j+1), gz(i,j+1) ]);
    if dx(1) * dx(2) ~= 0;
      return
    end
  end
end
