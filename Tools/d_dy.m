function gy = d_dy(Im,i,j,k,dx,stencil)

switch stencil
  case 'standard'
    gy  = (Im(:,i,j+1,k)   - Im(:,i,j-1,k)) / (2*dx);
     
  case 'filtered'
    gy  = (Im(:,i,j+1,k)   - Im(:,i,j-1,k)) / (2*dx);
    
    % Estimate by expansion in cross-component ...
    gpx = (Im(:,i+1,j+1,k) - Im(:,i+1,j-1,k)) / (2*dx);
    gmx = (Im(:,i-1,j+1,k) - Im(:,i-1,j-1,k)) / (2*dx);
    
    % Estimate by expansion in cross-component ...
    gpz = (Im(:,i,j+1,k+1) - Im(:,i,j-1,k+1)) / (2*dx);
    gmz = (Im(:,i,j+1,k-1) - Im(:,i,j-1,k-1)) / (2*dx);
    
    % Total Estimated gradient ...
    gy = (gy + 0.5 * (gpx + gmx) + 0.5 * (gpz + gmz)) / 3;
end