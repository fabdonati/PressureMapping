function gz = d_dz(Im,i,j,k,dx,stencil)

switch stencil
  case 'standard'
    gz  = (Im(:,i,j,k+1)   - Im(:,i,j,k-1)) / (2*dx);
    
  case 'filtered'
    gz  = (Im(:,i,j,k+1)   - Im(:,i,j,k-1)) / (2*dx);

    % Estimate by expansion in cross-component ...
    gpx = (Im(:,i+1,j,k+1) - Im(:,i+1,j,k-1)) / (2*dx);
    gmx = (Im(:,i-1,j,k+1) - Im(:,i-1,j,k-1)) / (2*dx);
    
    % Estimate by expansion in cross-component ...
    gpy = (Im(:,i,j+1,k+1) - Im(:,i,j+1,k-1)) / (2*dx);
    gmy = (Im(:,i,j-1,k+1) - Im(:,i,j-1,k-1)) / (2*dx);
    
    % Total Estimated gradient ...
    gz = (gz + 0.5 * (gpx + gmx) + 0.5 * (gpy + gmy)) / 3;
end