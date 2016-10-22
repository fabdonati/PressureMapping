function gx = d_dx(Im,i,j,k,dx,stencil)

switch stencil
  case 'standard'
    gx  = (Im(:,i+1,j,k)   - Im(:,i-1,j,k)) / (2*dx);
    
  case 'filtered'
    gx  = (Im(:,i+1,j,k)   - Im(:,i-1,j,k)) / (2*dx);
    
    % Estimate by expansion in cross-component ...
    gpy = (Im(:,i+1,j+1,k) - Im(:,i-1,j+1,k)) / (2*dx);
    gmy = (Im(:,i+1,j-1,k) - Im(:,i-1,j-1,k)) / (2*dx);
    
    % Estimate by expansion in cross-component ...
    gpz = (Im(:,i+1,j,k+1) - Im(:,i-1,j,k+1)) / (2*dx);
    gmz = (Im(:,i+1,j,k-1) - Im(:,i-1,j,k-1)) / (2*dx);
    
    % Total Estimated gradient ...
    gx = (gx + 0.5 * (gpy + gmy) + 0.5 * (gpz + gmz)) / 3;
end
