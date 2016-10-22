function val = mean3d( Im, i, j, k, stencil )
%  Function - mean3d
%
%  : This function computes the velocity field at (i,j,k)
%  : using stencils to reduce noise
% 
%  Inputs -
%      Im    - the image (3 x Nx x Ny x Nz)
%
%      i,j,k - index in the image to compute derivative
%
%  Outputs -
%
%      val - value of velocity (V)
%
%  by F Donati
%  2014

if nargin < 5
    stencil = 'standard';
end

switch stencil
  
  case 'standard'
    % Straight forward value ...
    val  = Im(:,i,j,k);
    
  case 'filtered'
    % Average of the value with stencil estimates of the value ...
    val1 = 0.5 * Im(:,i,j,k) +  0.25 * (Im(:,i+1,j,k) + Im(:,i-1,j,k));
    val2 = 0.5 * Im(:,i,j,k) +  0.25 * (Im(:,i,j+1,k) + Im(:,i,j-1,k));
    val3 = 0.5 * Im(:,i,j,k) +  0.25 * (Im(:,i,j,k+1) + Im(:,i,j,k-1));
    
    % Final result (the average of the 3 estimates) ...
    val  = (1/3) * ( val1 + val2 + val3 );
end