function val = mean2d(Im,i,j,stencil)
%  Function - MeanV
%
%  : This function computes the velocity field at (i,j)
%  : using stencils to reduce noise
%
%  Inputs -
%      Im    - the image (3 x Nx x Ny)
%
%      i,j - index in the image to compute derivative
%
%  Outputs -
%
%      val - value of velocity (V)
%
%  by F Donati
%  2014


switch stencil
    
    % Straight forward value
    case 'standard'
        val  = Im(:,i,j);
        
    case 'filtered'
        % Average of the value with stencil estimates of the value
        val1 = 0.5 * Im(:,i,j) +  0.25 * (Im(:,i+1,j) + Im(:,i-1,j));
        val2 = 0.5 * Im(:,i,j) +  0.25 * (Im(:,i,j+1) + Im(:,i,j-1));
        
        % Final result (the average of the 2 estimates)
        val  = (1/2) * ( val1 + val2);
end