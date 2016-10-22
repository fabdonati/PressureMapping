function val = colvec(vec)
%  Function - ApplyGaussianFilter
%
%  : This function ensures a vector is oriented by column
% 
%  by D Nordsletten
%  2014
%
  val = vec;
  if(size(vec,1) < size(vec,2))
     val = vec'; 
  end