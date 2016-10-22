function[ opts ] = get_problem_variables( opts )

if ~isfield(opts,'flow')
  opts.rho = 1060 * 1e-9;
  opts.mu = 0.004 * 1e-3;
end

if ~isfield(opts,'interp2d')
  opts.interp2d = 'nn';
end

if ~isfield(opts,'stencil')
  opts.interp2d = 'filtered';
end

if ~isfield(opts,'AdveFilter')
  opts.AdveFilter = 1;
end

if ~isfield(opts,'compliance')
  opts.compliance = 0;
end

[xx,yy,zz] = ndgrid(-1:1);
if ~isfield(opts,'se')
  opts.se  = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1.0;
end

bVTK            = 'mask.vtk';
opts.MeshName   = 'AORTA';
opts.dt         = 0.0408;
to = 2;
