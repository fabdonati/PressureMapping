function [] = save_cheart_velocity( V, opts, t )
  
if isfield(opts,'MeshName')
  name = opts.MeshName;
  Vel = generate_hexa_nodes( V );
  fid = fopen([name '_VEL-' num2str(t) '.D'],'w+');
  fprintf(fid,'%i %i\n',size(Vel,1),size(Vel,2));
  for i = 1:size(Vel,1)
    fprintf(fid,'%d %d %d\n',Vel(i,:));
  end
  fclose(fid);
end

