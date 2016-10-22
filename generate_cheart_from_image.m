function [] = generate_cheart_from_image ( W, B, opts )

if isfield(opts,'MeshName')
  name = opts.MeshName;
  Conn = generate_hexa_connectivity( W );
  fid = fopen([name '_FE.T'],'w+');
  fprintf(fid,'%i %i\n',size(Conn,1),max(max(Conn)));
  for j = 1:size(Conn,1),
    fprintf(fid,'%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n',Conn(j,:));
  end
  fclose(fid);
  
  Nodes = generate_hexa_nodes( W );
  fid = fopen([name '_FE.X'],'w+');
  fprintf(fid,'%i %i\n',size(Nodes,1),size(Nodes,2));
  for j = 1:size(Nodes,1)
    fprintf(fid,'%d %d %d\n',Nodes(j,:));
  end
  fclose(fid);
  
  Mask = generate_hexa_scalar( B );
 
  fid = fopen([name '_cMASK_FE.X'],'w+');
  gid = fopen([name '_vMASK_FE.X'],'w+');
  fprintf(fid,'%i %i\n',size(Conn,1),1);
  fprintf(gid,'%i %i\n',size(Conn,1),1);
  
  for i = 1:size(Conn,1)
    if all(Mask(Conn(i,:))) == 1
      fprintf(fid,'%i\n',1);
    else
      fprintf(fid,'%i\n',0);
    end
  end
  fclose(fid);
  
  for i = 1:size(Conn,1)
    if any(Mask(Conn(i,:))) == 1
      fprintf(gid,'%i\n',1);
    else
      fprintf(gid,'%i\n',0);
    end
  end
  fclose(gid);
  
  fid = fopen([name '_cMASK_FE.T'],'w+');
  gid = fopen([name '_vMASK_FE.T'],'w+');
  fprintf(fid,'%i %i\n',size(Conn,1),size(Conn,1));
  for i = 1:size(Conn,1)
    fprintf(fid,'%i\n',i);
  end
  fclose(fid);
  
  fprintf(gid,'%i %i\n',size(Conn,1),size(Conn,1));
  for i = 1:size(Conn,1)
    fprintf(gid,'%i\n',i);
  end
  fclose(gid);
else
  disp('Error! No output name for cheart mesh files defined!')
  return
end
