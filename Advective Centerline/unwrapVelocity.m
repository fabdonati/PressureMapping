function VelocityU = unwrapVelocity(Velocity)

Velocity( isnan( Velocity ) ) = 0;
maxVenc = max( max(Velocity(:)), abs(min(Velocity(:)) ) );

VelocityU(:,1,:,:,:) = unwrap( Velocity(:,1,:,:,:)/maxVenc*pi )/pi*maxVenc;
VelocityU(:,2,:,:,:) = unwrap( Velocity(:,2,:,:,:)/maxVenc*pi )/pi*maxVenc;
VelocityU(:,3,:,:,:) = unwrap( Velocity(:,3,:,:,:)/maxVenc*pi )/pi*maxVenc;