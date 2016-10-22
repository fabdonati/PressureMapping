function val = InertialBernoulli2(dist,vel,rho,dt)

p = 0;
for t = 1:size(vel,1)-1
    for i = 1:size(dist,1)
        p = p + (vel(t+1,i) - vel(t,i))/dt * dist(i)/1000;
    end
    val(t) = - 0.5 * rho * p;
    p = 0;
end