function val = AdvectiveBernoulli(vel,rho)

for t = 1:size(vel,1)
    val(t) = -0.5 * rho * (vel(t,end)^2 - vel(t,1)^2);
end