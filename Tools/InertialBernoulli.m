function[val] = InertialBernoulli(x,v,rho,dt,scheme,bDebug)

% Np is number of points on centreline
% Nt is number of time frames
% x is the centreline approximated points array with dim : Np x 3
% v is the projected velocity array with dim: Nt, Np

val = zeros(size(v,1),1);
p = 0;

switch scheme
    case 'centered2'
        for t = 2:size(v,1)-1
            for i = 1:size(x,1)-1
                p = p + 0.5 * ( (v(t+1,i) - v(t-1,i))/(2*dt) + (v(t+1,i+1) - v(t-1,i+1))/(2*dt) ) * norm(x(i+1,:) - x(i,:),2);
            end
            val(t) = - rho * p;
            p = 0;
        end
        for i = 1:size(x,1)-1
            val(1) = -rho * ( (v(2,i) - v(end,i))/(2*dt) + (v(2,i+1) - v(end,i+1))/(2*dt) ) * norm(x(i+1,:) - x(i,:),2);
            val(end) = -rho * ( (v(1,i) - v(end-1,i))/(2*dt) + (v(1,i+1) - v(end-1,i+1))/(2*dt) ) * norm(x(i+1,:) - x(i,:),2);
        end
        
        
    case 'forward2'
        for t = 1:size(v,1)-1
            for i = 1:size(x,1)-1
                p = p + 0.5 * ( (v(t+1,i) - v(t,i))/dt + (v(t+1,i+1) - v(t,i+1))/dt ) * norm(x(i+1,:) - x(i,:),2);
            end
            val(t) = - rho * p;
            p = 0;
        end
        for i = 1:size(x,1)-1
            val(end) = val(end) -rho * ( (v(1,i) - v(end,i))/dt + (v(1,i+1) - v(end,i+1))/dt ) * norm(x(i+1,:) - x(i,:),2);
        end
        
    case 'backward2'
        for t = 2:size(v,1)
            for i = 1:size(x,1)-1
                p = p + 0.5 * ( (v(t,i) - v(t-1,i))/dt + (v(t,i+1) - v(t-1,i+1))/dt ) * norm(x(i+1,:) - x(i,:),2);
            end
            val(t) = - rho * p;
            p = 0;
        end
        for i = 1:size(x,1)-1
            val(1) = val(1) -rho * ( (v(1,i) - v(end,i))/dt + (v(1,i+1) - v(end,i+1))/dt ) * norm(x(i+1,:) - x(i,:),2);
        end
end

if bDebug
    figure,
    for i = 1:size(x,1)-1
        for t = 1:size(v,1)-1
            plot(t,v(t+1,i)*v(t,i),'b*'),hold on
        end
    end
end