
function [tspan,x] = eul(tspan,x0,p)
%% EULER - Numerical ODE solver: Euler's method

  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;
  
  k=p(1); b=p(2);

  A=[0 1;-k -b];
  f=@(tspan, x) A*x;
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Step
    x(:,k) = x(:,k-1) + f(tspan(k-1),x(:,k-1))*dt;
    
  end
  x=x';
end