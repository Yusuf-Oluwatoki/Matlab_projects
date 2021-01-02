function ss = forced_spring_ss(p,data,s)
% sum-of-squares for forced spring-mass system.
% Function takes in the spring system parameters v and gamma
% as a vector p, data as the time and system values and s as the soltion
% method.
%
% 

time = data.ydata(:,1); %time from the previous system
Aobs = data.ydata(:,2); % system values
y0   = data.y0;         %initial values for the solution


[~,y]=forced_spring_ode(time,y0,p,s);  %function to solved the ode 
                                        %with new parameters p
Amodel = y(:,1);            %keep new system values

ss = sum((Aobs-Amodel).^2); %sum of squares error 
                            %of the new system and previous system
end

