function [value, isterminal, direction] = event_Secondary(t,X,prms)
value      = X(1) + prms.u - 1; % When x = 1-mu
isterminal = 1; % stops the integration
direction  = 0; % either direction 
end

