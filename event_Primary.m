function [value, isterminal, direction] = event_Primary(t,X,prms)
value      = X(1) + prms.u; % When x = -mu
isterminal = 1; % stops the integration
direction  = 0; % either direction 
end
