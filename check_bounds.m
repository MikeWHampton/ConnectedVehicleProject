function [safe_bound,non_conservative_bound] = check_bounds(hdwysim,velsim,kappa_min,kappa_max,hst_min,hst_max)

% This code helps you check if the safe and non-conservative bounds are
% satified

if hdwysim >= velsim/kappa_max+hst_min
    
    safe_bound =  'Safe bound is satisfied.';
    
else
    
    safe_bound = 'Safe bound is violated!!';
    
end

if hdwysim <= velsim/kappa_min+hst_max
    
    non_conservative_bound =  'Non-conservative bound is satisfied.'
    
else
    
    non_conservative_bound = 'Non-conservative bound is violated!!'
    
end

end

