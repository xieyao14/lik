function [mu]= search_mu_timeonly(yob, kernel_sum, dT, M, initial_mu) % proj different structure with X kernel
%%
left_mu = initial_mu;
right_mu = left_mu+0.1; 
s_left = sum(yob./(1-exp(-dT*(left_mu+kernel_sum)))-1, "all");
s_right = sum(yob./(1-exp(-dT*(right_mu+kernel_sum)))-1, "all");

%
while s_left<0 || s_right>0
    if s_left<0 % monotone decreasing
        left_mu = left_mu*0.5;
        right_mu = left_mu+0.1;
    elseif s_right>0
        left_mu = left_mu+0.1;
        right_mu = left_mu+0.1;
    end
    s_left = sum(yob./(1-exp(-dT*(left_mu+kernel_sum)))-1, "all");
    s_right = sum(yob./(1-exp(-dT*(right_mu+kernel_sum)))-1, "all");
end

% assert:

assert(s_left>0 && s_right<0)

%% bisection
while abs(left_mu-right_mu)>1e-5
    next_mu = (left_mu+right_mu)/2;
    res = sum(yob./(1-exp(-dT*(next_mu+kernel_sum)))-1, "all")*dT/M;
    
    if res>0
        left_mu = next_mu;
    else 
        right_mu = next_mu;
    end
end
mu = left_mu;












