function eff = kernel_timeonly_highrank(tau, lag)
%%
eff = zeros(size(tau));

base_const = 0.02;
eff = eff + base_const*exp(-(8.*(lag).^2/25));

for rk = 1:15
    eff = eff + 0.3*2^(-rk)*(0.6+(cos(2+((tau-9-lag)/15)*1.3*(rk+1)*pi)))*exp(-(8.*(lag).^2*rk^2/25));
    
end

end