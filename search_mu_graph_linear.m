function [mu]= search_mu_graph_linear(yob,ybarob, Pob, Lambar, h, M, initial_mu) % proj different structure with X kernel
%%
% Pob can be used to compute kernel sum, for example the kernel sum for
% node k1 from k2 can be computed by sum(P{k1,k2},1)
[V,~] = size(yob);
mu = zeros([V,1]);
% dldLambda = yob./Lambda + kron( ybarob.*(h./(exp(h*Lambar)-1)-1./Lambar) - h*(1-ybarob), ones(V,1) ); 
% dldLambda: (V,N*batch_sz)

constant_bar_part = sum(ybarob.*(h./(exp(h*Lambar)-1)-1./Lambar) - h*(1-ybarob),"all");

for k = 1:V
    %%
    left_mu = initial_mu(k);
    right_mu = left_mu+0.1; 

    % kernel_sum
    kernel_sum = zeros(1,length(Lambar));
    for k2 = 1:V
        kernel_sum(1,:) = kernel_sum(1,:) + sum(Pob{k,k2},1);
    end
    

    s_left = sum(yob(k,:)./(left_mu+kernel_sum), "all")+constant_bar_part;
    s_right = sum(yob(k,:)./(right_mu+kernel_sum), "all")+constant_bar_part;


    % 
    while s_left<0 || s_right>0
        if s_left<0 % monotone decreasing
            left_mu = left_mu*0.5;
            right_mu = left_mu+0.1;
        elseif s_right>0
            left_mu = left_mu+0.1;
            right_mu = left_mu+0.1;
        end
        s_left = sum(yob(k,:)./(left_mu+kernel_sum), "all")+constant_bar_part;
        s_right = sum(yob(k,:)./(right_mu+kernel_sum), "all")+constant_bar_part;
    end
    
    % assert:
    if isnan(s_left) || isnan(s_right)
%         fprintf("Unable to find root mu on data\n")
        mu = initial_mu;
        return;
    end

    assert(s_left>0 && s_right<0)
    
    %% bisection
    while abs(left_mu-right_mu)>1e-5
        next_mu = (left_mu+right_mu)/2;
        res = sum(yob(k,:)./(right_mu+kernel_sum), "all")+constant_bar_part;
        res = res*h/M;
        
        if res>0
            left_mu = next_mu;
        else 
            right_mu = next_mu;
        end
    end

    mu(k) = left_mu;
end












