function [P,E,eventJall]= network_proj(edata, X, N, d, K) % proj different structure with X kernel
%%
[batch_size,~] = size(edata);
Xbar = zeros(N,d,K*K); %including when event index is N
% [trncN,~,~] = size(X);
% Xbar(1:trncN,:,:) = X;
Xbar(1:N-1,:,:) = X;

P = cell(K,K);
E = cell(K,K);
eventJall = cell(K,1);

for k1 = 1:K
    for k2 = 1:K
        ik = k1+K*(k2-1);
        IJV_X = []; %total length is (total number of events on node k1*d, 
        eventJ = []; %before removing trianglar outindex
        for i = 1:batch_size
            event_ik_i = find(edata(i,(k2-1)*N+1:k2*N) > 0)';  %interval index of events
            eventJ = [eventJ; event_ik_i+N*(i-1)];

            num_event = numel(event_ik_i);
            Xbar_ik = Xbar(:,:,ik); 
            I = kron(event_ik_i , ones(d,1));
            J = I + kron(ones(num_event,1), (1:d)');
            V = reshape(reshape(Xbar_ik(event_ik_i,:),[num_event,d])', [d*num_event,1]);
            
            idxout = find(J > N); 
            I(idxout)=[]; 
            J(idxout)=[]; 
            V(idxout)=[]; 
            IJV_X = [IJV_X; [I,J+N*(i-1),V]]; 
            
        end
        P{k1,k2} = sparse(IJV_X(:,1),IJV_X(:,2),IJV_X(:,3),N,N*batch_size);
        E{k1,k2} = sparse(IJV_X(:,1),IJV_X(:,2),1,N,N*batch_size);
        eventJall{k2,1} = eventJ;
    end
end






