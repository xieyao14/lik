function [Pob,Eob,eventJall,eventJob]= network_proj_mod(edata, X, N, d, K) % proj different structure with X kernel
%%
Nob = N-d;
[batch_size,~] = size(edata);
Xbar = zeros(N,d,K*K); %including when event index is N
Xbar(1:N-1,:,:) = X;

Pob = cell(K,K);
Eob = cell(K,K);
eventJall = cell(K,1);
eventJob = cell(K,1);

for k1 = 1:K
    for k2 = 1:K
        ik = k1+K*(k2-1);
        IJV_X = []; %total length is (total number of events on node k1*d, 
        IJV_X_ob = [];
        eventJ = []; %before removing trianglar outindex
        obeventJ = [];
        for i = 1:batch_size

            event_ik_i = find(edata(i,(k2-1)*N+1:k2*N) > 0)';  %interval index of events
            event_ik_i_ob = event_ik_i(event_ik_i>d)-d;
            eventJ = [eventJ; event_ik_i+N*(i-1)];
            obeventJ = [obeventJ; event_ik_i_ob+Nob*(i-1)];

            num_event = numel(event_ik_i);
            Xbar_ik = Xbar(:,:,ik); 
            I = kron(event_ik_i , ones(d,1));
            J = I + kron(ones(num_event,1), (1:d)');
            V = reshape(reshape(Xbar_ik(event_ik_i,:),[num_event,d])', [d*num_event,1]);
            
            idxout = find(J > N | J <= d); 
            I(idxout)=[]; 
            J(idxout)=[]; 
            V(idxout)=[]; 
            IJV_X = [IJV_X; [I,J+N*(i-1),V]]; 
            IJV_X_ob = [IJV_X_ob; [I,J-d+Nob*(i-1),V]]; 
            
        end
        Pob{k1,k2} = sparse(IJV_X_ob(:,1),IJV_X_ob(:,2),IJV_X_ob(:,3),N,Nob*batch_size);
        Eob{k1,k2} = sparse(IJV_X_ob(:,1),IJV_X_ob(:,2),1,N,Nob*batch_size);
        eventJall{k2,1} = eventJ;
        eventJob{k2,1} = obeventJ;
    end
end






