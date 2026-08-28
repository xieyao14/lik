function [Pob,Eob,eventJall,eventJob]= network_proj_mod_speed_N_graph2(edata, X, N, d, K) % proj different structure with X kernel
%%
Nob = N-d;
[batch_size,~] = size(edata);
Xbar = zeros(N,d,K*K); %including when event index is N
Xbar(:,:,:) = X;

Pob = cell(K,K);
Eob = cell(K,K);
eventJob = cell(K,1);
eventJall = [];


for k2 = 1:K
    edata_k2 = edata(:,(k2-1)*N+1:k2*N);
    num_event_k2 = sum(edata_k2,"all");

    event_i = find(edata_k2' > 0);
    [event_k2,col] = ind2sub([N,batch_size],event_i);

    

    I = kron(event_k2 , ones(d,1));
    J = I + kron(ones(num_event_k2,1), (1:d)');

    idxout = find(J > N | J <= d); 

    event_k2_ob = event_k2(event_k2>d)-d;
    obeventJ = zeros([length(event_k2_ob),1]); % pre-allocate and delete extra rows in the end
    row_start = 1;
    incre_afterd = sum(edata_k2(:,d+1:end),2);

    for i = 1:batch_size
        if incre_afterd(i)>0
            row_end = row_start+incre_afterd(i)-1;

            obeventJ(row_start:row_end) = event_k2_ob(row_start:row_end)+Nob*(i-1); 

            row_start = row_end+1;
        else
            continue
        end
    end

    row_start = 1;
    incre_all = sum(edata_k2,2);
    
    for i = 1:batch_size
        if incre_all(i)>0
            row_end = row_start+incre_all(i)-1;

            J(d*(row_start-1)+1:d*row_end) = J(d*(row_start-1)+1:d*row_end) -d+Nob*(i-1);

            row_start = row_end+1;
        else
            continue
        end
    end

    
    I(idxout)=[]; 
    J(idxout)=[]; 

    V_graph = zeros([d*num_event_k2,K]);
    for k1 = 1:K
        ik = k1+K*(k2-1);
        Xbar_ik = Xbar(:,:,ik); 
        V_graph(:,k1) = reshape(reshape(Xbar_ik(event_k2,:),[num_event_k2,d])', [d*num_event_k2,1]);
    end
    V_graph(idxout,:) = [];

    IJV_X_ob = [I,J,V_graph]; 

    
    eventJob{k2,1} = obeventJ;

    for k1 = 1:K
        Eob{k1,k2} = sparse(IJV_X_ob(:,1),IJV_X_ob(:,2),1,N,Nob*batch_size);
        Pob{k1,k2} = sparse(IJV_X_ob(:,1),IJV_X_ob(:,2),IJV_X_ob(:,2+k1),N,Nob*batch_size);
    end
end 








