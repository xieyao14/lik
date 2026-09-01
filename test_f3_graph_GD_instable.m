% this is eliminate mu on graph with V = 5


clear all; rng(2024);
scriptDir = string(fileparts(mfilename("fullpath")));
addpath(scriptDir);
thein = fullfile(scriptDir, "Input") + string(filesep);
theout = fullfile(scriptDir, "Output") + string(filesep);
theplot = fullfile(scriptDir, "Plots") + string(filesep);
requiredFolders = [thein, theout, theplot];
for folder = requiredFolders
    if ~isfolder(folder), mkdir(folder); end
end
thedoc = "test_f3_graph_GD_instable";

%% parameters
V = 5;
mu_true = [0.3051; 0.3208; 0.2791; 0.3011; 0.3393]; % hardcode

ut = 8; % max time length T+d
ut_kernel = 4;
N = 32; %time discretization subinterval numbers
Nprime = 8; % memory depth
h = ut/(N+Nprime); % grid length

%% network graph
edges = [2,1;1,2;2,3;5,4;1,5;4,3;2,2;5,5;];
adj = zeros(V,V);
for i = 1:size(edges,1)
    adj(edges(i,1),edges(i,2))=1;
end
G = digraph(adj);

%%
if 1
    figure(1)
    p = plot(G,'Layout','circle','ArrowSize',5,'LineWidth',1);
    box off;
    axis off;
    p.NodeFontSize = 8;
    ax = gca;
    figH = gcf;
    set(figH, 'Units', 'points','OuterPosition', [0 0 157 180]) % standard size: 19.7 17.5
    exportgraphics(ax,strcat(theplot,thedoc,"_Graph",".pdf"));
end

%% discrete time kernel
disp("generating kernel...")
tic,
nt = (N+Nprime)*4;
ntau = (N+Nprime)*4;
t_grid = ((0.5:1:nt)/nt)'*ut_kernel;
tau_grid = ((0.5:1:nt)/nt)'*ut_kernel;

psi0 = zeros(nt, ntau, V^2);
for j = 1:ntau
    for k1 = 1:V
        for k2 = 1:V
            psi0(:,j,(k2-1)*V+k1) = truf_test141(k1,k2,t_grid, tau_grid(j));
        end
    end
end

downfactor = nt/(N+Nprime);
psi = zeros(nt/downfactor,ntau/downfactor,V^2);
for k = 1:V^2
    psi(:,:,k) = matrix_downsize(psi0(:,:,k), downfactor);
end
    

psi = psi(:,1:Nprime,:);
toc

%% kernel plot
figure(101);
Psimat = zeros(V*(N+Nprime),V*Nprime); % influence k2->k1
for k2 = 1:V
    rows = (k2-1)*(N+Nprime);
    for k1 = 1:V
        cols = (k1-1)*Nprime;
        Psimat(rows+1:rows+N+Nprime,cols+1:cols+Nprime) = psi(:,1:Nprime,k1+(k2-1)*V);
    end
end
fig = imagesc(Psimat);  
colorbar();
set(gca,'FontSize',18)
drawnow();
hold on;
for k1 = 1:V-1
    y = k1*(N+Nprime)*ones(2,1)+0.5;
    x = [0,V*Nprime]+0.5;
    p = plot(x,y,'LineWidth',1); % horizontal line
    p.Color = [0.2,0.2,0.2,0.4];
    y = [0,V*(N+Nprime)];
    x = k1*(Nprime)*ones(2,1)+0.5;
    p = plot(x,y,'LineWidth',1); % vertical line
    p.Color = [0.2,0.2,0.2,0.4];
end

tmp = G.Edges.EndNodes;
for i = 1:length(tmp)
    k2 = tmp(i,1);
    k1 = tmp(i,2);
    line([0,Nprime]+(k1-1)*(Nprime)*ones(1,2)+0.5, [Nprime,0]+(k2-1)*(N+Nprime)*ones(1,2)+0.5, 'Color', 'w','LineWidth',0.5);
    line([0,Nprime]+(k1-1)*(Nprime)*ones(1,2)+0.5, [N+Nprime,N]+(k2-1)*(N+Nprime)*ones(1,2)+0.5, 'Color', 'w','LineWidth',0.5);
end

xticks([[1:V]*Nprime-Nprime/2])
xticklabels({'1','2','3','4','5'})
xtickangle(0)
yticks([1:V]*(N+Nprime)-N/2)
yticklabels({'1','2','3','4','5'})
xlabel('{Node} $u$','fontsize',8,'interpreter','latex');
ylabel('{Node} $u^\prime$','fontsize',8,'interpreter','latex');
hold off;
% ax = gca;
% set(gca,'FontSize',8)
% figH = gcf;
% set(figH, 'Units', 'points','OuterPosition', [0 0 157 190]) % standard size: 19.7 17.5
% hold off;
% ax = gca;
% exportgraphics(ax,strcat(theplot,thedoc,"_TrueKer",".pdf"));

%% K Psi transfer
K = zeros([N+Nprime,N,V^2]);
% mask
a= ones(N+Nprime, Nprime);
psi_mask = triu(a, -N+1).*tril(a,0);
psi_mask = psi_mask(:,end:-1:1);
[im, jm, ~] = find(psi_mask);   %the nonzero indexing in the psi mask
indm = sub2ind( [N+Nprime, Nprime], im, jm);
jk = im+jm-Nprime;
for ik = 1:V^2
    psi_ik = psi(:,:,ik).*psi_mask;
    vpsi = psi_ik(indm); %[im,jm,vpsi] is the sparse repn of psi
    K(:,:,ik) = sparse(im,jk,vpsi, N+Nprime,N);
end

%% K visualization
figure(102);
Kmat = zeros(V*(N+Nprime),V*N);
for k2 = 1:V
    rows = (k2-1)*(N+Nprime);
    for k1 = 1:V
        cols = (k1-1)*N;
        Kmat(rows+1:rows+N+Nprime,cols+1:cols+N) = K(:,:,k1+(k2-1)*V);
    end
end
fig = imagesc(Kmat);  
colorbar();
set(gca,'FontSize',18)
drawnow();
hold on;



for k1 = 1:V-1
    y = k1*(N+Nprime)*ones(2,1)+0.5;
    x = [0,V*N]+0.5;
    p = plot(x,y,'LineWidth',1); % horizontal line
    p.Color = [0.2,0.2,0.2,0.4];
    y = [0,V*(N+Nprime)];
    x = k1*(N)*ones(2,1)+0.5;
    p = plot(x,y,'LineWidth',1); % vertical line
    p.Color = [0.2,0.2,0.2,0.4];
end


%% generate data
Phi_func = @(x,xbar) (1-exp(-xbar))./xbar.*x;


M = 50000; %40000;
y_ob = false(M, (Nprime+N)*V);
y_ob_3D = false(M,Nprime+N,V); 
% for i = 1:M 
% reshape(y_ob_3D(i,:,:),[1,(Nprime+N)*V])
% end
lambda_true = zeros(M,N*V);
lambda_true_3D = zeros(M,N,V);
lambar_true = zeros(M,N*V);
lambar_true_3D = zeros(M,N,V);

disp('generating trajectories...')
tic,
for i = 1: Nprime+N
    if i <= Nprime
        
        lambdat = zeros([M,V]);
        for k1 = 1:V
            for k2 = 1:V
                ik = k1+V*(k2-1);
                ypre = y_ob_3D(:, 1:i-1,k2); 
                kernelt = K(1:i-1,i,ik);
                lambdat(:,k1) = lambdat(:,k1) + sum(bsxfun(@times, kernelt, ypre' ),1)';
            end
            lambdat(:,k1) = lambdat(:,k1) + mu_true(k1);
        end
        lambart = sum(lambdat,2);
        lambart = kron(lambart, ones(1,V));

        pyt = Phi_func( h*max(0,lambdat), h*max(0,lambart));
%         y_ob_3D(:, i, :) = (rand(M,V) <  pyt);

        pyt_extend = [pyt, 1-sum(pyt, 2)];% (M,V+1)
        now_event = mnrnd(1,pyt_extend);
        y_ob_3D(:, i, :) = now_event(:,1:V);
  
    else
        t = i-Nprime; % the most historical time
        lambdat = zeros([M,V]);
        for k1 = 1:V
            for k2 = 1:V
                ik = k1+V*(k2-1);
                ypre = y_ob_3D(:, t:t+Nprime-1,k2); 
                kernelt = K(t:t+Nprime-1,t,ik);
                lambdat(:,k1) = lambdat(:,k1) + sum(bsxfun(@times, kernelt, ypre' ),1)';
            end
            lambdat(:,k1) = lambdat(:,k1) + mu_true(k1);
        end
        if min(lambdat,[],"all") < 0
            warning( sprintf('min lamdba = %6.4f\n',min(lambdat,[],"all")))
        end
        lambart = sum(lambdat,2);
        

        lambart = kron(lambart, ones(1,V));

        pyt = Phi_func( h*max(0,lambdat), h*max(0,lambart));
        pyt_extend = [pyt, 1-sum(pyt, 2)];% (M,V+1)
        now_event = mnrnd(1,pyt_extend);
        y_ob_3D(:, i, :) = now_event(:,1:V);
%         y_ob_3D(:, i, :) = (rand(M,V) <  pyt);

        lambda_true_3D(:,t,:) = lambdat; 
        lambar_true_3D(:,t,:) = lambart; 
        
    end
end

for i = 1:M 
    y_ob(i,:) = reshape(y_ob_3D(i,:,:),[1,(Nprime+N)*V]);
    lambda_true(i,:) = reshape(lambda_true_3D(i,:,:),[1,N*V]);
    lambar_true(i,:) = reshape(lambar_true_3D(i,:,:),[1,N*V]);
end
toc


event = y_ob;
save(strcat(thein,thedoc,"_Ber.mat"), "event");

%%
tmp = min(lambda_true,[],2);
nega_idx = find(tmp<=0)
lambda_true(nega_idx,:) = [];
[min_lam_true,i_min]= min(min(lambda_true,[],2))
max_lam_true= max(lambda_true(:))

figure(4);clf;
imagesc(lambda_true); colorbar;
title(sprintf('true lambda, min=%5.4f, max=%5.4f',min_lam_true,max_lam_true )); 

figure(5),clf;
imagesc(y_ob);
title('y observed')

figure(6), clf; hold on;
plot(lambda_true(3:5:20,:)')
plot(lambda_true(i_min,:)','.-')
grid on;



% lambda in (0.04,2.32)




%% read in data
% event = load(strcat(thein,"test76_graph_LIK_10","_Ber_clean.mat"));
% event_data = event.event;
% event_data = event_data(1:M,:);
event(nega_idx,:) = [];

%% split data
event_data = event(1:M,:);
ntr = M*0.8; %32000; %16000;
nte = min(500,M-ntr);

tmp = randperm(M);
idx_tr = tmp(1:ntr);
idx_te = tmp(ntr+1:ntr+nte);


event_data_tr = event_data(idx_tr,:);
event_data_te = event_data(idx_te,:);
            

% disp('compute true lambda...')
% tic,
% % compute projected data
% [Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N(event_data, psi(:,1:Nprime,:), N+Nprime, Nprime, V);
% 
% % compute Lambda, size (V, N*M)
% Lambda_true = zeros(V,N*M);
% for k1 = 1:V
%     for k2 = 1:V
%         Lambda_true(k1,:) = Lambda_true(k1,:) + sum(Pob{k1,k2},1);
%     end
%     Lambda_true(k1,:) = Lambda_true(k1,:)+mu_true(k1);
% end
% Lambar_true = sum(Lambda_true,1); % size (1, N*batch_size)
% 
% lam_true = zeros(M, N*V); %[M, N*V]
% lambar_true = zeros(M, N*V); %[M, N*V]
% for b = 1:M
%     lam_true(b,:) = reshape((Lambda_true(:,(b-1)*N+1:b*N))',[N*V,1]);
%     lambar_true(b,:) = kron(ones(1,V),Lambar_true(1,(b-1)*N+1:b*N));
% end

lam_true_val = lambda_true(idx_te,:);
lambar_true_val = lambar_true(idx_te,:);
p_val_true = (1-exp(-h.*lambar_true_val))./lambar_true_val.*lam_true_val;
% toc


%% 
warm_start = 0;
% zero init

if warm_start
    disp("warm start\n")
    X = load(strcat(theout,"test135_graph_fixmuGD_Epoch40.mat"));
    X = X.X;
    theta_psi = X;
else
    theta_psi = zeros(Nprime+N, Nprime, V^2); %initial value
    X = zeros(Nprime+N, Nprime, V^2); %initial matrix
end

% rank 1 init
%theta_psi = psi*v1(:,1)*v1(:,1)';

%% kernel recovery
use_VI = 0;
% barrier hyper parameter
min_b = 0.03; %0.03;
B_weight = 0.1; % 0.1
mode = "eliminate";

%etaK = 0.1; %learning rate
smooth_weight = 0.1; % todo: study the relationship btw lr and smooth_weight
% smooth_weight should be related to dimension


if use_VI
    label = "VI";
    lr_schedule = [0.4*ones(50,1)];
else
    label = "GD";
    lr_schedule = [0.4*ones(50,1)];
end

bs_schedule = [800*ones(50,1)] ;

num_epoch = numel(lr_schedule );

nll_all = zeros(num_epoch,1);
err1_kernel_all = zeros(num_epoch,1);
p_mae_all = zeros(num_epoch,2);
mu_mae_all = zeros(num_epoch,1);
mu_all = zeros(num_epoch,V);

if mode == "cheat"
    mu = mu_true;
elseif mode == "eliminate"
    mu = zeros(V,1); % initial base
    for ik = 1:V
        mu(ik) = mean(mean(event_data_tr(:,(ik-1)*(N+Nprime)+1:ik*(N+Nprime)),2))/h;
    end
end
%%

for iepoch = 1:num_epoch

    etaK = lr_schedule(iepoch);
    batch_size = bs_schedule(iepoch);
    num_batch = floor(ntr/batch_size);

    disp('training...')
    tic,

    count_violates = 0; 

    nll_sum = 0;
    count_type2 = 0;

    idtr_epoch = randperm(ntr); % scrable training data
    for ibatch = 1: num_batch
%         ibatch 
        idtr_batch = idtr_epoch((ibatch-1)*batch_size+1:ibatch*batch_size);
        %%
        edata = event_data_tr(idtr_batch,:); 
             
        %% compute projected data
%         tic
%         [Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N(edata, X, N+Nprime, Nprime, V);
%         toc
% 
%         %%
%         tic
%         [Pob_g, Eob_g, eventJall,eventJob_g]= network_proj_mod_speed_N_graph(edata, X, N+Nprime, Nprime, V);
%         toc

        %%
        [Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N_graph2(edata, X, N+Nprime, Nprime, V);
        

        %% compute Lambda, size (V, N*batch_size)
        
        Lambda = zeros(V,N*batch_size);
        for k1 = 1:V
            for k2 = 1:V
                Lambda(k1,:) = Lambda(k1,:) + sum(Pob{k1,k2},1);
            end
            Lambda(k1,:) = Lambda(k1,:)+mu(k1);
        end
        Lambar = sum(Lambda,1); % size (1, N*batch_size), sum over all locations

        % decide if Lambda violates buffer
        lam_batch = zeros(batch_size, N*V); %[batch_size, N]
        for b = 1:batch_size
            lam_batch(b,:) = reshape((Lambda(:,(b-1)*N+1:b*N))',[N*V,1]);
        end

%         lam_batch_speed = lam_batch;


        min_lam_batch = min( lam_batch, [],2);
        id_type1 = find( min_lam_batch < min_b );
        batch_sz1 = numel(id_type1);
        if  batch_sz1 > 0
            min_min_lam_batch = min( min_lam_batch);
            warning(sprintf('has violation, min lam = %6.4f\n', min_min_lam_batch));
%             assert( min_min_lam_batch > 0 ); %still ask min to be >0
            if min_min_lam_batch < 0
                warning('min lam less than zero!')
            end
        end

        id_type2 = (1:batch_size);
        id_type2(id_type1) = [];

        id_type1_squeeze = zeros([1,N*batch_sz1]); % the index of type1 events in ybarob
        for bb = 1:batch_sz1
            id = id_type1(bb);
            id_type1_squeeze((bb-1)*N+1:bb*N) = (id-1)*N+1:id*N;
        end
        id_type1_K = kron(ones(V,1),id_type1_squeeze); % the index of type1 events in yob

        batch_sz2 = batch_size-batch_sz1;
        count_violates = count_violates + batch_sz1;


        % compute log-likelihood on type2
        yob = zeros(V,N*batch_size); % reshape of edata, boolean matrix size (V, N*batch_size)
        for k1 = 1:V
            yob(k1,eventJob{k1,1}) = 1;
        end
        ybarob = sum(yob,1);


        yob_type1 = yob(:,id_type1_squeeze);
        yob_type2 = yob;
        yob_type2(:,id_type1_squeeze) = [];

        ybarob_type1 = ybarob(id_type1_squeeze);
        ybarob_type2 = ybarob;
        ybarob_type2(id_type1_squeeze) = [];

        Lambda_type1 = Lambda(:,id_type1_squeeze);
        Lambda_type2 = Lambda;
        Lambda_type2(:,id_type1_squeeze) = [];

        Lambar_type1 = Lambar(id_type1_squeeze);
        Lambar_type2 = Lambar;
        Lambar_type2(id_type1_squeeze) = [];

        lklhA1 = sum(ybarob_type2.*log( (1-exp(-h*Lambar_type2))./Lambar_type2 ),"all");
        lklhA2 = sum(yob_type2.*log(Lambda_type2),"all");
        lklhB = sum((1-ybarob_type2).*Lambar_type2,"all")*h;
        likelihood = (lklhA1+lklhA2-lklhB);

        nll_sum = nll_sum-likelihood;
        count_type2 = count_type2 +batch_sz2; %number of trajectories that is added into nll
        


        %% compute grad field for F on type2
  
        if use_VI % negative vector field to minimize the loss
            dldLambda = -(kron( (1-exp(-h*Lambar))./Lambar, ones(V,1)).*Lambda - yob); % (V,N*batch_sz2) 
        else % negative vector field to minimize the loss
            dldLambda = yob./Lambda + kron( ybarob.*(h./(exp(h*Lambar)-1)-1./Lambar) - h*(1-ybarob), ones(V,1) ); % (V,N*batch_sz2)
        end
     

        %%
%         dldlam_batch = zeros(batch_sz2, N*1); %[batch_size, N]
%         for b = 1:batch_sz2
%             dldlam_batch(b,:) = reshape((dldLambda(:,(b-1)*N+1:b*N))',[N*1,1]);
%         end
        
        dldX = zeros(N+Nprime,Nprime,V^2); % gardient on kernel matrix 
        diag_ind = -(1:Nprime);
        % compute grad from dK
        for k1 = 1:V
            tmp = dldLambda(k1,:); % k1-th row of dl/dLambda
            tmp_reshaped = reshape(tmp, [1,N,batch_size]);
            tmp_reshaped = tmp_reshaped(:,:,id_type2);
            for k2 = 1:V
               
                %% compute grad from dK
                Eobk1k2 = full(Eob{k1,k2});
                E_reshaped = reshape(Eobk1k2, [N+Nprime,N,batch_size]);
                E_reshaped = E_reshaped(:,:,id_type2);
                
                dldK_k1k2 = sum(E_reshaped.*tmp_reshaped,3);

%                 tmp1 = dldK_k1k2; 


%                 %% compute grad from dK
%                 dldK_k1k2 = zeros(N+Nprime,N);
%                 for ib = id_type2
%                     % sum over batch
%                     dldK_k1k2 = dldK_k1k2 + Eobk1k2(:,N*(ib-1)+1:N*ib).*kron(tmp(N*(ib-1)+1:N*ib), ones(N+Nprime,1));
%                 end
%                 tmp2 = dldK_k1k2;
                
           
                
                %% clear dldK_k1k2_alongbatch;
                
                dldK_k1k2_extend = zeros(N+Nprime,N+Nprime);
                dldK_k1k2_extend(:,Nprime+1:N+Nprime) = dldK_k1k2;
                

                %% reshape dK to dPhi
                
                dldX_k1k2 = spdiags(dldK_k1k2_extend',diag_ind);
                dldX(:,:,k1+V*(k2-1)) = dldX_k1k2/batch_sz2; 
                
            end                
        end
       
        



        %% compute grad field for barrier on type1
        
        dBdX = zeros(N+Nprime,Nprime,V^2);
        
    
        if batch_sz1 > 0
            if 0
                dBdLambda = -min_b./(Lambda);
            else
                dBdLambda = (Lambda-min_b)/(min_b*0.1);
            end
            
            dBdLambda = dBdLambda.*(Lambda < min_b);
            for k1 = 1:V
                tmp = dBdLambda(k1,:); % k1-th row of dB/dLambda
                for k2 = 1:V
                    Eobk1k2 = Eob{k1,k2};
                    % compute grad from dK

                    %%% old
                    % dBdK_k1k2 = zeros(N+Nprime,N);
                    % for ib = id_type1
                    %     % sum over batch
                    %     dBdK_k1k2 = dBdK_k1k2 + Eobk1k2(:,N*(ib-1)+1:N*ib).*kron(tmp(N*(ib-1)+1:N*ib), ones(N+Nprime,1));
                    % end
                    E_reshaped = reshape( ...
                        full(Eobk1k2), [N+Nprime, N, batch_size]);

                    tmp_reshaped = reshape( ...
                        tmp, [1, N, batch_size]);

                    dBdK_k1k2 = sum( ...
                        E_reshaped(:,:,id_type1) .* ...
                        tmp_reshaped(:,:,id_type1), 3);

                    % clear dBdK_k1k2_alongbatch
                    % reshape dK to dPhi
                    dBdK_k1k2_extend = zeros(N+Nprime,N+Nprime);
                    dBdK_k1k2_extend(:,Nprime+1:N+Nprime) = dBdK_k1k2;

                    dBdX_k1k2 = spdiags(dBdK_k1k2_extend',diag_ind);
                    dBdX(:,:,k1+V*(k2-1)) = dBdX_k1k2; 
                end                
            end
            dBdX = dBdX*B_weight; 
        end

        

        
       
       
        
        %% update kernel 
        dX = dldX-dBdX;  
        X = X + etaK*dX; %update X

 
        

        % tSVD to enforce low-rank: only apply after burn in
        if (0)
            Xconc = zeros((N+Nprime)*V^2, Nprime); %initial matrix
            for k = 1:V^2
                Xconc(((N+Nprime)*(k-1)+1):(N+Nprime)*k,:) = X(:,:,k);
            end

            [u,s,v] = svd(Xconc,'econ');
            svec = diag(s);
%             rX = max(1,sum(svec>tau));
            rX = 3;
            Xconc = u(:,1:rX)*s(1:rX,1:rX)*v(:,1:rX)';

            for k = 1:V^2
                X(:,:,k) = Xconc(((N+Nprime)*(k-1)+1):(N+Nprime)*k,:);
            end

        end

        %% update mu under eliminate mode
        if mode == "eliminate"
            mu_proposed = search_mu_graph_linear(yob, ybarob, Pob, Lambar, h, batch_size, mu);
            mu = 0.9*mu+0.1*mu_proposed;
            mu_all(iepoch,:) = mu;
        end

    end

    % smoothify kernel
    Lt = make_1d_laplacian(Nprime+N);
    Ltau = make_1d_laplacian(Nprime);

    for ik = 1:V^2
        Xk = X(:,:,ik);
        Xk = Xk- etaK*smooth_weight*(Lt*Xk);


        Xk = Xk';
        Xk = Xk- etaK*smooth_weight*(Ltau*Xk);
        Xk = Xk';

        X(:,:,ik) = Xk;
    end
    theta_psi = X;


    if mod(iepoch,10)==0
        
        Xmat = zeros(V*(N+Nprime),V*Nprime);
        for k2 = 1:V
            rows = (k2-1)*(N+Nprime);
            for k1 = 1:V
                cols = (k1-1)*Nprime;
                Xmat(rows+1:rows+N+Nprime,cols+1:cols+Nprime) = X(:,:,k1+(k2-1)*V);
            end
        end


        figure(9),clf;
        fig = imagesc(Xmat); 

        colorbar();
        fig.AlphaData = 0.8;
        set(gca,'FontSize',18)
        hold on;
        for k1 = 1:V-1
            y = k1*(N+Nprime)*ones(2,1);
            x = [0,V*Nprime];
            p = plot(x,y,'LineWidth',2);
            p.Color = [0.2,0.2,0.2,0.4];
            y = [0,V*(N+Nprime)];
            x = k1*(Nprime)*ones(2,1);
            p = plot(x,y,'LineWidth',2);
            p.Color = [0.2,0.2,0.2,0.4];
        end
        title(sprintf('estimated kernel'));
        set(gca,'XTick',[], 'YTick', []);
        hold off;
        drawnow();
    end
       
    %% kernel error on submatrix    
    true_norm = zeros(V,V,3);
    for k1 = 1:V
        for k2 = 1:V
            ik = k1+V*(k2-1);
            psi_sub = psi(Nprime+1:N,1:Nprime,ik);
            true_norm(k1,k2,:) = [norm(psi_sub,1);norm(psi_sub,2);norm(psi_sub,"inf")];
        end
    end
    
    [I,J] = find(true_norm(:,:,1)>0);
    idx_edges = [I,J];
    [I,J] = find(true_norm(:,:,1)==0);
    idx_nonedges = [I,J];
    l_edges = size(idx_edges,1);
    l_nonedges = size(idx_nonedges,1);

    %%
    true_norm_all = zeros([l_edges,3]);
    absolute_error_all =  zeros([l_edges,3]);
    for i = 1:l_edges
        k1 = idx_edges(i,1);
        k2 = idx_edges(i,2);
        ik = k1+V*(k2-1);
        kernel_diff = X(Nprime+1:N,1:Nprime,ik)-psi(Nprime+1:N,1:Nprime,ik);
        kernel_diff = reshape(kernel_diff, [Nprime*(N-Nprime),1]);
        psi_vec = reshape(psi(Nprime+1:N,1:Nprime,ik), [Nprime*(N-Nprime),1]);
        true_norm_all(i,:) = [norm(psi_vec,1);norm(psi_vec,2);norm(psi_vec,"inf")];
        absolute_error_all(i,:) = [norm(kernel_diff,1);norm(kernel_diff,2);norm(kernel_diff,"inf")];
    end
    err_12inf = mean(absolute_error_all,1)./mean(true_norm_all,1);
    err1 = err_12inf(2);
    err1_kernel_all(iepoch) = err1;
    

    % epoch average nll
    nll_epoch =nll_sum/count_type2;
    nll_all(iepoch) = nll_epoch;

    if mod(iepoch,10)==0

        figure(19),clf;
        plot( 1:iepoch, nll_all(1:iepoch), '.-');
        grid on; title('tr nll'); set(gca,'FontSize',15);
        drawnow();

 
    end

    % mu error
    errmu = norm(mu-mu_true,2)/norm(mu_true,2);
    mu_mae_all(iepoch) = errmu;

    toc,

    %% prediction error on test traj
    disp('testing...')
    tic,
    % compute projected data
    [Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N_graph2(event_data_te, X, N+Nprime, Nprime, V);
    
    % compute Lambda, size (V, N*nte)
    Lambda_val = zeros(V,N*nte);
    for k1 = 1:V
        for k2 = 1:V
            Lambda_val(k1,:) = Lambda_val(k1,:) + sum(Pob{k1,k2},1);
        end
        Lambda_val(k1,:) = Lambda_val(k1,:)+mu(k1);
    end
    Lambar_val = sum(Lambda_val,1); % size (1, N*nte)
    
    lam_val = zeros(nte, N*V); %[M, N*V]
    lambar_val = zeros(nte, N*V); %[M, N*V]
    for b = 1:nte
        lam_val(b,:) = reshape((Lambda_val(:,(b-1)*N+1:b*N))',[N*V,1]);
        lambar_val(b,:) = kron(ones(1,V),Lambar_val(1,(b-1)*N+1:b*N));
    end
    
    hatp_val = (1-exp(-h.*lambar_val))./lambar_val.*lam_val;







    min_lam_val = min(lam_val,[],2);
    if sum( min_lam_val < min_b ) > 0
        warining('has violation on val set');
%         assert( min(min_lam_val) > 0);
    end
    
    mae_lam_val = mean(abs( lam_val - lam_true_val),2);
    l1_lam_val_true = mean(abs(lam_true_val), 2);
    [~,idmae] = sort( mae_lam_val,'descend');

    mae_p_val = mean(abs( hatp_val - p_val_true),2);
    l1_p_val_true = mean( abs(p_val_true), 2);
    relmae_p_val = mae_p_val./l1_p_val_true;
    avg_p_mae = mean(relmae_p_val);
    worst_p_mae = max(relmae_p_val );

    p_mae_all(iepoch,:) = [avg_p_mae,worst_p_mae];
    toc
            
    if mod(iepoch,10)==0


        p_mae_all_GD_small = load(strcat(theout,thedoc,"GD","learning_rate_",string(0.2),"_pmae.mat"));
        p_mae_all_GD_small = p_mae_all_GD_small.p_mae_all;
        
        p_mae_all_VI = load(strcat(theout,thedoc,"VI","learning_rate_",string(0.4),"_pmae.mat"));
        p_mae_all_VI = p_mae_all_VI.p_mae_all;
        
        figure(20),clf;
        plot( 1:iepoch, p_mae_all(1:iepoch,1), '.-');
        hold on;
        plot(1:iepoch, p_mae_all_GD_small(1:iepoch,1), '-.');
        
        plot(1:iepoch, p_mae_all_VI(1:iepoch,1), '--');
        grid on; title('prediction error'); 
        set(gca,'FontSize',15);
        drawnow();
    end

    %%
    fprintf('epoch %d, count of viol=%d, nll_tr=%6.4f; err1=%6.4f, errp=%6.4f, %6.4f, errmu = %6.4f\n', ...
        iepoch, count_violates, nll_epoch, err1, avg_p_mae, worst_p_mae, errmu);

%     if mod(iepoch,10) == 0
%         save(strcat(theout,thedoc,label,"_Epoch",string(iepoch),".mat"), "X");
%         save(strcat(theout,thedoc,label,"_Epoch",string(iepoch),"_mu.mat"), "mu");
%     end



   
end

%%
save(strcat(theout,thedoc,label,"learning_rate_",string(etaK),"_pmae.mat"), "p_mae_all");




%% kernel
save(strcat(theout,thedoc,label,".mat"), "X");

%% kernel plot
Xmat = zeros(V*(N+Nprime),V*Nprime);
for k2 = 1:V
    rows = (k2-1)*(N+Nprime);
    for k1 = 1:V
        cols = (k1-1)*Nprime;
        Xmat(rows+1:rows+N+Nprime,cols+1:cols+Nprime) = X(:,:,k1+(k2-1)*V);
    end
end


figure(12),clf;
fig = imagesc(Xmat); 

colorbar();
fig.AlphaData = 0.8;
set(gca,'FontSize',18)
hold on;
for k1 = 1:V-1
    y = k1*(N+Nprime)*ones(2,1);
    x = [0,V*Nprime];
    p = plot(x,y,'LineWidth',2);
    p.Color = [0.2,0.2,0.2,0.4];
    y = [0,V*(N+Nprime)];
    x = k1*(Nprime)*ones(2,1);
    p = plot(x,y,'LineWidth',2);
    p.Color = [0.2,0.2,0.2,0.4];
end
% title(sprintf('estimated kernel'));
set(gca,'XTick',[], 'YTick', []);
tmp = G.Edges.EndNodes;
for i = 1:length(tmp)
    k2 = tmp(i,1);
    k1 = tmp(i,2);
    line([0,Nprime]+(k1-1)*(Nprime)*ones(1,2)+0.5, [Nprime,0]+(k2-1)*(N+Nprime)*ones(1,2)+0.5, 'Color', 'w','LineWidth',0.5);
    line([0,Nprime]+(k1-1)*(Nprime)*ones(1,2)+0.5, [N+Nprime,N]+(k2-1)*(N+Nprime)*ones(1,2)+0.5, 'Color', 'w','LineWidth',0.5);
end

xticks([[1:V]*Nprime-Nprime/2])
xticklabels({'1','2','3','4','5'})
xtickangle(0)
yticks([1:V]*(N+Nprime)-N/2)
yticklabels({'1','2','3','4','5'})
xlabel('{Node} $u$','fontsize',8,'interpreter','latex');
ylabel('{Node} $u^\prime$','fontsize',8,'interpreter','latex');
hold off;
ax = gca;
set(gca,'FontSize',8)
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 157 190]) % standard size: 19.7 17.5
hold off;
ax = gca;
% exportgraphics(ax,strcat(theplot,thedoc,"_EstKer",label,".pdf"));

%% kernel error
true_norm = zeros(V,V,3);
for k1 = 1:V
    for k2 = 1:V
        ik = k1+V*(k2-1);
        psi_sub = psi(Nprime+1:N,1:Nprime,ik);
        true_norm(k1,k2,:) = [norm(psi_sub,1);norm(psi_sub,2);norm(psi_sub,"inf")];
    end
end

[I,J] = find(true_norm(:,:,1)>0);
idx_edges = [I,J];
[I,J] = find(true_norm(:,:,1)==0);
idx_nonedges = [I,J];
l_edges = size(idx_edges,1);
l_nonedges = size(idx_nonedges,1);

% label = "GD"
% 
% X = load(strcat(theout,thedoc,label,".mat"));
% X = X.X;

% true_norm_all = zeros([l_edges,3]);
% absolute_error_all =  zeros([l_edges,3]);
for i = 1:l_nonedges
    k1 = idx_nonedges(i,1);
    k2 = idx_nonedges(i,2);
    ik = k1+V*(k2-1);
    X(:,:,ik) = 0;
end

psi_vec = reshape(psi(Nprime+1:N,1:Nprime,:),[],1);
kernel_diff = X-psi;
kernel_diff = reshape(kernel_diff(Nprime+1:N,1:Nprime,:),[],1);
true_norm_all = [norm(psi_vec,1);norm(psi_vec,2);norm(psi_vec,"inf")];
absolute_error_all = [norm(kernel_diff,1);norm(kernel_diff,2);norm(kernel_diff,"inf")];
kernel_rela_err = absolute_error_all./true_norm_all;
fprintf(['-----Kernel Error Table----- \n' ...
    'Line 1: Rel l1, Rel l2, Rel l inf\n'])
kernel_rela_err'
save(strcat(theout,thedoc,label,"_EstKerRelErr.mat"), "kernel_rela_err");

%% mu
save(strcat(theout,thedoc,label,"_mu.mat"), "mu");

mu_all_err_L2 = zeros([num_epoch,1]);
for i=1:num_epoch
    mu_all_err_L2(i) = norm(mu_true'-mu_all(i,:),2)/norm(mu_true,2);
end

figure(7),clf;

box on;
plot(1:iepoch, mu_all_err_L2(1:iepoch),'.-','LineWidth',1);
grid on;

xlabel('{Epoch}','interpreter','latex');
ylabel('{L2 relative error of $\mu$}','interpreter','latex');
ylim([0,0.2])
set(gca,'FontSize',8);
% lgd = legend('train','Location','best');
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,label,"_mu_all",".pdf")); 
save(strcat(theout,thedoc,label,"_mu_all.mat"), "mu_all");

%% mu error
mudiff = mu-mu_true;
mu_error = [norm(mudiff,1);norm(mudiff,2);norm(mudiff,"inf")];
mu_rela_err = mu_error./[norm(mu_true,1);norm(mu_true,2);norm(mu_true,"inf")];
fprintf(['-----Mu Error Table----- \n' ...
    'Line 1: Rel l1, Rel l2, Rel l inf\n'])
mu_rela_err'
save(strcat(theout,thedoc,label,"_EstMuRelErr.mat"), "mu_rela_err");

%% likelihood plot
figure(8),clf;

box on;
plot(1:iepoch, nll_all(1:iepoch),'.-','LineWidth',1);
grid on;

xlabel('{Epoch}','interpreter','latex');
ylabel('{Negative log-likelihood}','interpreter','latex');
set(gca,'FontSize',8);
lgd = legend('train','Location','best');
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,label,"_TrainLogLike",".pdf")); 
save(strcat(theout,thedoc,label,"_TrainLogLike",".mat"), "nll_all");

%% prob prediction error
X = load(strcat(theout,thedoc,label,".mat"));
X = X.X;

[Pob_te, Eob_te, eventJall_te,eventJob_te]= network_proj_mod_speed_N_graph2(event_data_te, X, N+Nprime, Nprime, V);
% compute Lambda, size (1, N*batch_size)



Lambda_te = zeros(V,N*nte);
for k1 = 1:V
    for k2 = 1:V
        Lambda_te(k1,:) = Lambda_te(k1,:) + sum(Pob_te{k1,k2},1);
    end
    Lambda_te(k1,:) = Lambda_te(k1,:)+mu(k1);
end
Lambar_te = sum(Lambda_te,1); % size (1, Nob*batch_size), sum over all locations
Lambar_te = kron(Lambar_te, ones(V,1));
Prob_te = (1-exp(-h.*Lambar_te))./Lambar_te.*Lambda_te;

% reshape Lambda to pair the shape of data y
prob_te = zeros(N*V, nte);
for b = 1:nte
    prob_te(:,b) = reshape((Prob_te(:,(b-1)*N+1:b*N))',[N*V,1]);
end

[Pob_te_tru, Eob_te, eventJall_te,eventJob_te]= network_proj_mod_speed_N_graph2(event_data_te, psi, N+Nprime, Nprime, V);
% compute Lambda, size (1, N*batch_size)
Lambda_te_tru = zeros(V,N*nte);
for k1 = 1:V
    for k2 = 1:V
        Lambda_te_tru(k1,:) = Lambda_te_tru(k1,:) + sum(Pob_te_tru{k1,k2},1);
    end
    Lambda_te_tru(k1,:) = Lambda_te_tru(k1,:)+mu_true(k1);
end

Lambar_te_tru = sum(Lambda_te_tru,1); % size (1, Nob*batch_size), sum over all locations
Lambar_te_tru = kron(Lambar_te_tru, ones(V,1));
Prob_te_tru = (1-exp(-h.*Lambar_te_tru))./Lambar_te_tru.*Lambda_te_tru;

% reshape Lambda to pair the shape of data y
prob_te_tru = zeros(N*V, nte);
for b = 1:nte
    prob_te_tru(:,b) = reshape((Prob_te_tru(:,(b-1)*N+1:b*N))',[N*V,1]);
end


prob_diff = prob_te-prob_te_tru;



proberror = zeros([nte,6]); 
for l = 1:nte
    proberror(l,1) = norm(prob_diff(:,l),1)/N; % MAE
    proberror(l,2) = sqrt(norm(prob_diff(:,l),2).^2/N); % RMSE
    proberror(l,3) = norm(prob_diff(:,l),"inf"); % worst

    proberror(l,4) = norm(prob_diff(:,l),1)/norm(prob_te_tru(:,l),1); % Relative l1
    proberror(l,5) = norm(prob_diff(:,l),2)/norm(prob_te_tru(:,l),2); % Relative l2
    proberror(l,6) = norm(prob_diff(:,l),"inf")/norm(prob_te_tru(:,l),"inf"); % Relative l inf
end



fprintf(['-----Error Table----- \n' ...
    'Line 1: MAE, RMSE, Worst; Rel l1, Rel l2, Rel l inf\n'])
mean(proberror,1)
save(strcat(theout,thedoc,label,"_ProbPredErr",".mat"), "proberror");


return;

%%
figure(8),clf;

nll_stable = load(strcat(theout,"test_f3_graph_eliminatemu","GD","_TrainLogLike",".mat"), "nll_all");
nll_stable = nll_stable.nll_all;



nll_stable_VI = load(strcat(theout,"test_f3_graph_eliminatemu","VI","_TrainLogLike",".mat"), "nll_all");
nll_stable_VI = nll_stable_VI.nll_all;


nll_all = load(strcat(theout,thedoc,"GD","_TrainLogLike",".mat"));
nll_all = nll_all.nll_all;




range = 5:50;
box on;
plot(range, nll_stable(range),'--','LineWidth',1,'Color',[0,0,1,0.8] );
grid on;
hold on;
plot(range, nll_stable_VI(range),'-.','LineWidth',1,'Color',[0.8,0.5,1,0.8] );
plot(range, nll_all(range),'.-','LineWidth',1,'Color',[1,0,0,0.8] );


xlabel('{Epoch}','interpreter','latex');
ylabel('{Negative log-likelihood}','interpreter','latex');
set(gca,'FontSize',8);
% lgd = legend('GD: $\gamma_k=0.2,k\le 50$, $\gamma_k=0.1$ otherwise', ...
%     'VI: $\gamma_k=0.4,k\le 50$, $\gamma_k=0.2$ otherwise', ...
%     'GD: $\gamma_k=0.4,k\le 50$, $\gamma_k=0.2$ otherwise', ...
%     'interpreter','latex','Location','best');

lgd = legend('GD: $\gamma_k=0.2$', ...
    'VI: $\gamma_k=0.4$', ...
    'GD: $\gamma_k=0.4$', ...
    'interpreter','latex','Location','best');
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,label,"_GDinstability",".pdf")); 



%%
figure(9),clf;

p_mae_all_GD_small = load(strcat(theout,thedoc,"GD","learning_rate_",string(0.2),"_pmae.mat"));
p_mae_all_GD_small = p_mae_all_GD_small.p_mae_all;

p_mae_all_GD_large = load(strcat(theout,thedoc,"GD","learning_rate_",string(0.4),"_pmae.mat"));
p_mae_all_GD_large = p_mae_all_GD_large.p_mae_all;


p_mae_all_VI = load(strcat(theout,thedoc,"VI","learning_rate_",string(0.4),"_pmae.mat"));
p_mae_all_VI = p_mae_all_VI.p_mae_all;




range = 5:50;
box on;
plot(range, p_mae_all_GD_small(range),'--','LineWidth',1,'Color',[0,0,1,0.8] );
grid on;
hold on;
plot(range, p_mae_all_VI(range),'-.','LineWidth',1,'Color',[0.8,0.5,1,0.8] );
plot(range, p_mae_all_GD_large(range),'.-','LineWidth',1,'Color',[1,0,0,0.8] );


xlabel('{Epoch}','interpreter','latex');
ylabel('{$\ell_1$ relative errors}','interpreter','latex');
ylim([0.06,0.12])
set(gca,'FontSize',8);
lgd = legend('GD: $\gamma_k=0.2,k\le 50$, $\gamma_k=0.1$ otherwise', ...
    'VI: $\gamma_k=0.4,k\le 50$, $\gamma_k=0.2$ otherwise', ...
    'GD: $\gamma_k=0.4,k\le 50$, $\gamma_k=0.2$ otherwise', ...
    'interpreter','latex','Location','best');
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,label,"_GDinstability_ProbPred",".pdf")); 

