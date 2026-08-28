% this is eliminate mu on graph with V = 5


%% Initialization
clearvars;
close all;
clc;

rng(2024, "twister");

%% Locate project directories relative to this script
scriptDir = string(fileparts(mfilename("fullpath")));

if strlength(scriptDir) == 0
    scriptDir = string(pwd);
end

% Add the main code directory and helper functions
addpath(scriptDir);

functionsDir = fullfile(scriptDir, "Functions");
if isfolder(functionsDir)
    addpath(functionsDir);
end

%% Input, output, and plot directories
thein   = fullfile(scriptDir, "Input")  + string(filesep);
theout  = fullfile(scriptDir, "Output") + string(filesep);
theplot = fullfile(scriptDir, "Plots")  + string(filesep);

% Create output directories if needed
if ~isfolder(theout)
    mkdir(theout);
end

if ~isfolder(theplot)
    mkdir(theplot);
end

%% Check graph parameter files
requiredFiles = ["peak.mat", "freq.mat"];

for filename = requiredFiles
    filepath = fullfile(scriptDir, filename);

    if ~isfile(filepath)
        error("Required graph parameter file is missing: %s", filepath);
    end
end

% Some graph functions load peak.mat and freq.mat from the current folder
originalDir = pwd;
restoreDir = onCleanup(@() cd(originalDir));
cd(scriptDir);

%% Experiment name
thedoc = "test_f3_graph_GLM";

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

%uhub

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

theta_psi = zeros(Nprime+N, Nprime, V^2); %initial value
X = zeros(Nprime+N, Nprime, V^2); %initial matrix



%% kernel recovery
use_GLMI = 1;
if use_GLMI
    label = "GLMI";
else
    label = "GLMS";
end




lr_schedule = [0.4*ones(50,1); 0.2*ones(50,1); 0.2*ones(50,1)];
bs_schedule = [800*ones(50,1); 800*ones(100,1)] ;

num_epoch = numel(lr_schedule );

nll_all = zeros(num_epoch,1);
p_mae_all = zeros(num_epoch,2);
mu_mae_all = zeros(num_epoch,1);

mu = zeros(V,1); % initial base
for ik = 1:V
    mu(ik) = mean(mean(event_data_tr(:,(ik-1)*(N+Nprime)+1:ik*(N+Nprime)),2))/h;
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
        idtr_batch = idtr_epoch((ibatch-1)*batch_size+1:ibatch*batch_size);
        edata = event_data_tr(idtr_batch,:); 


        
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


        id_type2 = (1:batch_size);



        % compute log-likelihood on type2
        yob = zeros(V,N*batch_size); % reshape of edata, boolean matrix size (V, N*batch_size)
        for k1 = 1:V
            yob(k1,eventJob{k1,1}) = 1;
        end
        ybarob = sum(yob,1);


        if strcmp(label, "GLMI")
            prob = Lambda;
            low = (prob<0);
            high = (prob>1);
            prob(low) = 1e-7;
            prob(high) = 1;
            prob_bar = min(1-1e-7,sum(prob,1));
        else
            prob = 1./(1+exp(-Lambda));
            prob_bar = min(1-1e-7,sum(prob,1));
        end

    
        lklhA1 = sum(yob.*log(prob),"all");
        lklhA2 = sum((1-ybarob).*log(1-prob_bar),"all");

        likelihood = (lklhA1+lklhA2);

        nll_sum = nll_sum-likelihood;
        
        


        %% compute grad field for F
  
        if use_GLMI % negative vector field to minimize the loss
            tmp = Lambda;
            low = (tmp<0);
            high = (tmp>1);
            tmp(low) = 0;
            tmp(high) = 1;
            dldLambda = -(tmp-yob); % (1,Nob*batchsize)
        else % negative vector field to minimize the loss
            dldLambda = -(1./(1+exp(-Lambda))-yob); % (1,Nob*batchsize) 
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
            for k2 = 1:V
               
                %% compute grad from dK
                Eobk1k2 = full(Eob{k1,k2});
                E_reshaped = reshape(Eobk1k2, [N+Nprime,N,batch_size]);
                
                dldK_k1k2 = sum(E_reshaped.*tmp_reshaped,3);
                
           
                
                %% clear dldK_k1k2_alongbatch;
                
                dldK_k1k2_extend = zeros(N+Nprime,N+Nprime);
                dldK_k1k2_extend(:,Nprime+1:N+Nprime) = dldK_k1k2;
                

                %% reshape dK to dPhi
                
                dldX_k1k2 = spdiags(dldK_k1k2_extend',diag_ind);
                dldX(:,:,k1+V*(k2-1)) = dldX_k1k2/batch_size; 
                
            end                
        end

        % gradient for mu from likelihood
        dldmu_list =  zeros([V,batch_size]);
        for ib = 1:batch_size
            dldmu_list(:,ib) = mean(dldLambda(:,N*(ib-1)+1:N*ib),2);
        end
        dldmu = mean(dldmu_list,2);
       
        



        

        

        
       
       
        
        %% update kernel 
        dX = dldX;  
        X = X + etaK*dX; %update X
        mu = mu+etaK*dldmu;



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

        figure(19),clf;
        subplot(121),
        plot( 1:iepoch-1, nll_all(1:iepoch-1), '.-');
        grid on; title('tr nll'); set(gca,'FontSize',15);
        subplot(122),
        plot( 1:iepoch-1, p_mae_all(1:iepoch-1), '.-');
        grid on; title('te L2 loss'); set(gca,'FontSize',15);
        drawnow();
    end
       
    
    

    % epoch average nll
    nll_epoch =nll_sum/ntr;
    nll_all(iepoch) = nll_epoch;

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

    if strcmp(label, 'GLMI')
        tmp = lam_val;
        low = (tmp<0);
        high = (tmp>1);
        tmp(low) = 1e-7;
        tmp(high) = 1;
        prob_te = tmp; 
        prob_bar_te = min(1-1e-7,sum(prob_te,1));
    elseif strcmp(label, 'GLMS')
        prob_te = 1./(1+exp(-lam_val)); 
        prob_bar_te = min(1-1e-7,sum(prob_te,1));
    end

    hatp_val = prob_te;
       
    

    





    

    mae_p_val = mean(abs( hatp_val - p_val_true),2);
    l1_p_val_true = mean( abs(p_val_true), 2);
    relmae_p_val = mae_p_val./l1_p_val_true;
    avg_p_mae = mean(relmae_p_val);
    worst_p_mae = max(relmae_p_val );

    p_mae_all(iepoch,:) = [avg_p_mae,worst_p_mae];
    toc
            

    %%
    fprintf('epoch %d, count of viol=%d, nll_tr=%6.4f; errp=%6.4f, %6.4f, errmu = %6.4f\n', ...
        iepoch, count_violates, nll_epoch, avg_p_mae, worst_p_mae, errmu);




   
end

%% kernel and mu
save(strcat(theout,thedoc,label,".mat"), "X");
save(strcat(theout,thedoc,label,"_mu.mat"), "mu");

return;




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

if strcmp(label, 'GLMI')
    tmp = Lambda_te;
    low = (tmp<0);
    high = (tmp>1);
    tmp(low) = 1e-7;
    tmp(high) = 1;
    Prob_te = tmp; 
elseif strcmp(label, 'GLMS')
    Prob_te = 1./(1+exp(-Lambda_te)); 
end

% Prob_te = (1-exp(-h.*Lambar_te))./Lambar_te.*Lambda_te;

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


%% error table
label = "GLMI"



proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
proberror = proberror.proberror;
mean_proberr = mean(proberror,1)
std_proberr = std(proberror,1)

label = "GLMS"


proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
proberror = proberror.proberror;
mean_proberr = mean(proberror,1)
std_proberr = std(proberror,1)



