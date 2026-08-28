%% Initialization
clearvars;
close all;
clc;

rng(2024, "twister");

%% Locate project directories relative to this script
scriptDir = string(fileparts(mfilename("fullpath")));

% Add functions stored beside this script
addpath(scriptDir);

% Also add the Functions subdirectory, if present
functionsDir = fullfile(scriptDir, "Functions");
if isfolder(functionsDir)
    addpath(functionsDir);
end

%% Input, output, and plot directories
thein   = fullfile(scriptDir, "Input")  + string(filesep);
theout  = fullfile(scriptDir, "Output") + string(filesep);
theplot = fullfile(scriptDir, "Plots")  + string(filesep);

% Create directories if necessary
requiredFolders = [thein, theout, theplot];

for folder = requiredFolders
    if ~isfolder(folder)
        mkdir(folder);
    end
end

%% Experiment name
thedoc = "test_f0_timeonly_stationary_HPE";

%% kenrel Psi
N = 32; %320;
Nprime = 16; %80;

T = 20; %max time length T
h = T/N;
tau_max = h*Nprime;

nt = (N+Nprime);
ntau = Nprime;
ut_kernel = 20; 
tau_grid = ((0.5:1:ntau)/nt)'*ut_kernel;

psi1 = 0.01*exp(-(8.*(tau_grid ).^2/25))...
    + 0.25*exp(-(6.*(tau_grid ).^2/25)).*sin(-0.05*pi+tau_grid/2*pi);

psi1= psi1';

figure(1),clf;
plot(tau_grid, psi1, '.-' )
grid on;

psi_full = ones(nt,1)*psi1;

%% make stationary
psi = psi_full;

% mask
a= ones(N+Nprime, Nprime);
psi_mask = triu(a, -N+1).*tril(a,0);
psi_mask = psi_mask(:,end:-1:1);
psi = psi.*psi_mask;


[im, jm, ~] = find(psi_mask);   %the nonzero indexing in the psi mask
indm = sub2ind( [N+Nprime, Nprime], im, jm);
vpsi = psi(indm); %[im,jm,vpsi] is the sparse repn of psi

figure(2),clf;
subplot(121), imagesc(psi);
title('Psi'); colorbar();
subplot(122), spy(psi_mask);
title('Psi mask'); 

% graph laplacian
Lt = make_1d_laplacian(Nprime+N)/h^2;
Ltau = make_1d_laplacian(Nprime)/h^2;



%% full K
% no mask
a= ones(N+Nprime, Nprime);
full_mask = triu(a, -N+1);
full_mask = full_mask(:,end:-1:1);
[imf, jmf, ~] = find(full_mask);
vpsi_full = psi_full(sub2ind( [N+Nprime, Nprime], imf, jmf));
K_full = sparse(imf,imf+jmf,vpsi_full, N+Nprime,N+Nprime);

%% K and Psi
jk = im+jm-Nprime;
K = sparse(im,jk,vpsi, N+Nprime,N);

% from K to psi
indk = sub2ind( [N+Nprime, N], im, jk);
vk = K(indk); %[im,jk,vk] is the sparse repn of psi

jm2 =jk+Nprime-im;
assert( norm(jm - jm2)<1e-15);

%psi2 = sparse(im,jk+Nprime-im, vk, N+Nprime, Nprime );
psi2 = sparse(im,jm, vk, N+Nprime, Nprime );
assert( norm(psi - psi2)<1e-15);

figure(3),clf;
subplot(121), imagesc(K);
title('K'); colorbar();
subplot(122), spy(K);



%% generate traj
phi_func = @(x) 1-exp(-x);

mu_true = 0.2; %0.125; 
    %mu = 0.2 for N= 320

M = 10000; 
y_ob = false(M, Nprime+N);

eta_ob = false(M, Nprime+N, N);
lambda_true = zeros(M,N);

disp('generating trajectories...')
tic,
for i = 1: Nprime+N
    
    if i <= Nprime
        
        if (0)
            y_ob(:, i ) = (rand( M, 1) < mu_true); % no use kernel
        else
            ypre = y_ob(:, 1:i-1);
            kernelt = K(1:i-1,i);
            lambdat =sum(bsxfun(@times, kernelt, ypre' ),1)'+ mu_true;
            pyt = phi_func( h*max(0,lambdat));
            y_ob(:, i) = (rand(M,1) <  pyt);
        end
    else
        t = i-Nprime;
        ypre = y_ob(:, t:t+Nprime-1);
        eta_ob(:, t:t+Nprime-1, t)= ypre;
        kernelt = K(t:t+Nprime-1,t);
        lambdat =sum(bsxfun(@times, kernelt, ypre' ),1)'+ mu_true;
        if min(lambdat) < 0
            warning( sprintf('min lamdba = %6.4f\n',min(lambdat)))
        end
       
        lambda_true(:,t) = lambdat; 
        pyt = phi_func( h*max(0,lambdat)); %todo: give warning if true lambda<0
        y_ob(:, i) = (rand(M,1) <  pyt);
    end
end
toc


%
[min_lam_true,i_min]= min(min(lambda_true,[],2));
max_lam_true= max(lambda_true(:));

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

%return;

%% split data
event_data = y_ob;

ntr = 4800; %32000; %16000;
nte = min(500,M -ntr);

tmp = randperm(M);
idx_tr = tmp(1:ntr);
idx_te = tmp(ntr+1:ntr+nte);


event_data_tr = event_data(idx_tr,:);
mu = mean(mean(event_data_tr,2))/h;
event_data_te = event_data(idx_te,:);
            
lam_true_val = lambda_true(idx_te,:);
p_val_true =  phi_func(h* lam_true_val);



%%
dT = ut_kernel/(N+Nprime);
event_data_tensor = event_data;
cont_event_data = cell(M,1);
tplt = dT:dT:ut_kernel;
for m = 1:M
      idx = logical(event_data_tensor(m,:));
      cont_event_data{m} = [cont_event_data{m},[tplt(idx);ones(size(tplt(idx)))]];
end

%% split data
event_data_tr = cont_event_data(idx_tr);
event_data_te = cont_event_data(idx_te);


%% training hyperparameters
num_epoch = 90;



batch_size = 400;
batch_size_schedule = batch_size*ones(num_epoch,1);
eta_schedule = 0.02*ones(num_epoch,1);

iepoch = 0;
likelihood_epoch = zeros(num_epoch,1);



V = 1;
AB = cell(3,1);
AB{1} = 1;
AB{2} = 1;
AB{3} = mu;

while iepoch < num_epoch  % epoch -> batch -> fold
    iepoch = iepoch+1;
    A = AB{1};
    beta = AB{2};
    fprintf('\n--- epoch %d ---\n', iepoch);

    batch_size = batch_size_schedule(iepoch);
    num_batch = floor(ntr/batch_size);
    idxepoch = randperm(ntr, num_batch*batch_size); %scramble the idx

    eta = eta_schedule(iepoch);

    likelihood_batch = [];
    for ibatch = 1:num_batch
        %% load batch training data
        idx_batch = idxepoch((ibatch-1)*batch_size+1:ibatch*batch_size);
        edata = event_data_tr(idx_batch); 
        
        % compute Lambda, size (1, Nob*batch_size)
        Lambda = cell(batch_size,1);
        for ib = 1:batch_size
            [~,num_events] = size(edata{ib});
            lam_now = edata{ib};
            lam_now(1,:) = 0;
            %%
            for ie = 1:num_events
        
                k1 = ceil(lam_now(2,ie));
                ti = edata{ib}(1,ie);
                for je = 1:(ie-1)
                    k2 = ceil(lam_now(2,je)); 
                    tj = edata{ib}(1,je);
                    if ti<=tj
                        continue
                    end
                    lam_now(1,ie) = lam_now(1,ie) + A(k2,k1)*beta*exp(-beta*(ti-tj));
                end
                lam_now(1,ie) = lam_now(1,ie) + mu(k1);
            end
            Lambda{ib} = lam_now;
        end

        % compute log-likelihood
        lklh1 = 0;
        lklh2 = 0;
        lklh3 = 0;
        for ib = 1:batch_size
            [~,num_events] = size(edata{ib});
            if num_events == 0
                continue
            end
            lklh1 = lklh1 + sum(log(Lambda{ib}(1,:)),'all');
            lhlh2 = lklh2 + sum(mu,'all')*ut_kernel;

            ui_mask = ceil(Lambda{ib}(2,:));
            A_u_ui = A(:,ui_mask);
            tmp = kron((1-exp(-beta*(ut_kernel-edata{ib}(1,:)))),ones(V,1));
            lklh3 = lklh3 + sum(A_u_ui.*tmp,'all');
        end

        likelihood = (lklh1-lklh2-lklh3)/batch_size;

        likelihood_batch = [likelihood_batch; likelihood];

        %% gradient from likelihood
        dldA = zeros(V,V);
        dldB = 0;
        dldmu = zeros(V,1);

        % compute gradient
        for ib = 1:batch_size
            %%
            [~,num_events] = size(edata{ib});
            if num_events == 0
                continue
            end
            lam_now = Lambda{ib};

            for ie = 1:num_events
                k1 = ceil(lam_now(2,ie)); 
                ti = edata{ib}(1,ie);
                for je = 1:(ie-1)
                    k2 = ceil(lam_now(2,je)); 
                    tj = edata{ib}(1,je);
                    if ti<=tj
                        continue
                    end
                    dldA(k2,k1) = dldA(k2,k1) + beta*exp(-beta*(ti-tj))/Lambda{ib}(1,ie);
                    dldB = dldB + A(k2,k1)*( exp(-beta*(ti-tj)) - beta*(ti-tj)*exp(-beta*(ti-tj)))/Lambda{ib}(1,ie);
                end
                dldmu(k1) = dldmu(k1)+1/Lambda{ib}(1,ie);
            end

            dldmu = dldmu-ut_kernel;


            for je = 1:num_events
                k2 = ceil(lam_now(2,je)); 
                tj = edata{ib}(1,je);
                dldA(k2,:) = dldA(k2,:) - (1-exp(-beta*(ut_kernel-tj)));

                
            end
            ui_mask = ceil(Lambda{ib}(2,:));
            A_u_ui = A(:,ui_mask);
            tmp = kron((beta*exp(-beta*(ut_kernel-edata{ib}(1,:)))),ones(V,1));
            dldB = dldB - sum(A_u_ui.*tmp,'all');

            
%                 dldA
%                 dldB
        end

        dldA = dldA./batch_size;
        dldB = dldB/batch_size;
        dldmu = dldmu/batch_size;




        %% update kernel matrix

        dA = dldA;

        dB = dldB;
            
        A = A + eta*dA;
        A = max(A,0);

        beta = beta + eta*dB; 

        mu(1:V-1) = mu(1:V-1)+eta*dldmu(1:V-1);
        mu = max(mu,1e-4);

        AB{1}=A;
        AB{2}=beta;
        AB{3}=mu;
        
        %%
        likelihood_ibatch = likelihood_batch(ibatch);
        fprintf('batch %d likelihood: %6.4e\n',ibatch,likelihood_ibatch)

    end
    %% training likelihood
    likelihood_epoch(iepoch) = mean(likelihood_batch);  
    if 1
        figure(19),clf;
        hold on;
        plot(1:iepoch, likelihood_epoch(1:iepoch),'x-');
        grid on;
        xlabel('epoch')
        title('tr likelihood');
        drawnow();
    end
end

nll_all = -likelihood_epoch;
save(strcat(theout,thedoc,"_nllall.mat"), "nll_all");
save(strcat(theout,thedoc,".mat"), "AB");

return;


%%
nll_all = load(strcat(theout,thedoc,"_nllall.mat"));
nll_all = nll_all.nll_all;
figure(19),clf;

box on;
plot( 1:iepoch, nll_all, '.-');
grid on;


xlabel('{Epoch}','interpreter','latex');
ylabel('{Negative log-likelihood}','interpreter','latex');
set(gca,'FontSize',8);
lgd = legend('train','Location','best');
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,"__TrainLogLike.pdf"))

%%
nt = (N+Nprime);
ntau = Nprime;
tau_grid = ((0.5:1:ntau)/nt)'*ut_kernel;

psi0 = zeros(nt, ntau, V^2);
for j = 1:ntau
    for k1 = 1:V
        for k2 = 1:V
            a = A(k1,k2);
            psi0(:,j,(k2-1)*V+k1) = a*exp(-beta.*(tau_grid(j)))*ones(nt,1);
        end
    end
end

downfactor = 1;
psi_learn = zeros(nt/downfactor,ntau/downfactor,V^2);
for k = 1:V^2
    psi_learn(:,:,k) = matrix_downsize(psi0(:,:,k), downfactor);
end
X = psi_learn;
imagesc(X)

%%
event_data_te = event_data(idx_te,:);
[Pob_te, Eob_te, eventJall_te,eventJob_te]= network_proj_mod_speed_N_graph2(event_data_te, X, N+Nprime, Nprime, 1);
% compute Lambda, size (1, N*batch_size)
Lambda_te = sum(Pob_te{1,1},1)+mu;

% reshape Lambda to pair the shape of data y
lambda_te = zeros(N*1, nte);
for b = 1:nte
    lambda_te(:,b) = reshape((Lambda_te(:,(b-1)*N+1:b*N))',[N*1,1]);
end
prob_te = (1-exp(-h*lambda_te));

[Pob_te_tru, Eob_te, eventJall_te,eventJob_te]= network_proj_mod_speed_N_graph2(event_data_te, psi, N+Nprime, Nprime, 1);
% compute Lambda, size (1, N*batch_size)
Lambda_te_tru = sum(Pob_te_tru{1,1},1)+mu_true;

% reshape Lambda to pair the shape of data y
lambda_te_tru = zeros(N*1, nte);
for b = 1:nte
    lambda_te_tru(:,b) = reshape((Lambda_te_tru(:,(b-1)*N+1:b*N))',[N*1,1]);
end
prob_te_tru = (1-exp(-h*lambda_te_tru));


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
std(proberror,1)
save(strcat(theout,thedoc,"_ProbPredErr",".mat"), "proberror");




