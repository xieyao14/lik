clear all; rng(2024);
addpath /Users/leviathaniety/Documents/MATLAB/Functions;
thein = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Input/";
theout = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Output/";
theplot = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Plots/";
thedoc = "test_f1_timeonly_highdim_GLM";

%% kenrel Psi
N = 32*10; %320;
Nprime = 8*10; %80;

T = 20; %max time length T
h = T/(N); 

downfactor  = 4;
nt = (N+Nprime)*downfactor;
ntau = Nprime*downfactor;
ut_kernel = 20;
t_grid = ((0.5:1:nt)/nt)'*ut_kernel;

psi0 = zeros(nt, ntau);
for j = 1:ntau
    psi0(:,j) = kernel_timeonly_highrank(t_grid, t_grid(j)); 
end
psi_full = matrix_downsize(psi0, downfactor);
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

%% full K
% no mask
a= ones(N+Nprime, Nprime);
full_mask = triu(a, -N+1);
full_mask = full_mask(:,end:-1:1);
[imf, jmf, ~] = find(full_mask);
vpsi_full = psi_full(sub2ind( [N+Nprime, Nprime], imf, jmf));
K_full = sparse(imf,imf+jmf,vpsi_full, N+Nprime,N+Nprime);

%% K and Psi
% [i,j,v]= find(psi);
% K = sparse(i,i+j,v, N+Nprime,N+Nprime);
% K = K(:,Nprime+1:end);
% K = from_psi_to_K(psi, N, Nprime);

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

%% check kernel rank
psi_sub = psi(Nprime: N, :);
[~,s1,v1]= svd(psi_sub, "econ");
s_true = diag(s1);

svd_thres = .2; %target rank depends on the accuracy of recovery
                %with finite sample, the recovered psi before truncation is like psi plus noise   

%%

%% generate traj
phi_func = @(x) 1-exp(-x);

mu_true = 0.2; %0.125; 
    %mu = 0.2 for N= 320

M = 40000; %40000;
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



%% split data
event_data = y_ob;

ntr = 16000; %32000; %16000;
nte = min(500,M -ntr);

tmp = randperm(M);
idx_tr = tmp(1:ntr);
idx_te = tmp(ntr+1:ntr+nte);


event_data_tr = event_data(idx_tr,:);
event_data_te = event_data(idx_te,:);
            
lam_true_val = lambda_true(idx_te,:);
p_val_true =  phi_func(h* lam_true_val);


%% 

% zero init
theta_psi = zeros(Nprime+N, Nprime); %initial value
X = zeros(Nprime+N, Nprime); %initial matrix

% rank 1 init
%theta_psi = psi*v1(:,1)*v1(:,1)';

%% kernel recovery
use_GLMI = 0;
if use_GLMI
    label = "GLMI";
else
    label = "GLMS";
end
warm_start = 0;

mode = "eliminate";



lr_schedule = [0.4*ones(100,1), 0.2*ones(100,1), 0.2*ones(100,1) ];
bs_schedule = [400*ones(100,1), 400*ones(100,1), 400*ones(100,1) ];

num_epoch = numel(lr_schedule );

nll_all = zeros(num_epoch,1);
err1_kernel_all = zeros(num_epoch,1);
p_mae_all = zeros(num_epoch,2);
mu_mae_all = zeros(num_epoch,1);

if mode == "cheat"
    mu = mu_true;
elseif mode == "eliminate"
    mu = mean(mean(event_data_tr,2))/h;
end

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
        %%
        edata = event_data_tr(idtr_batch,:); 
             
        % compute projected data
        [Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N_graph2(edata, X, N+Nprime, Nprime, 1);
        
        % compute Lambda, size (1, Nob*batch_size)
        Lambda = sum(Pob{1,1},1)+mu;

        % decide if Lambda violates buffer
        lam_batch = zeros(batch_size, N*1); %[batch_size, N]
        for b = 1:batch_size
            lam_batch(b,:) = reshape((Lambda(:,(b-1)*N+1:b*N))',[N*1,1]);
        end

        % compute log-likelihood
        yob = zeros(1,N*batch_size); % reshape of edata, boolean matrix size (1, N*batch_size)
        yob(1,eventJob{1,1}) = 1;
        



        if strcmp(label, "GLMI")
            prob = Lambda;
            low = (prob<0);
            high = (prob>1);
            prob(low) = 1e-7;
            prob(high) = 1;
        else
            prob = 1./(1+exp(-Lambda));
        end

    
        lklhA1 = sum(yob.*log(prob),"all");
        lklhA2 = sum((1-yob).*log(1-prob),"all");

        likelihood = (lklhA1+lklhA2);

        nll_sum = nll_sum-likelihood;


        %% compute grad field for F
        dldX = zeros(N+Nprime,Nprime); % gardient on kernel matrix 
        
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


%         dldlam_batch = zeros(batch_sz2, N*1); %[batch_size, N]
%         for b = 1:batch_sz2
%             dldlam_batch(b,:) = reshape((dldLambda(:,(b-1)*N+1:b*N))',[N*1,1]);
%         end

        id_type2 = (1:batch_size);

        diag_ind = -(1:Nprime);

        % compute grad from dK
        Eob11 = Eob{1,1};
        dldK_k1k2 = zeros(N+Nprime,N);
        for ib = id_type2
            % sum over batch type 2
            dldK_k1k2 = dldK_k1k2 + Eob11(:,N*(ib-1)+1:N*ib).*kron(dldLambda(N*(ib-1)+1:N*ib), ones(N+Nprime,1));
        end
        % clear dldK_k1k2_alongbatch;
        dldK_k1k2_extend = zeros(N+Nprime,N+Nprime);
        dldK_k1k2_extend(:,Nprime+1:N+Nprime) = dldK_k1k2;

        % reshape dK to dPsi
        dldX_k1k2 = spdiags(dldK_k1k2_extend',diag_ind);
        dldX(:,:) = dldX_k1k2/batch_size; 
        
        % gradient for mu from likelihood
        dldmu_list =  zeros([1,batch_size]);
        for ib = 1:batch_size
            dldmu_list(ib) = mean(dldLambda(1,N*(ib-1)+1:N*ib),2);
        end
        dldmu = mean(dldmu_list,2);



        
        
       
       
        
        %% update kernel 
        dX = dldX;  
        X = X + etaK*dX; %update X
        mu = mu+etaK*dldmu;

        

        



    end

    
    theta_psi = X;

    if mod(iepoch,10)==0
        figure(9),clf;
        imagesc(theta_psi); colorbar();
        title(sprintf('epoch %d',iepoch));
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
       
    %kernel error on submatrix    
    theta_psi_sub = theta_psi(Nprime: N, :);
    err1 = norm(theta_psi_sub - psi_sub, 'fro')/ norm(psi_sub,'fro');
    err1_kernel_all(iepoch) = err1;
    

    % epoch average nll
    nll_epoch =nll_sum/ntr;
    nll_all(iepoch) = nll_epoch;

    % mu error
    errmu = abs(mu-mu_true)/mu_true;
    mu_mae_all(iepoch) = errmu;

    toc,

    %% prediction error on test traj
    disp('testing...')
    tic,
    [Pob_te, Eob_te, eventJall_te,eventJob_te]= network_proj_mod_speed_N_graph2(event_data_te, X, N+Nprime, Nprime, 1);

    % compute Lambda, size (1, N*batch_size)
    Lambda_te = sum(Pob_te{1,1},1)+mu;
    
    % reshape Lambda to pair the shape of data y
    lam_val = zeros(nte,N*1);
    for b = 1:nte
        lam_val(b,:) = reshape((Lambda_te(:,(b-1)*N+1:b*N))',[N*1,1]);
    end
    if strcmp(label, "GLMI")
        prob = lam_val;
        low = (prob<0);
        high = (prob>1);
        prob(low) = 1e-7;
        prob(high) = 1;
    else
        prob = 1./(1+exp(-lam_val));
    end
    hatp_val = prob;







    
    
    

    mae_p_val = mean(abs( hatp_val - p_val_true),2);
    l1_p_val_true = mean( abs(p_val_true), 2);
    relmae_p_val = mae_p_val./l1_p_val_true;
    avg_p_mae = mean(relmae_p_val);
    worst_p_mae = max(relmae_p_val );

    p_mae_all(iepoch,:) = [avg_p_mae,worst_p_mae];
    toc
            

    %%
    fprintf('epoch %d, count of viol=%d, nll_tr=%6.4f; err1=%6.4f, errmu = %6.4f, errp=%6.4f, %6.4f\n', ...
        iepoch, count_violates, nll_epoch, err1, errmu, avg_p_mae, worst_p_mae );



   
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
[Pob_te, Eob_te, eventJall_te,eventJob_te]= network_proj_mod_speed_N_graph2(event_data_te, X, N+Nprime, Nprime, 1);
% compute Lambda, size (1, N*batch_size)
Lambda_te = sum(Pob_te{1,1},1)+mu;

% reshape Lambda to pair the shape of data y
lambda_te = zeros(N*1, nte);
for b = 1:nte
    lambda_te(:,b) = reshape((Lambda_te(:,(b-1)*N+1:b*N))',[N*1,1]);
end
if strcmp(label, "GLMI")
    prob = lambda_te;
    low = (prob<0);
    high = (prob>1);
    prob(low) = 1e-7;
    prob(high) = 1;
else
    prob = 1./(1+exp(-lambda_te));
end
prob_te = prob;


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