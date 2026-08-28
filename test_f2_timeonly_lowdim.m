clear all; rng(2024);
addpath /Users/leviathaniety/Documents/MATLAB/Functions;
thein = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Input/";
theout = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Output/";
theplot = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Plots/";
thedoc = "test_f2_timeonly_lowhdim";

%% kernel Psi
N = 32*1; %320;
Nprime = 8*1; %80;

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
use_VI = 0;

% barrier hyper parameter
min_b = 0.01; %0.03;
B_weight = 0.1; % 0.1
smooth_weight = 0.2; % 1*h^2
mode = "eliminate";


%lr_schedule = [0.4*ones(50,1), 0.2*ones(50,1), 0.1*ones(50,1) ];
if use_VI
    label = "VI";
    lr_schedule = [0.4*ones(100,1), 0.2*ones(100,1), 0.2*ones(100,1) ]/10;
else
    label = "GD";
    lr_schedule = [0.2*ones(100,1), 0.1*ones(100,1), 0.1*ones(100,1) ]/10;
end

bs_schedule = [400*ones(100,1), 400*ones(100,1), 400*ones(100,1) ];

num_epoch = numel(lr_schedule );

nll_all = zeros(num_epoch,1);
err1_kernel_all = zeros(num_epoch,1);
p_mae_all = zeros(num_epoch,2);
mu_mae_all = zeros(num_epoch,1);
mu_all = zeros(num_epoch,1);

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
        id_type1_squeeze = zeros([1,N*batch_sz1]);
        for bb = 1:batch_sz1
            id = id_type1(bb);
            id_type1_squeeze((bb-1)*N+1:bb*N) = (id-1)*N+1:id*N;
        end
        batch_sz2 = batch_size-batch_sz1;
        
        count_violates = count_violates + batch_sz1;


        % compute log-likelihood on type2
        yob = sparse(1,eventJob{1,1},1,1,N*batch_size);
        yob_type1 = yob(id_type1_squeeze);
        yob_type2 = yob;
        yob_type2(id_type1_squeeze) = [];

        Lambda_type1 = Lambda(id_type1_squeeze);
        Lambda_type2 = Lambda;
        Lambda_type2(id_type1_squeeze) = [];

        lklhA1 = sum(yob_type2.*log( (1-exp(-h*Lambda_type2))./Lambda_type2) ,"all");
        lklhA2 = sum(yob_type2.*log(Lambda_type2),"all");
        lklhB = sum((1-yob_type2).*Lambda_type2,"all")*h;
        likelihood = (lklhA1+lklhA2-lklhB);

        nll_sum = nll_sum-likelihood;
        count_type2 = count_type2 +batch_sz2; %number of trajectories that is added into nll


        %% compute grad field for F on type2
        dldX = zeros(N+Nprime,Nprime); % gardient on kernel matrix 
        if use_VI % negative vector field to minimize the loss
            dldLambda =- kron((1-(exp(-h*Lambda))-yob), ones(1,1)); % (1,Nob*batch_sz2) 
        else % negative vector field to minimize the loss
            dldLambda = yob.*(h./(exp(h*Lambda)-1)) - h*(1-yob); % (1,Nob*batch_sz2)
        end


%         dldlam_batch = zeros(batch_sz2, N*1); %[batch_size, N]
%         for b = 1:batch_sz2
%             dldlam_batch(b,:) = reshape((dldLambda(:,(b-1)*N+1:b*N))',[N*1,1]);
%         end



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
        dldX(:,:) = dldX_k1k2/batch_sz2; 
        



        %% compute grad field for barrier on type1
        dBdX = zeros(N+Nprime,Nprime);
        if batch_sz1 > 0

            if 0
                dBdLambda = -min_b./(Lambda);
            else
                dBdLambda = (Lambda-min_b)/(min_b*0.1);
            end
            
%             dBdLambda = -min_b./(Lambda);
            dBdLambda = dBdLambda.*(Lambda < min_b);
    
            % compute grad from dK
            dBdK_k1k2 = zeros(N+Nprime,N);
            for ib = id_type1
                % sum over batch
                dBdK_k1k2 = dBdK_k1k2 + Eob11(:,N*(ib-1)+1:N*ib).*kron(dBdLambda(N*(ib-1)+1:N*ib), ones(N+Nprime,1));
            end
            % clear dBdK_k1k2_alongbatch
            dBdK_k1k2_extend = zeros(N+Nprime,N+Nprime);
            dBdK_k1k2_extend(:,Nprime+1:Nprime+N) = dBdK_k1k2;
    
            % reshape dK to dPsi
            dBdX_k1k2 = spdiags(dBdK_k1k2_extend',diag_ind);
            dBdX(:,:) = dBdX_k1k2; 
            dBdX = dBdX*B_weight; %remove last row of dX
        end


        
       
       
        
        %% update kernel 
        dX = dldX-dBdX;  
        X = X + etaK*dX; %update X

        

        % tSVD to enforce low-rank: only apply after burn in
        if (0)
            tmp = X;
            [~,s2,v2]= svd(tmp(Nprime: N, :), "econ");
            s2 = diag(s2);
            if s2(1) < svd_thres %all less than thres
                rkpsi = Nprime;
            else
                rkpsi = sum(s2 > svd_thres);
                assert( rkpsi > 0);
                if rkpsi < Nprime
                    v2proj = v2(:,1:rkpsi);
                    theta_psi = (theta_psi*v2proj)*v2proj';
                end
            end
        end

        %% update mu under eliminate mode
        if mode == "eliminate"
            mu_proposed = search_mu_timeonly(yob, sum(Pob{1,1},1), h, batch_size, mu);
            mu = 0.9*mu+0.1*mu_proposed;
            mu_all(iepoch) = mu;
        end

    end

    % smoothify kernel
    Lt = make_1d_laplacian(Nprime+N);
    Ltau = make_1d_laplacian(Nprime);

    X = X- etaK*smooth_weight*(Lt*X);

    X = X';
    X = X- etaK*smooth_weight*(Ltau*X);
    X = X';
    theta_psi = X;

    if mod(iepoch,10)==0
        figure(9),clf;
        imagesc(theta_psi); colorbar();
        title(sprintf('epoch %d',iepoch));
        drawnow();
    end
       
    %kernel error on submatrix    
        % todo: compute the vector on active region, N'*N size
    theta_psi_sub = theta_psi(Nprime: N, :);
    err1 = norm(theta_psi_sub - psi_sub, 'fro')/ norm(psi_sub,'fro');
    err1_kernel_all(iepoch) = err1;
    

    % epoch average nll
    nll_epoch =nll_sum/count_type2;
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
    hatp_val = (1-exp(-h*lam_val));







    min_lam_val = min(lam_val,[],2);
    if sum( min_lam_val < min_b ) > 0
        warning('has violation on val set');
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
            

    %%
    fprintf('epoch %d, count of viol=%d, nll_tr=%6.4f; err1=%6.4f, errmu = %6.4f, errp=%6.4f, %6.4f\n', ...
        iepoch, count_violates, nll_epoch, err1, errmu, avg_p_mae, worst_p_mae );



   
end

%% kernel
save(strcat(theout,thedoc,label,".mat"), "X");

%% kernel plot
% label = "GD"
% 
% X = load(strcat(theout,thedoc,label,".mat"));
% X = X.X;

figure(12),clf;
fig = imagesc(reshape(X(:,1:Nprime,1),[N+Nprime,Nprime]), 'YData', [-Nprime, N]); 
line([0,Nprime]+0.5, [0,-Nprime], 'Color', 'w','LineWidth',1);
line([0,Nprime]+0.5, [N,N-Nprime], 'Color', 'w','LineWidth',1);
colorbar();
fig.AlphaData = 0.8;
set(gca,'FontSize',8)
ax = gca;
xlabel('{Time delay} $j-i$','fontsize',8,'interpreter','latex');
ylabel('{Time index} $i$','fontsize',8,'interpreter','latex');
drawnow();
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 157 180]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,"_EstKer",label,".pdf"));  

%% kernel error
tmpX = X(Nprime+1:N,1:Nprime);
tmppsi = psi(Nprime+1:N,1:Nprime);
tmpdiff = tmpX-tmppsi;
if 1
    tmppsi = reshape(tmppsi, [(N-Nprime)*Nprime,1]);
    tmpdiff = reshape(tmpdiff, [(N-Nprime)*Nprime,1]);
end
kernel_error = [norm(tmpdiff,1);norm(tmpdiff,2);norm(tmpdiff,"inf")];
kernel_rela_err = kernel_error./[norm(tmppsi,1);norm(tmppsi,2);norm(tmppsi,"inf")];
fprintf(['-----Kernel Error Table----- \n' ...
    'Line 1: Rel l1, Rel l2, Rel l inf\n'])
kernel_rela_err'
save(strcat(theout,thedoc,label,"_EstKerRelErr.mat"), "kernel_rela_err");

%% mu
save(strcat(theout,thedoc,label,"_mu.mat"), "mu");

figure(7),clf;

box on;
plot(1:iepoch, mu_all(1:iepoch),'.-','LineWidth',1);
grid on;

xlabel('{Epoch}','interpreter','latex');
ylabel('{$\mu$}','interpreter','latex');
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
save(strcat(theout,thedoc,label,"_ProbPredErr",".mat"), "proberror");


return;

%% true kernel plot
figure(2),clf;
fig = imagesc(reshape(psi(:,1:Nprime,1),[N+Nprime,Nprime]), 'YData', [-Nprime, N]); 
line([0,Nprime]+0.5, [0,-Nprime], 'Color', 'w','LineWidth',1);
line([0,Nprime]+0.5, [N,N-Nprime], 'Color', 'w','LineWidth',1);
colorbar();
fig.AlphaData = 0.8;
set(gca,'FontSize',8)
ax = gca;
xlabel('{Time delay} $j-i$','fontsize',8,'interpreter','latex');
ylabel('{Time index} $i$','fontsize',8,'interpreter','latex');
drawnow();
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 157 180]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,"_TrueKer.pdf"));  

%% error table
label = "GD"

X = load(strcat(theout,thedoc,label,".mat"));
X = X.X;
mu = load(strcat(theout,thedoc,label,"_mu.mat"));
mu = mu.mu;

kernel_rela_err = load(strcat(theout,thedoc,label,"_EstKerRelErr.mat"));
kernel_rela_err = kernel_rela_err.kernel_rela_err

mu_rela_err = load(strcat(theout,thedoc,label,"_EstMuRelErr.mat"));
mu_rela_err = mu_rela_err.mu_rela_err

proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
proberror = proberror.proberror;
mean_proberr = mean(proberror,1)
std_proberr = std(proberror,1)

label = "VI"

X = load(strcat(theout,thedoc,label,".mat"));
X = X.X;
mu = load(strcat(theout,thedoc,label,"_mu.mat"));
mu = mu.mu;

kernel_rela_err = load(strcat(theout,thedoc,label,"_EstKerRelErr.mat"));
kernel_rela_err = kernel_rela_err.kernel_rela_err

mu_rela_err = load(strcat(theout,thedoc,label,"_EstMuRelErr.mat"));
mu_rela_err = mu_rela_err.mu_rela_err

proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
proberror = proberror.proberror;
mean_proberr = mean(proberror,1)
std_proberr = std(proberror,1)


%% prob prediction from VI and GD
tmp1 = isfile(strcat(theout,thedoc,"GD",".mat"));
tmp2 = isfile(strcat(theout,thedoc,"VI",".mat"));


if tmp1 && tmp2
    idx_te= randperm(nte, nte); 
    batch_size0 = 4;
    tmp = idx_te(1:batch_size0);
    tmp = [7         102    181    480 ]
    
    % LIKGD
    X = load(strcat(theout,thedoc,"GD",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"GD","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall,eventJob]= network_proj_mod_speed_N_graph2(event_data_te(tmp,:), X, N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_GD = (1-exp(-h.*teLambda));

    % LIKVI
    X = load(strcat(theout,thedoc,"VI",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"VI","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall,eventJob]= network_proj_mod_speed_N_graph2(event_data_te(tmp,:), X, N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_VI = (1-exp(-h.*teLambda));


    [truP, truE, eventJall,eventJob]= network_proj_mod_speed_N_graph2(event_data_te(tmp,:), psi, N+Nprime, Nprime, 1);
    truLambda = sum(truP{1,1},1)+mu_true;
    truprob = (1-exp(-h.*truLambda));

    figure(13),clf;
    ax = tiledlayout(floor(batch_size0/2),2);

    for b = 1:batch_size0
        nexttile
        hold on;    
        grid on;
        box on;
        jminvis = N*(b-1)+1;
        jmaxvis = N*b;
        eventJall1 = eventJob{1,1};
        eventJallvis = eventJall1(logical((eventJall1 < jmaxvis).*(eventJall1>= jminvis)));

        teLambdavis = teprob_GD(1,jminvis:jmaxvis);
        fig = plot(teLambdavis,'-.', 'linewidth', 1.5, 'Color', [1,0,0,0.5]);

        teLambdavis = teprob_VI(1,jminvis:jmaxvis);
        plot(teLambdavis,'--', 'linewidth', 1.5, 'Color', [0,0,1,0.5]);
    
        truLambdavis = truprob(1,jminvis:jmaxvis);
        plot(truLambdavis,'-', 'linewidth', 1.5, 'Color', [0 0 0 0.3]);
        plot(eventJallvis-jminvis+1,truprob(eventJallvis),'+');

        xlabel('{Time index} $j$','interpreter','latex');
        ylabel('{Predicted probability}','interpreter','latex');

%         ylim([0,0.07]);

        xticks(0:8:N);
        xticklabels(0:8:N);
        lgd = legend('TULIK-GD','TULIK-VI','true','Location','best');
        lgd.NumColumns = 1;
        set(gca,'FontSize',8);

    end


%     figH = gcf;
%     set(figH, 'Units', 'points','OuterPosition', [0 0 470 320]);
%     exportgraphics(ax,strcat(theplot,thedoc,"_ProbPred",".pdf"));
end


%% min lambda on training/testing data
label = "GD"

X = load(strcat(theout,thedoc,label,".mat"));
X = X.X;
mu = load(strcat(theout,thedoc,label,"_mu.mat"));
mu = mu.mu;

abs(mu-mu_true)/mu


tmpX = X(Nprime+1:N,1:Nprime);
tmppsi = psi(Nprime+1:N,1:Nprime);
tmpdiff = tmpX-tmppsi;
tmppsi = reshape(tmppsi,[(N-Nprime)*Nprime,1]);
tmpdiff = reshape(tmpdiff,[(N-Nprime)*Nprime,1]);
% if 1
%     tmppsi = reshape(tmppsi, [(N-Nprime)*Nprime,1]);
%     tmpdiff = reshape(tmpdiff, [(N-Nprime)*Nprime,1]);
% end
kernel_error = [norm(tmpdiff,1);norm(tmpdiff,2);norm(tmpdiff,"inf")];
kernel_rela_err = kernel_error./[norm(tmppsi,1);norm(tmppsi,2);norm(tmppsi,"inf")];
fprintf(['-----Kernel Error Table----- \n' ...
    'Line 1: Rel l1, Rel l2, Rel l inf\n'])
kernel_rela_err'


%%
[Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N_graph2(event_data_tr, X, N+Nprime, Nprime, 1);
        
Lambda = sum(Pob{1,1},1)+mu;

% decide if Lambda violates buffer
lam_batch = zeros(ntr, N*1); %[batch_size, N]
for b = 1:ntr
    lam_batch(b,:) = reshape((Lambda(:,(b-1)*N+1:b*N))',[N*1,1]);
end

min_lam_tr = min( lam_batch, [],2);
tmp_tr = sort(min_lam_tr);
fprintf("minimum training lambda\n")
tmp_tr(1:5)'

[Pob, Eob, eventJall,eventJob]= network_proj_mod_speed_N_graph2(event_data_te, X, N+Nprime, Nprime, 1);
        
Lambda = sum(Pob{1,1},1)+mu;

% decide if Lambda violates buffer
lam_batch = zeros(nte, N*1); %[batch_size, N]
for b = 1:nte
    lam_batch(b,:) = reshape((Lambda(:,(b-1)*N+1:b*N))',[N*1,1]);
end

min_lam_te = min( lam_batch, [],2);
tmp_te = sort(min_lam_te);
fprintf("minimum testing lambda\n")
tmp_te(1:5)'



%% prob prediction from VI and GD
tmp1 = isfile(strcat(theout,thedoc,"GD",".mat"));
tmp2 = isfile(strcat(theout,thedoc,"VI",".mat"));


if tmp1 && tmp2
    idx_te= randperm(nte, nte); 
    batch_size0 = 4;
    tmp = idx_te(1:batch_size0);
    tmp = [ 175 181  465  480 ]
    
    % LIKGD
    X = load(strcat(theout,thedoc,"GD",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"GD","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall]= network_proj(event_data_te(tmp,:), X(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_GD = (1-exp(-h.*teLambda));

    % LIKVI
    X = load(strcat(theout,thedoc,"VI",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"VI","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall]= network_proj(event_data_te(tmp,:), X(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_VI = (1-exp(-h.*teLambda));


    [truP, truE, eventJall]= network_proj(event_data_te(tmp,:), psi(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    truLambda = sum(truP{1,1},1)+mu_true;
    truprob = (1-exp(-h.*truLambda));

    figure(13),clf;
    ax = tiledlayout(floor(batch_size0/2),2);

    for b = 1:batch_size0
        nexttile
        hold on;    
        grid on;
        box on;
        jminvis = (N+Nprime)*(b-1)+1;
        jmaxvis = (N+Nprime)*b;
        eventJall1 = eventJall{1,1};
        eventJallvis = eventJall1(logical((eventJall1 < jmaxvis).*(eventJall1>= jminvis)));

        teLambdavis = teprob_GD(1,jminvis:jmaxvis);
        fig = plot(teLambdavis,'-.', 'linewidth', 1.5, 'Color', [1,0,0,0.5]);

        teLambdavis = teprob_VI(1,jminvis:jmaxvis);
        plot(teLambdavis,'--', 'linewidth', 1.5, 'Color', [0,0,1,0.5]);
    
        truLambdavis = truprob(1,jminvis:jmaxvis);
        plot(truLambdavis,'-', 'linewidth', 1.5, 'Color', [0 0 0 0.3]);
        plot(eventJallvis-jminvis+1,truprob(eventJallvis),'+');

        xlabel('{Time index} $j$','interpreter','latex');
        ylabel('{Predicted probability}','interpreter','latex');

        ylim([0,0.6]);
        xlim([1,N+Nprime]);

        xticks(0:8:N+Nprime);
        xticklabels(-Nprime:8:N);
        lgd = legend('TULIK-GD','TULIK-VI','true','Location','best');
        lgd.NumColumns = 1;
        set(gca,'FontSize',8);

    end


    figH = gcf;
    set(figH, 'Units', 'points','OuterPosition', [0 0 470 320]);
    exportgraphics(ax,strcat(theplot,thedoc,"_ProbPred",".pdf"));
end


%% adjusted prob prediction from VI and GD
tmp1 = isfile(strcat(theout,thedoc,"GD",".mat"));
tmp2 = isfile(strcat(theout,thedoc,"VI",".mat"));


if tmp1 && tmp2
    idx_te= randperm(nte, nte); 
    batch_size0 = 3;
    tmp = idx_te(1:batch_size0);
    tmp = [ 181  465  480 ]
    
    % LIKGD
    X = load(strcat(theout,thedoc,"GD",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"GD","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall, eventJob]= network_proj_mod(event_data_te(tmp,:), X(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_GD = (1-exp(-h.*teLambda));

    % LIKVI
    X = load(strcat(theout,thedoc,"VI",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"VI","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall,eventJob]= network_proj_mod(event_data_te(tmp,:), X(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_VI = (1-exp(-h.*teLambda));


    [truP, truE, eventJall, eventJob]= network_proj_mod(event_data_te(tmp,:), psi(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    truLambda = sum(truP{1,1},1)+mu_true;
    truprob = (1-exp(-h.*truLambda));

    figure(13),clf;
    ax = tiledlayout(1,3);

    for b = 1:batch_size0
        nexttile
        hold on;    
        grid on;
        box on;
        jminvis = (N)*(b-1)+1;
        jmaxvis = (N)*b;
        eventJall1 = eventJob{1,1};
        eventJallvis = eventJall1(logical((eventJall1 < jmaxvis).*(eventJall1>= jminvis)));

        teLambdavis = teprob_GD(1,jminvis:jmaxvis);
        fig = plot(teLambdavis,'-.', 'linewidth', 1.5, 'Color', [1,0,0,0.5]);

        teLambdavis = teprob_VI(1,jminvis:jmaxvis);
        plot(teLambdavis,'--', 'linewidth', 1.5, 'Color', [0,0,1,0.5]);
    
        truLambdavis = truprob(1,jminvis:jmaxvis);
        plot(truLambdavis,'-', 'linewidth', 1.5, 'Color', [0 0 0 0.3]);
        plot(eventJallvis-jminvis+1,truprob(eventJallvis),'+');

        xlabel('{Time index} $j$','interpreter','latex');
        ylabel('{Predicted probability}','interpreter','latex');

        ylim([0,0.6]);
        xlim([1,N]);

        xticks(0:8:N);
        xticklabels(0:8:N);
        if b == 1
            lgd = legend('TULIK-GD','TULIK-VI','true','Location','best');
            lgd.NumColumns = 1;
        end
        set(gca,'FontSize',8);

    end


    figH = gcf;
    set(figH, 'Units', 'points','OuterPosition', [0 0 470 200]);
    exportgraphics(ax,strcat(theplot,thedoc,"_ProbPred",".pdf"));
end

%% adjusted prob prediction from VI and GD
tmp1 = isfile(strcat(theout,thedoc,"GD",".mat"));
tmp2 = isfile(strcat(theout,thedoc,"VI",".mat"));


if tmp1 && tmp2
    idx_te= randperm(nte, nte); 
    batch_size0 = 3;
    tmp = idx_te(1:batch_size0);
    tmp = [ 181  465  480 ]
    
    % LIKGD
    X = load(strcat(theout,thedoc,"GD",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"GD","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall, eventJob]= network_proj_mod(event_data_te(tmp,:), X(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_GD = (1-exp(-h.*teLambda));

    % LIKVI
    X = load(strcat(theout,thedoc,"VI",".mat"));
    X = X.X;
    mu = load(strcat(theout,thedoc,"VI","_mu.mat"));
    mu = mu.mu;

    [teP, teE, eventJall,eventJob]= network_proj_mod(event_data_te(tmp,:), X(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    teLambda = sum(teP{1,1},1)+mu;
    teprob_VI = (1-exp(-h.*teLambda));


    [truP, truE, eventJall, eventJob]= network_proj_mod(event_data_te(tmp,:), psi(1:N+Nprime-1,:), N+Nprime, Nprime, 1);
    truLambda = sum(truP{1,1},1)+mu_true;
    truprob = (1-exp(-h.*truLambda));

    figure(13),clf;
    ax = tiledlayout(1,3);

    for b = 1:batch_size0
        nexttile
        hold on;    
        grid on;
        box on;
        jminvis = (N)*(b-1)+1;
        jmaxvis = (N)*b;
        eventJall1 = eventJob{1,1};
        eventJallvis = eventJall1(logical((eventJall1 < jmaxvis).*(eventJall1>= jminvis)));

        teLambdavis = teprob_GD(1,jminvis:jmaxvis);
        fig = plot(teLambdavis,'-.', 'linewidth', 1.5, 'Color', [1,0,0,0.5]);

        teLambdavis = teprob_VI(1,jminvis:jmaxvis);
        plot(teLambdavis,'--', 'linewidth', 1.5, 'Color', [0,0,1,0.5]);
    
        truLambdavis = truprob(1,jminvis:jmaxvis);
        plot(truLambdavis,'-', 'linewidth', 1.5, 'Color', [0 0 0 0.3]);

        if length(eventJallvis)>0
            plot(eventJallvis-jminvis+1, 0*ones(1,length(eventJallvis)),'.', 'MarkerSize',10, 'Color', 'red');
        end
%         plot(eventJallvis-jminvis+1,truprob(eventJallvis),'+');

        xlabel('{Time index} $j$','interpreter','latex');
        ylabel('{Predicted probability}','interpreter','latex');

        ylim([0,0.6]);
        xlim([1,N]);

        xticks(0:8:N);
        xticklabels(0:8:N);
        if b == 1
            lgd = legend('TULIK-GD','TULIK-VI','true','events','Location','best');
            lgd.NumColumns = 1;
        end
        set(gca,'FontSize',8);

    end


    figH = gcf;
    set(figH, 'Units', 'points','OuterPosition', [0 0 470 200]);
    exportgraphics(ax,strcat(theplot,thedoc,"_ProbPred_Adjusted",".pdf"));
end

%% 


nll_all = load(strcat(theout,thedoc,"GD","_TrainLogLike",".mat"));
nll_all_GD = nll_all.nll_all;
nll_all = load(strcat(theout,thedoc,"VI","_TrainLogLike",".mat"));
nll_all_VI = nll_all.nll_all;

figure(8),clf;

box on;
plot(1:num_epoch, nll_all_GD,'-.', 'linewidth', 1, 'Color', [1,0,0,0.5]);
grid on;
hold on;
plot(1:num_epoch, nll_all_VI,'--', 'linewidth', 1, 'Color', [0,0,1,0.5]);

xlabel('{Epoch}','interpreter','latex');
ylabel('{Negative log-likelihood}','interpreter','latex');
set(gca,'FontSize',8);
lgd = legend('TULIK-GD','TULIK-VI','Location','best');
lgd.NumColumns = 1;
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,"_TrainLogLike",".pdf")); 


%%
mu_all = load(strcat(theout,thedoc,"GD","_mu_all.mat"));
mu_all_GD = mu_all.mu_all;

mu_all = load(strcat(theout,thedoc,"VI","_mu_all.mat"));
mu_all_VI = mu_all.mu_all;


figure(7),clf;

box on;
plot(1:num_epoch, mu_true*ones([num_epoch,1]),'-', 'linewidth', 1.5, 'Color', [0 0 0 0.3]);
grid on;
hold on;
plot(1:num_epoch, mu_all_GD,'-.', 'linewidth', 1, 'Color', [1,0,0,0.5]);
plot(1:num_epoch, mu_all_VI,'--', 'linewidth', 1, 'Color', [0,0,1,0.5]);

xlabel('{Epoch}','interpreter','latex');
ylabel('{$\mu$}','interpreter','latex');
ylim([0.19,0.28])
lgd = legend('true','TULIK-GD','TULIK-VI','Location','best');
lgd.NumColumns = 1;
set(gca,'FontSize',8);
ax = gca;
figH = gcf;
set(figH, 'Units', 'points','OuterPosition', [0 0 235 235]) % standard size: 19.7 17.5
exportgraphics(ax,strcat(theplot,thedoc,"_mu_all",".pdf")); 