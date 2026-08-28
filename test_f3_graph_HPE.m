clear all; rng(2024);
addpath /Users/leviathaniety/Documents/MATLAB/mosek/9.3/toolbox/r2015a/;
addpath /Users/leviathaniety/Documents/MATLAB/Functions;
thein = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Input/";
theout = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Output/";
theplot = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/Plots/";
thedoc = "test_f3_graph_HPE";


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
mu = zeros(V,1); % initial base
for ik = 1:V
    mu(ik) = mean(mean(event_data_tr(:,(ik-1)*(N+Nprime)+1:ik*(N+Nprime)),2))/h;
end

event_data_te = event_data(idx_te,:);
            



lam_true_val = lambda_true(idx_te,:);
lambar_true_val = lambar_true(idx_te,:);
p_val_true = (1-exp(-h.*lambar_true_val))./lambar_true_val.*lam_true_val;



event_data_tensor = zeros(M,V,N);
for k1 = 1:V
    s = (k1-1)*N+1;
    e = k1*N;
    event_data_tensor(:,k1,:) = event_data(:,s:e);
end

%%
cont_event_data = cell(M,1);
tplt = h:h:ut;
for m = 1:M
    for k1 = 1:V
          idx = logical(event_data_tensor(m,k1,:));
          cont_event_data{m} = [cont_event_data{m},[tplt(idx);ones(size(tplt(idx)))*k1]];
    end
end

%% split data
event_data_tr = cont_event_data(idx_tr);
event_data_te = cont_event_data(idx_te);


%% training hyperparameters
num_epoch = 150;



batch_size = 800;
batch_size_schedule = batch_size*ones(num_epoch,1);
eta_schedule = 0.02*ones(num_epoch,1);

iepoch = 0;
likelihood_epoch = zeros(num_epoch,1);



AB = cell(3,1);
AB{1} = ones(V);
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
            lhlh2 = lklh2 + sum(mu,'all')*ut;

            ui_mask = ceil(Lambda{ib}(2,:));
            A_u_ui = A(:,ui_mask);
            tmp = kron((1-exp(-beta*(ut-edata{ib}(1,:)))),ones(V,1));
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

            dldmu = dldmu-ut;


            for je = 1:num_events
                k2 = ceil(lam_now(2,je)); 
                tj = edata{ib}(1,je);
                dldA(k2,:) = dldA(k2,:) - (1-exp(-beta*(ut-tj)));

                
            end
            ui_mask = ceil(Lambda{ib}(2,:));
            A_u_ui = A(:,ui_mask);
            tmp = kron((beta*exp(-beta*(ut-edata{ib}(1,:)))),ones(V,1));
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
nt = (N+Nprime)*5;
ntau = Nprime*5;
t_grid = ((0.5:1:nt)/nt)'*ut_kernel;
tau_grid = ((0.5:1:nt)/nt)'*ut_kernel;

psi0 = zeros(nt, ntau, V^2);
for j = 1:ntau
    for k1 = 1:V
        for k2 = 1:V
            a = A(k1,k2);
            psi0(:,j,(k2-1)*V+k1) = a*exp(-beta.*(t_grid - tau_grid(j)));
        end
    end
end

downfactor = 5;
psi_learn = zeros(nt/downfactor,ntau/downfactor,V^2);
for k = 1:V^2
    psi_learn(:,:,k) = matrix_downsize(psi0(:,:,k), downfactor);
end
X = psi_learn;

%%
event_data_te = event_data(idx_te,:);
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
std(proberror,1)
save(strcat(theout,thedoc,"_ProbPredErr",".mat"), "proberror");




