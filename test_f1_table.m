%%

theout = "/Users/leviathaniety/Dropbox (GaTech)/PROJ-PointProcessWithUncertainty/Codes/CodePackage/cluster_output/";

table_show = zeros([6,9]);
table_full = zeros([3,9,10]);

PID_list = 0:9;

for process_id = PID_list

    thedoc = strcat("test_f1_timeonly_highdim_PID",num2str(process_id));
    
    label = "VI"

    mu_rela_err = load(strcat(theout,thedoc,label,"_EstMuRelErr.mat"));
    mu_rela_err = mu_rela_err.mu_rela_err;
    table_full(:,1,process_id+1) = mu_rela_err;
    
    kernel_rela_err = load(strcat(theout,thedoc,label,"_EstKerRelErr.mat"));
    kernel_rela_err = kernel_rela_err.kernel_rela_err;
    table_full(:,3,process_id+1) = kernel_rela_err;
    
    proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
    proberror = proberror.proberror;
    mean_proberr = mean(proberror,1);
    table_full(:,5,process_id+1) = mean_proberr(4:6);

    label = "GD"
    

    mu_rela_err = load(strcat(theout,thedoc,label,"_EstMuRelErr.mat"));
    mu_rela_err = mu_rela_err.mu_rela_err;
    table_full(:,2,process_id+1) = mu_rela_err;
    
    kernel_rela_err = load(strcat(theout,thedoc,label,"_EstKerRelErr.mat"));
    kernel_rela_err = kernel_rela_err.kernel_rela_err;
    table_full(:,4,process_id+1) = kernel_rela_err;
    
    
    proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
    proberror = proberror.proberror;
    mean_proberr = mean(proberror,1);
    table_full(:,6,process_id+1) = mean_proberr(4:6);


    thedoc = strcat("test_f1_timeonly_highdim_GLM_PID",num2str(process_id));

    label = "GLMI"


    proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
    proberror = proberror.proberror;
    mean_proberr = mean(proberror,1);
    table_full(:,7,process_id+1) = mean_proberr(4:6);

    
    label = "GLMS"
    
    
    proberror = load(strcat(theout,thedoc,label,"_ProbPredErr",".mat"));
    proberror = proberror.proberror;
    mean_proberr = mean(proberror,1);
    table_full(:,8,process_id+1) = mean_proberr(4:6);


    thedoc = strcat("test_f1_timeonly_highdim_HPE_PID",num2str(process_id));

    proberror = load(strcat(theout,thedoc,"_ProbPredErr",".mat"));
    proberror = proberror.proberror;
    mean_proberr = mean(proberror,1);
    table_full(:,9,process_id+1) = mean_proberr(4:6);
    
end


%%
table_show([1,3,5],:) = mean(table_full,3);
table_show([2,4,6],:) = std(table_full,0,3);
table_show = table_show*100;
table_show = round(table_show,2);

norm_list = ["\\ell_1", "\\ell_1", "\\ell_2", "\\ell_2", "\\ell_\\infty", "\\ell_\\infty"];

for row = 1:6
    table_row = '';
    if mod(row,2)==1
        table_row = strcat('\\multirow{2}{*} {$',norm_list(row),'$}');
        for col = 1:9
            table_row = strcat(table_row,' & ',num2str(table_show(row,col)));
        end
        table_row = strcat(table_row,'\\\\', '\n');
    else
        for col = 1:9
            table_row = strcat(table_row,' & ',' (', num2str(table_show(row,col)),')');
        end
        table_row = strcat(table_row,'\\\\','[3pt]', '\n');
    end
    table_row = strcat(table_row, '\n');
    fprintf(table_row);
end





