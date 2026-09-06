function [table_show, table_full, report] = test_f1_table(options)
%TEST_F1_TABLE Aggregate paired replicas for the N=320, N'=80 experiment.
%
% Each replica is identified by a number (replica_01 through replica_10 by
% default). A replica contributes to the table only when all TULIK-VI,
% TULIK-GD, GLM-L, GLM-S, and HP-E metrics are present and valid. Missing or
% malformed artifacts are reported instead of stopping the aggregation.

if nargin < 1
    options = struct();
end
validateattributes(options, {'struct'}, {'scalar'});

script_dir = string(fileparts(mfilename("fullpath")));
output_dir = string(table_option(options, "output_dir", ...
    fullfile(script_dir, "Output")));
replica_ids = table_option(options, "replica_ids", 1:10);
print_latex = logical(table_option(options, "print_latex", true));
save_summary = logical(table_option(options, "save_summary", false));
require_complete = logical(table_option(options, "require_complete", false));
summary_path = string(table_option(options, "summary_path", ...
    fullfile(output_dir, "test_f1_table_summary.mat")));

validateattributes(replica_ids, {'numeric'}, ...
    {'vector', 'integer', 'nonnegative', 'finite', 'nonempty'});
validateattributes(print_latex, {'logical'}, {'scalar'});
validateattributes(save_summary, {'logical'}, {'scalar'});
validateattributes(require_complete, {'logical'}, {'scalar'});
replica_ids = reshape(replica_ids, 1, []);
assert(numel(unique(replica_ids)) == numel(replica_ids), ...
    "Replica IDs must be unique.");

column_names = ["TULIK-VI mu", "TULIK-GD mu", ...
    "TULIK-VI kernel", "TULIK-GD kernel", ...
    "TULIK-VI prediction", "TULIK-GD prediction", ...
    "GLM-L prediction", "GLM-S prediction", "HP-E prediction"];
norm_names = ["l1", "l2", "linf"];

% Columns: destination, display name, filename template, variable, type.
artifact_specs = {
    1, column_names(1), "test_f1_timeonly_highdim_replica_%02d_TULIK_VI_EstMuRelErr.mat",  "mu_rela_err",     "vector";
    2, column_names(2), "test_f1_timeonly_highdim_replica_%02d_TULIK_GD_EstMuRelErr.mat",  "mu_rela_err",     "vector";
    3, column_names(3), "test_f1_timeonly_highdim_replica_%02d_TULIK_VI_EstKerRelErr.mat", "kernel_rela_err", "vector";
    4, column_names(4), "test_f1_timeonly_highdim_replica_%02d_TULIK_GD_EstKerRelErr.mat", "kernel_rela_err", "vector";
    5, column_names(5), "test_f1_timeonly_highdim_replica_%02d_TULIK_VI_ProbPredErr.mat",  "proberror",       "prediction";
    6, column_names(6), "test_f1_timeonly_highdim_replica_%02d_TULIK_GD_ProbPredErr.mat",  "proberror",       "prediction";
    7, column_names(7), "test_f1_timeonly_highdim_replica_%02d_GLM_L_ProbPredErr.mat",     "proberror",       "prediction";
    8, column_names(8), "test_f1_timeonly_highdim_replica_%02d_GLM_S_ProbPredErr.mat",     "proberror",       "prediction";
    9, column_names(9), "test_f1_timeonly_highdim_replica_%02d_HPE_ProbPredErr.mat",       "proberror",       "prediction";
    };

num_replicas = numel(replica_ids);
table_full = nan(3, 9, num_replicas);
max_issue_count = num_replicas*size(artifact_specs, 1);
missing_artifacts = repmat(struct( ...
    'replica_id', nan, 'metric', "", 'path', "", 'reason', ""), ...
    max_issue_count, 1);
issue_count = 0;

for replica_index = 1:num_replicas
    replica_id = replica_ids(replica_index);
    for spec_index = 1:size(artifact_specs, 1)
        column = artifact_specs{spec_index, 1};
        metric = artifact_specs{spec_index, 2};
        filename = sprintf(artifact_specs{spec_index, 3}, replica_id);
        artifact_path = fullfile(output_dir, filename);
        variable_name = artifact_specs{spec_index, 4};
        metric_type = artifact_specs{spec_index, 5};

        [values, issue] = load_table_metric( ...
            artifact_path, variable_name, metric_type);
        table_full(:, column, replica_index) = values;

        if strlength(issue) > 0
            issue_count = issue_count + 1;
            missing_artifacts(issue_count) = struct( ...
                'replica_id', replica_id, ...
                'metric', metric, ...
                'path', string(artifact_path), ...
                'reason', issue);
        end
    end
end
missing_artifacts = missing_artifacts(1:issue_count);

complete_mask = reshape(all(isfinite(table_full), [1, 2]), 1, []);
complete_replica_ids = replica_ids(complete_mask);
incomplete_replica_ids = replica_ids(~complete_mask);

table_show = nan(6, 9);
if any(complete_mask)
    paired_results = table_full(:, :, complete_mask);
    table_show([1, 3, 5], :) = mean(paired_results, 3);
    table_show([2, 4, 6], :) = std(paired_results, 0, 3);
    table_show = round(100*table_show, 2);
end

if all(complete_mask)
    status = "complete";
else
    status = "incomplete";
end

report = struct();
report.status = status;
report.output_dir = output_dir;
report.requested_replica_ids = replica_ids;
report.complete_replica_ids = complete_replica_ids;
report.incomplete_replica_ids = incomplete_replica_ids;
report.num_requested = num_replicas;
report.num_complete = numel(complete_replica_ids);
report.column_names = column_names;
report.norm_names = norm_names;
report.filename_patterns = string(artifact_specs(:, 3)).';
report.missing_artifacts = missing_artifacts;
report.aggregation = "paired complete replicas only";
report.scale_factor = 100;
report.summary_path = summary_path;

print_aggregation_report(report);
if print_latex
    if isempty(complete_replica_ids)
        fprintf("LaTeX rows were not printed because no complete replicas are available.\n");
    else
        print_table_latex(table_show);
    end
end

if save_summary
    summary_dir = string(fileparts(summary_path));
    if strlength(summary_dir) > 0 && ~isfolder(summary_dir)
        mkdir(summary_dir);
    end
    save(summary_path, "table_show", "table_full", "report");
    fprintf("Saved aggregation summary: %s\n", summary_path);
end

if require_complete && status ~= "complete"
    error("TULIK:IncompleteTable2Inputs", ...
        "Only %d of %d requested replicas are complete.", ...
        report.num_complete, report.num_requested);
end
end

function value = table_option(options, name, default_value)
if isfield(options, name)
    value = options.(name);
else
    value = default_value;
end
end

function [values, issue] = load_table_metric(path, variable_name, metric_type)
values = nan(3, 1);
issue = "";

if ~isfile(path)
    issue = "missing file";
    return;
end

try
    loaded = load(path, variable_name);
catch exception
    issue = "load failed: " + string(exception.message);
    return;
end

if ~isfield(loaded, variable_name)
    issue = "missing variable '" + variable_name + "'";
    return;
end

raw = loaded.(variable_name);
if ~isnumeric(raw) || ~isreal(raw) || any(~isfinite(raw), "all")
    issue = "metric must be finite, real, and numeric";
    return;
end

switch metric_type
    case "vector"
        if numel(raw) ~= 3
            issue = "expected a three-element relative-error vector";
            return;
        end
        values = reshape(raw, 3, 1);
    case "prediction"
        if ~ismatrix(raw) || isempty(raw) || size(raw, 2) < 6
            issue = "expected a nonempty prediction-error matrix with at least six columns";
            return;
        end
        values = mean(raw(:, 4:6), 1).';
    otherwise
        error("TULIK:UnknownTableMetricType", ...
            "Unknown table metric type: %s", metric_type);
end
end

function print_aggregation_report(report)
fprintf("Table 2 aggregation: %d/%d complete paired replicas.\n", ...
    report.num_complete, report.num_requested);

for replica_id = report.incomplete_replica_ids
    issues = report.missing_artifacts( ...
        [report.missing_artifacts.replica_id] == replica_id);
    labels = string({issues.metric});
    fprintf("  Replica %02d incomplete: %s\n", replica_id, strjoin(labels, ", "));
end
end

function print_table_latex(table_show)
norm_labels = ["\ell_1", "\ell_2", "\ell_\infty"];
fprintf("LaTeX rows (means and sample standard deviations, scaled by 100):\n");

for norm_index = 1:3
    mean_row_index = 2*norm_index - 1;
    std_row_index = mean_row_index + 1;
    mean_row = "\multirow{2}{*}{$" + norm_labels(norm_index) + "$}";
    std_row = "";

    for column = 1:9
        if norm_index > 1 && column <= 2
            mean_value = "--";
            std_value = "";
        else
            mean_value = sprintf("%.2f", table_show(mean_row_index, column));
            std_value = "(" + sprintf("%.2f", table_show(std_row_index, column)) + ")";
        end
        mean_row = mean_row + " & " + mean_value;
        std_row = std_row + " & " + std_value;
    end

    mean_row = mean_row + " \\";
    std_row = std_row + " \\[3pt]";
    fprintf("%s\n%s\n", mean_row, std_row);
end
end
