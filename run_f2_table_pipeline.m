function pipeline_report = run_f2_table_pipeline(options)
%RUN_F2_TABLE_PIPELINE Run and aggregate the current f2 method scripts.
%
% The pipeline owns replica identifiers, seeds, output names, resumability,
% and aggregation. The existing method files are read and executed without
% being edited. Legacy-only adjustments are limited to run configuration,
% removal of an unused allocation, and continuation into existing metric
% sections located after an early return.

if nargin < 1
    options = struct();
end
validateattributes(options, {'struct'}, {'scalar'});

script_dir = string(fileparts(mfilename("fullpath")));
config = parse_pipeline_options(options, script_dir);

if ~isfolder(config.output_dir)
    mkdir(config.output_dir);
end
if ~isfolder(config.plot_dir)
    mkdir(config.plot_dir);
end

previous_figure_visibility = get(groot, "defaultFigureVisible");
figure_cleanup = onCleanup(@() set(groot, ...
    "defaultFigureVisible", previous_figure_visibility));
set(groot, "defaultFigureVisible", "off");

num_runs = numel(config.replica_ids)*numel(config.methods);
run_records = repmat(struct( ...
    'replica_id', nan, 'seed', nan, 'method', "", 'status', "", ...
    'elapsed_seconds', nan, 'message', "", 'required_files', strings(0, 1)), ...
    num_runs, 1);
record_index = 0;

fprintf("f2 pipeline: %d replicas, %d methods, %d epochs per method.\n", ...
    numel(config.replica_ids), numel(config.methods), config.num_epochs);

for replica_id = config.replica_ids
    seed = config.base_seed + replica_id - 1;

    for method = config.methods
        record_index = record_index + 1;
        required_files = required_method_files( ...
            config.output_dir, replica_id, method);
        run_records(record_index).replica_id = replica_id;
        run_records(record_index).seed = seed;
        run_records(record_index).method = method;
        run_records(record_index).required_files = required_files;

        if ~config.overwrite && all(isfile(required_files))
            run_records(record_index).status = "skipped";
            run_records(record_index).elapsed_seconds = 0;
            run_records(record_index).message = "all required outputs already exist";
            fprintf("Replica %02d | %-8s | skipped (outputs exist)\n", ...
                replica_id, method);
            continue;
        end

        fprintf("Replica %02d | %-8s | running with seed %d...\n", ...
            replica_id, method, seed);
        timer = tic;
        try
            captured_output = run_current_method( ...
                method, replica_id, seed, config, script_dir);
            elapsed_seconds = toc(timer);
            missing_files = required_files(~isfile(required_files));

            if isempty(missing_files)
                status = "completed";
                message = "";
            else
                status = "failed";
                message = "method finished without creating: " + ...
                    strjoin(missing_files, ", ");
            end

            run_records(record_index).status = status;
            run_records(record_index).elapsed_seconds = elapsed_seconds;
            run_records(record_index).message = message;
            fprintf("Replica %02d | %-8s | %s in %.1f s\n", ...
                replica_id, method, status, elapsed_seconds);

            if config.verbose && strlength(captured_output) > 0
                fprintf("%s\n", captured_output);
            end
            if status == "failed" && ~config.continue_on_error
                error("TULIK:MissingMethodOutputs", "%s", message);
            end
        catch exception
            elapsed_seconds = toc(timer);
            run_records(record_index).status = "failed";
            run_records(record_index).elapsed_seconds = elapsed_seconds;
            run_records(record_index).message = string(exception.message);
            fprintf(2, "Replica %02d | %-8s | failed in %.1f s: %s\n", ...
                replica_id, method, elapsed_seconds, exception.message);
            if ~config.continue_on_error
                rethrow(exception);
            end
        end

        close all force;
    end
end

aggregation_options = struct( ...
    'output_dir', config.output_dir, ...
    'replica_ids', config.replica_ids, ...
    'print_latex', config.print_latex, ...
    'save_summary', config.save_summary, ...
    'require_complete', false);
[table_show, table_full, aggregation_report] = ...
    test_f2_table(aggregation_options);

pipeline_report = struct();
pipeline_report.config = config;
pipeline_report.runs = run_records;
pipeline_report.table_show = table_show;
pipeline_report.table_full = table_full;
pipeline_report.aggregation = aggregation_report;
pipeline_report.num_completed = sum([run_records.status] == "completed");
pipeline_report.num_skipped = sum([run_records.status] == "skipped");
pipeline_report.num_failed = sum([run_records.status] == "failed");

fprintf("f2 pipeline finished: %d completed, %d skipped, %d failed.\n", ...
    pipeline_report.num_completed, pipeline_report.num_skipped, ...
    pipeline_report.num_failed);
clear figure_cleanup;
end

function config = parse_pipeline_options(options, script_dir)
allowed_options = ["replica_ids", "methods", "base_seed", ...
    "num_trajectories", "num_train", "num_test", "num_epochs", ...
    "batch_size", "output_dir", "plot_dir", "overwrite", ...
    "continue_on_error", "print_latex", "save_summary", "verbose"];
unknown_options = setdiff(string(fieldnames(options)), allowed_options);
assert(isempty(unknown_options), "Unknown pipeline option(s): %s", ...
    strjoin(unknown_options, ", "));

config = struct();
config.replica_ids = pipeline_option(options, "replica_ids", 1:10);
config.methods = replace(upper(string(pipeline_option(options, "methods", ...
    ["TULIK_VI", "TULIK_GD", "GLM_L", "GLM_S", "HPE"]))), "-", "_");
config.base_seed = pipeline_option(options, "base_seed", 2024);
config.num_trajectories = pipeline_option(options, "num_trajectories", 40000);
config.num_train = pipeline_option(options, "num_train", 16000);
config.num_test = pipeline_option(options, "num_test", 500);
config.num_epochs = pipeline_option(options, "num_epochs", 300);
config.batch_size = pipeline_option(options, "batch_size", 400);
config.output_dir = string(pipeline_option(options, "output_dir", ...
    fullfile(script_dir, "Output")));
config.plot_dir = string(pipeline_option(options, "plot_dir", ...
    fullfile(script_dir, "Plots")));
config.overwrite = logical(pipeline_option(options, "overwrite", false));
config.continue_on_error = logical(pipeline_option( ...
    options, "continue_on_error", true));
config.print_latex = logical(pipeline_option(options, "print_latex", true));
config.save_summary = logical(pipeline_option(options, "save_summary", false));
config.verbose = logical(pipeline_option(options, "verbose", false));

validateattributes(config.replica_ids, {'numeric'}, ...
    {'vector', 'integer', 'positive', 'finite', 'nonempty'});
config.replica_ids = reshape(config.replica_ids, 1, []);
assert(numel(unique(config.replica_ids)) == numel(config.replica_ids), ...
    "Replica IDs must be unique.");

valid_methods = ["TULIK_VI", "TULIK_GD", "GLM_L", "GLM_S", "HPE"];
config.methods = reshape(config.methods, 1, []);
assert(~isempty(config.methods) && all(ismember(config.methods, valid_methods)), ...
    "Methods must be selected from: %s", strjoin(valid_methods, ", "));
assert(numel(unique(config.methods)) == numel(config.methods), ...
    "Methods must be unique.");

validateattributes(config.base_seed, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'});
validateattributes(config.num_trajectories, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
assert(config.num_trajectories >= 18, ...
    "num_trajectories must be at least 18 for the legacy diagnostic plots.");
validateattributes(config.num_train, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(config.num_test, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(config.num_epochs, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<=', 300});
validateattributes(config.batch_size, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
assert(config.num_trajectories >= config.num_train + config.num_test, ...
    "num_trajectories must cover the requested training and test sets.");
assert(config.batch_size <= config.num_train && ...
    mod(config.num_train, config.batch_size) == 0, ...
    "batch_size must divide num_train exactly.");
assert(config.base_seed + max(config.replica_ids) - 1 <= 2^32 - 1, ...
    "Replica seed exceeds MATLAB's supported range.");

logical_fields = ["overwrite", "continue_on_error", "print_latex", ...
    "save_summary", "verbose"];
for field = logical_fields
    validateattributes(config.(field), {'logical'}, {'scalar'});
end
assert(isscalar(config.output_dir) && strlength(config.output_dir) > 0, ...
    "output_dir must be a nonempty string scalar.");
assert(isscalar(config.plot_dir) && strlength(config.plot_dir) > 0, ...
    "plot_dir must be a nonempty string scalar.");
end

function value = pipeline_option(options, name, default_value)
if isfield(options, name)
    value = options.(name);
else
    value = default_value;
end
end

function captured_output = run_current_method( ...
    method, replica_id, seed, config, script_dir)
switch method
    case {"TULIK_VI", "TULIK_GD"}
        captured_output = run_current_tulik( ...
            method, replica_id, seed, config, script_dir);
    case {"GLM_L", "GLM_S"}
        captured_output = run_current_glm( ...
            method, replica_id, seed, config, script_dir);
    case "HPE"
        captured_output = run_current_hpe( ...
            replica_id, seed, config, script_dir);
    otherwise
        error("TULIK:UnknownF2Method", "Unknown f2 method: %s", method);
end
end

function captured_output = run_current_tulik( ...
    method, replica_id, seed, config, script_dir)
source = fileread(fullfile(script_dir, "test_f2_timeonly_lowdim.m"));
source = configure_output_directories(source, config, "TULIK");
source = replace_once(source, "clearvars;", ...
    "% Workspace clearing disabled by run_f2_table_pipeline.", ...
    "TULIK workspace setup");
source = replace_once(source, 'rng(2024, "twister");', ...
    sprintf('rng(%d, "twister");', seed), "TULIK random seed");
source = configure_legacy_common( ...
    source, replica_id, config, "TULIK");

use_vi = method == "TULIK_VI";
source = replace_once(source, "use_VI = 0;", ...
    sprintf('use_VI = %d;', use_vi), "TULIK method selector");
source = configure_epoch_schedule(source, config, "TULIK");
source = truncate_before_return(source, 1, "TULIK"); %#ok<NASGU>

captured_output = evalc('eval(source);');
end

function captured_output = run_current_glm( ...
    method, replica_id, seed, config, script_dir)
source = fileread(fullfile(script_dir, "test_f2_timeonly_lowdim_GLM.m"));
source = configure_output_directories(source, config, "GLM");
source = replace_once(source, "clear all; rng(2024);", ...
    sprintf('rng(%d);', seed), "GLM random seed");
source = configure_legacy_common(source, replica_id, config, "GLM");

use_glm_l = method == "GLM_L";
source = replace_once(source, "use_GLMI = 1;", ...
    sprintf('use_GLMI = %d;', use_glm_l), "GLM method selector");
source = replace_once(source, 'label = "GLMI";', 'label = "GLM_L";', ...
    "GLM-L output label");
source = replace_once(source, 'label = "GLMS";', 'label = "GLM_S";', ...
    "GLM-S output label");
source = replace_expected(source, 'strcmp(label, "GLMI")', ...
    'use_GLMI', 3, "GLM internal method checks");
source = configure_epoch_schedule(source, config, "GLM");
source = continue_after_first_return(source, "GLM"); %#ok<NASGU>

captured_output = evalc('eval(source);');
end

function captured_output = run_current_hpe( ...
    replica_id, seed, config, script_dir)
source = fileread(fullfile(script_dir, "test_f2_timeonly_lowdim_HPE.m"));
source = configure_output_directories(source, config, "HPE");
source = replace_once(source, "clear all; rng(2024);", ...
    sprintf('rng(%d);', seed), "HP-E random seed");
source = configure_legacy_common(source, replica_id, config, "HPE");
source = replace_once(source, "num_epoch = 300;", ...
    sprintf('num_epoch = %d;', config.num_epochs), "HP-E epoch count");
source = replace_once(source, "batch_size = 400;", ...
    sprintf('batch_size = %d;', config.batch_size), ...
    "HP-E batch size"); %#ok<NASGU>

captured_output = evalc('eval(source);');
end

function source = configure_output_directories(source, config, family)
if family == "TULIK"
    output_marker = ...
        'theout  = fullfile(scriptDir, "Output") + string(filesep);';
    plot_marker = ...
        'theplot = fullfile(scriptDir, "Plots")  + string(filesep);';
else
    output_marker = ...
        'theout = fullfile(scriptDir, "Output") + string(filesep);';
    plot_marker = ...
        'theplot = fullfile(scriptDir, "Plots") + string(filesep);';
end

source = replace_once(source, output_marker, ...
    sprintf('theout = %s + string(filesep);', ...
    matlab_string_literal(config.output_dir)), family + " output directory");
source = replace_once(source, plot_marker, ...
    sprintf('theplot = %s + string(filesep);', ...
    matlab_string_literal(config.plot_dir)), family + " plot directory");
end

function source = configure_legacy_common( ...
    source, replica_id, config, family)
source = replace_once(source, "M = 40000; %40000;", ...
    sprintf('M = %d;', config.num_trajectories), family + " trajectory count");
source = replace_once(source, "ntr = 16000; %32000; %16000;", ...
    sprintf('ntr = %d;', config.num_train), family + " training count");
source = replace_once(source, "nte = min(500,M -ntr);", ...
    sprintf('nte = min(%d, M-ntr);', config.num_test), family + " test count");
source = replace_once(source, ...
    "eta_ob = false(M, Nprime+N, N);", ...
    "% Unused eta_ob allocation omitted by run_f2_table_pipeline.", ...
    family + " unused eta allocation");
source = replace_once(source, ...
    "eta_ob(:, t:t+Nprime-1, t)= ypre;", ...
    "% Unused eta_ob assignment omitted by run_f2_table_pipeline.", ...
    family + " unused eta assignment");

switch family
    case "TULIK"
        old_document = 'thedoc = "test_f2_timeonly_lowhdim";';
        new_document = sprintf( ...
            'thedoc = "test_f2_timeonly_lowdim_replica_%02d_TULIK_";', ...
            replica_id);
    case "GLM"
        old_document = 'thedoc = "test_f2_timeonly_lowhdim_GLM";';
        new_document = sprintf( ...
            'thedoc = "test_f2_timeonly_lowdim_replica_%02d_";', ...
            replica_id);
    case "HPE"
        old_document = 'thedoc = "test_f2_timeonly_lowdim_HPE";';
        new_document = sprintf( ...
            'thedoc = "test_f2_timeonly_lowdim_replica_%02d_HPE";', ...
            replica_id);
    otherwise
        error("TULIK:UnknownF2Family", ...
            "Unknown f2 method family: %s", family);
end
source = replace_once(source, old_document, new_document, ...
    family + " output document name");
end

function source = configure_epoch_schedule(source, config, family)
source = replace_once(source, ...
    'bs_schedule = [400*ones(100,1), 400*ones(100,1), 400*ones(100,1) ];', ...
    sprintf('bs_schedule = %d*ones(size(lr_schedule));', config.batch_size), ...
    family + " batch schedule");
source = replace_once(source, 'num_epoch = numel(lr_schedule );', ...
    sprintf(['num_epoch = min(%d, numel(lr_schedule));\n' ...
    'lr_schedule = lr_schedule(1:num_epoch);\n' ...
    'bs_schedule = bs_schedule(1:num_epoch);'], config.num_epochs), ...
    family + " epoch limit");
end

function source = replace_once(source, old_text, new_text, description)
source_count = count(string(source), string(old_text));
if source_count ~= 1
    error("TULIK:LegacyAdapterMismatch", ...
        "Expected one %s marker but found %d.", description, source_count);
end
source = strrep(source, char(old_text), char(new_text));
end

function source = replace_expected( ...
    source, old_text, new_text, expected_count, description)
source_count = count(string(source), string(old_text));
if source_count ~= expected_count
    error("TULIK:LegacyAdapterMismatch", ...
        "Expected %d %s markers but found %d.", ...
        expected_count, description, source_count);
end
source = strrep(source, char(old_text), char(new_text));
end

function source = truncate_before_return(source, return_number, description)
[starts, ~] = regexp(source, '(?m)^[ \t]*return;[ \t]*\r?$', 'start', 'end');
if numel(starts) < return_number
    error("TULIK:LegacyAdapterMismatch", ...
        "%s script has fewer than %d return statements.", ...
        description, return_number);
end
source = source(1:starts(return_number)-1);
end

function source = continue_after_first_return(source, description)
[starts, ends] = regexp(source, ...
    '(?m)^[ \t]*return;[ \t]*\r?$', 'start', 'end');
if numel(starts) < 2
    error("TULIK:LegacyAdapterMismatch", ...
        "%s script must contain at least two return statements.", description);
end
source = [source(1:starts(1)-1), ...
    sprintf('%% Continued by run_f2_table_pipeline.\n'), ...
    source(ends(1)+1:starts(2)-1)];
end

function literal = matlab_string_literal(value)
value = string(value);
assert(isscalar(value) && ~contains(value, '"'), ...
    "Paths containing double quotes are not supported.");
literal = ['"', char(value), '"'];
end

function files = required_method_files(output_dir, replica_id, method)
prefix = string(fullfile(output_dir, sprintf( ...
    'test_f2_timeonly_lowdim_replica_%02d_', replica_id)));

switch method
    case "TULIK_VI"
        stem = prefix + "TULIK_VI";
        suffixes = ["_EstMuRelErr.mat", "_EstKerRelErr.mat", "_ProbPredErr.mat"];
    case "TULIK_GD"
        stem = prefix + "TULIK_GD";
        suffixes = ["_EstMuRelErr.mat", "_EstKerRelErr.mat", "_ProbPredErr.mat"];
    case "GLM_L"
        stem = prefix + "GLM_L";
        suffixes = "_ProbPredErr.mat";
    case "GLM_S"
        stem = prefix + "GLM_S";
        suffixes = "_ProbPredErr.mat";
    case "HPE"
        stem = prefix + "HPE";
        suffixes = "_ProbPredErr.mat";
    otherwise
        error("TULIK:UnknownF2Method", "Unknown f2 method: %s", method);
end
files = reshape(stem + suffixes, [], 1);
end
