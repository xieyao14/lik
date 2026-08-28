# Hawkes Processes with Event Time Uncertainty

This repository contains MATLAB code for the synthetic experiments in the paper *Hawkes Processes with Event Time Uncertainty*. The proposed model is referred to as the **Time Uncertain Latent Influence Kernel (TULIK)** point process. Two estimation procedures are implemented:

- **TULIK-VI**: variational-inequality update;
- **TULIK-GD**: gradient-descent update.

The code currently covers the synthetic time-only and graph experiments. The real-data experiments on sepsis and Atlanta burglary data are not included in the present repository snapshot.

## Experiment organization

The experiment prefixes have the following meanings:

| Prefix | Experiment | Discretization |
| --- | --- | --- |
| `f0` | Stationary time-only process | $N=32$, $N'=16$ |
| `f1` | High-dimensional nonstationary time-only process | $N=320$, $N'=80$ |
| `f2` | Low-dimensional nonstationary time-only process | $N=32$, $N'=8$ |
| `f3` | Nonstationary point process on a five-node graph | $V=5$, $N=32$, $N'=8$ |

## Paper-to-code map

| Paper result | Experiment | TULIK driver | Baseline drivers | Aggregation script |
| --- | --- | --- | --- | --- |
| Figures 4 and 5, Table 1, Figure A.1 | Nonstationary time-only, $N=32$, $N'=8$ | `test_f2_timeonly_lowdim.m` | `test_f2_timeonly_lowdim_GLM.m`, `test_f2_timeonly_lowdim_HPE.m` | `test_f2_table.m` |
| Figures 6 and 7, Table 2 | Nonstationary time-only, $N=320$, $N'=80$ | `test_f1_timeonly_highdim.m` | `test_f1_timeonly_highdim_GLM.m`, `test_f1_timeonly_highdim_HPE.m` | `test_f1_table.m` |
| Figures A.2 and A.3, Table A.1 | Stationary time-only, $N=32$, $N'=16$ | `test_f0_timeonly_stationary.m` | `test_f0_timeonly_stationary_GLM.m`, `test_f0_timeonly_stationary_HPE.m` | `test_f0_table.m` |
| Figures 8, 9, and 10, Table 3, Figure A.4 | Five-node graph experiment | `test_f3_graph.m` | `test_f3_graph_GLM.m`, `test_f3_graph_HPE.m` | `test_f3_table.m` |
| Figure A.5 | GD instability on graph data | `test_f3_graph_GD_instable.m` | Uses stable VI and GD outputs from `test_f3_graph.m` | Not applicable |

Figures 1 through 3 and Tables 4 and A.2 are conceptual illustrations or data-definition tables and do not correspond to synthetic experiment drivers.

The following real-data results do not have corresponding code in the current snapshot:

| Paper result | Dataset |
| --- | --- |
| Figures 11 and 12, Table 5, Figure A.6 | Sepsis-associated derangements data |
| Figure 13, Table 6, Figure A.7 | Atlanta burglary data |

## Main experiment drivers

### Stationary time-only experiment (`f0`)

`test_f0_timeonly_stationary.m` imposes a time-invariant kernel by estimating a single kernel row and replicating it across time. It generates:

- the stationary kernel comparison in Figure A.2;
- probability predictions in Figure A.3;
- kernel, baseline-intensity, and prediction-error files used in Table A.1;
- training likelihood and baseline-intensity trajectories.

The corresponding baseline scripts are:

- `test_f0_timeonly_stationary_GLM.m` for GLM-L and GLM-S;
- `test_f0_timeonly_stationary_HPE.m` for the exponential Hawkes-process baseline.

`test_f0_table.m` aggregates ten independent repetitions and prints the rows of Table A.1.

### High-dimensional time-only experiment (`f1`)

`test_f1_timeonly_highdim.m` implements the nonstationary experiment with $N=320$ and $N'=80$. It generates the kernel estimates in Figure 6, the probability predictions in Figure 7, and the TULIK entries in Table 2.

After training, the script applies truncated SVD to the estimated kernel. The manuscript uses thresholds $0.6$ for TULIK-VI and $0.8$ for TULIK-GD, producing rank-three estimates.

The baseline and table scripts are:

- `test_f1_timeonly_highdim_GLM.m`;
- `test_f1_timeonly_highdim_HPE.m`;
- `test_f1_table.m`.

### Low-dimensional time-only experiment (`f2`)

`test_f2_timeonly_lowdim.m` implements the nonstationary experiment with $N=32$ and $N'=8$. It generates:

- the true and estimated kernels in Figure 4;
- the prediction curves in Figure 5;
- the error quantities in Table 1;
- the training dynamics in Figure A.1.

The baseline and table scripts are:

- `test_f2_timeonly_lowdim_GLM.m`;
- `test_f2_timeonly_lowdim_HPE.m`;
- `test_f2_table.m`.

Additional `f2` scripts are development or sensitivity experiments and are not directly reported in the current paper:

- `test_f2_timeonly_lowdim_stationary.m`: fits a stationary kernel restriction to the nonstationary simulated data;
- `test_f2_timeonly_lowdim_largerntr.m`: uses a larger training set;
- `test_f2_timeonly_lowdim_ablation.m`: varies barrier and smoothness parameters.

### Graph experiment (`f3`)

`test_f3_graph.m` implements the five-node graph experiment. It generates:

- the directed graph in Figure 8;
- the true and recovered kernel tensors in Figure 9;
- the probability curves in Figure 10;
- the TULIK entries in Table 3;
- the likelihood and baseline-intensity convergence plots in Figure A.4.

The post-training SVD thresholds are $0.8$ for TULIK-VI and $1.0$ for TULIK-GD, producing rank-two estimates.

The baseline and table scripts are:

- `test_f3_graph_GLM.m`;
- `test_f3_graph_HPE.m`;
- `test_f3_table.m`.

`test_f3_graph_GD_instable.m` reruns graph GD with the larger VI learning rate and produces Figure A.5. The plotting section also requires the stable VI and GD outputs from `test_f3_graph.m`.

## Shared functions

| Function | Purpose |
| --- | --- |
| `kernel_timeonly_highrank.m` | Constructs the nonstationary time-only influence kernel associated with equation (80). |
| `truf_test141.m` | Constructs the graph influence kernel associated with equation (81). |
| `matrix_downsize.m` | Averages a fine kernel grid onto the discrete experimental grid. |
| `make_1d_laplacian.m` | Constructs the one-dimensional graph Laplacian used by the temporal smoothness penalty. |
| `network_proj.m` | Constructs event-history projection matrices on the full time interval. |
| `network_proj_mod.m` | Constructs projection matrices restricted to the observed prediction interval. |
| `network_proj_mod_speed_N_graph2.m` | Vectorized projection routine used by the main time-only and graph training loops. |
| `search_mu_timeonly.m` | Computes the time-only baseline-intensity update by bisection. |
| `search_mu_graph_linear.m` | Computes nodewise baseline-intensity updates for graph data. |

## Running the synthetic experiments

### 1. Prepare the directory structure

The scripts write data, numerical results, and plots to separate directories. A recommended repository structure is:

```text
.
|-- Input/
|-- Output/
|-- Plots/
|-- README.md
|-- test_f0_*.m
|-- test_f1_*.m
|-- test_f2_*.m
|-- test_f3_*.m
`-- supporting MATLAB functions
```

Replace the user-specific absolute paths near the beginning of each legacy script with paths defined relative to the script location. For example:

```matlab
scriptDir = string(fileparts(mfilename("fullpath")));
thein = fullfile(scriptDir, "Input") + string(filesep);
theout = fullfile(scriptDir, "Output") + string(filesep);
theplot = fullfile(scriptDir, "Plots") + string(filesep);
```

### 2. Generate TULIK-VI and TULIK-GD results

Each main TULIK script uses the flag `use_VI` to select the estimator. Run the script twice:

```matlab
use_VI = 1;  % TULIK-VI
```

and

```matlab
use_VI = 0;  % TULIK-GD
```

The scripts use `label = "VI"` and `label = "GD"` when naming their outputs.

### 3. Generate GLM-L and GLM-S results

Each GLM script uses `use_GLMI` to select the link function. In the current code naming convention:

```matlab
use_GLMI = 1;  % GLM-L, stored with label "GLMI"
use_GLMI = 0;  % GLM-S, stored with label "GLMS"
```

Run each GLM driver once for each setting.

### 4. Generate HP-E results

Run the corresponding `_HPE.m` driver once for each experiment family.

### 5. Aggregate repeated runs

Tables 1, 2, 3, and A.1 report means and standard deviations over ten repetitions. The table scripts expect output files indexed by process identifiers `0:9`. Each repeated experiment should use a distinct random seed and include its process identifier in every output filename.

After all repetitions and methods have completed, run the corresponding table script. It reads the saved relative errors, computes means and standard deviations, multiplies them by 100, and prints LaTeX table rows.

## Output files

The proposed-method drivers save the following main objects:

| Filename suffix | Contents |
| --- | --- |
| `VI.mat` or `GD.mat` | Estimated kernel parameter `X` |
| `_mu.mat` | Estimated baseline intensity `mu` |
| `_EstKerRelErr.mat` | Relative kernel errors under $\ell_1$, $\ell_2$, and $\ell_\infty$ norms |
| `_EstMuRelErr.mat` | Relative baseline-intensity errors |
| `_ProbPredErr.mat` | Absolute and relative probability-prediction errors |
| `_TrainLogLike.mat` | Negative log-likelihood over epochs |
| `_mu_all.mat` | Baseline-intensity trajectory over epochs |

The table scripts use the last three relative-error components of `proberror`, corresponding to relative $\ell_1$, $\ell_2$, and $\ell_\infty$ prediction errors.

## Current reproduction notes

The code-to-experiment correspondence is clear, but the current repository snapshot should be reconciled with the manuscript before claiming exact reproduction.

1. **Time-only kernel definition.** Equation (80) in the manuscript contains 13 summands and no separate base-kernel term. The current `kernel_timeonly_highrank.m` contains 15 summands and adds a `0.02` exponential term.

2. **Time discretization.** The manuscript's time-only horizon has total length 20, suggesting $h=20/(N+N')$. The current `f1` and `f2` drivers use `h=20/N` while constructing the kernel on an $(N+N')$-interval grid.

3. **Low-dimensional hyperparameters.** The manuscript reports smoothness weight $\delta_s=0.08$ and learning rates $0.4/0.2$ for VI and $0.2/0.1$ for GD. The current `f2` driver uses `smooth_weight = 0.2` and divides the learning-rate schedules by 10.

4. **Graph time scale and smoothness.** The manuscript uses the horizon $[-0.8,3.2]$, with total length 4 and $h=0.1$. The current graph driver sets `ut=8`, giving `h=0.2`, while its kernel grid spans length 4. The manuscript reports $\delta_s=0.004$, whereas the current graph driver uses `smooth_weight = 0.1`.

5. **Repeated-run naming.** The table scripts expect ten process-indexed output sets, but the local experiment drivers currently fix `rng(2024)` and do not add a process identifier to output filenames. Cluster-ready drivers or a common process-ID wrapper are needed.

6. **Graph parameters.** `truf_test141.m` loads `peak.mat` and `freq.mat`. These files must be included to reproduce the graph kernel.

7. **MATLAB function filenames.** MATLAB requires a function file to have the same name as its declared function. Remove download suffixes such as `(1)`, `(2)`, or `(4)` from helper filenames before running the code.

8. **Figure A.5 plotting.** The plotting section of `test_f3_graph_GD_instable.m` appears after a `return` statement. Execute that section separately or remove the early return after the required stable and unstable outputs have been generated.

## MATLAB requirements

The code requires MATLAB with functionality for sparse matrices, graph objects, SVD, and figure export. The Statistics and Machine Learning Toolbox may be required by some sampling or evaluation operations, depending on the MATLAB release.

The graph experiment also requires:

```text
peak.mat
freq.mat
```

## Reproduction status

| Component | Status |
| --- | --- |
| Synthetic time-only experiment drivers | Included |
| Synthetic graph experiment drivers | Included, but `peak.mat` and `freq.mat` are required |
| GLM-L, GLM-S, and HP-E baselines | Included |
| Ten-repetition table aggregation | Included, but matching process-indexed drivers are needed |
| Sepsis experiment | Not included |
| Atlanta burglary experiment | Not included |
| Exact agreement between all manuscript and code settings | Requires resolution of the items above |

## Citation

If you use this code, please cite the accompanying paper. Complete citation information will be added after publication.

