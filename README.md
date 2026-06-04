Open-source companion code for the paper:

**Optimal Measurement of Traffic Hysteresis under Traffic Oscillations: A Binary Integer Programming Approach**

This repository implements methods for measuring traffic hysteresis loops from vehicle trajectory data. The proposed method represents trajectory points as a directed graph and selects wave propagation paths with a binary integer programming model. The selected paths define measurement regions for Edie's generalized traffic variables: flow, density, and speed.

## What This Repository Contains

- A binary integer programming implementation of the proposed measurement method.
- Baseline methods used in the paper: parallelogram measurement and traditional rectangular measurement.
- Analytical validation logic for homogeneous AV platoons under single-frequency oscillation.
- Curated NGSIM and OpenACC/OpenACC-derived input instances used in the numerical study.
- Plotting and summary utilities for flow-density diagrams, measurement-region plots, speed-variation summaries, loop magnitudes, center points, and solver logs.

## Method Map

| Paper concept | Code |
| --- | --- |
| Proposed minimum-cost network flow / binary integer programming method | `Solver.py` |
| Parallelogram benchmark | `Benchmark.py` |
| Traditional rectangular measurement baseline | `RectangleMethod.py` |
| Analytical single-frequency AV validation | `AnalyticalMethod.py` |
| Vehicle trajectory input readers | `ReadInput.py` |
| Shared experiment parameters and slot-window utilities | `CommonDefines.py` |
| Plotting and hysteresis metrics | `PlotHelper.py` |
| Command-line experiment driver | `main.py` |

## Dependencies

Python 3.10+ is recommended. The required Python packages are listed in `requirements.txt`:

- `numpy`
- `pandas`
- `matplotlib`
- `networkx`
- `plotly`
- `gurobipy`

The optimization model requires `gurobipy` and a valid Gurobi license. The paper experiments used **Gurobi 11.0.1** with a 300-second time limit. The model implementation surfaces missing solver configuration, infeasibility, and other solver failures directly.

## How To Run

The commands below are written for Windows PowerShell from the repository root. The `python main.py ...` commands are portable and can also be used on macOS or Linux after activating an equivalent virtual environment.

### 1. Get the code

```powershell
git clone <repository-url>
cd Hysteresis
```

If you already have the repository, open PowerShell in the existing checkout:

```powershell
cd D:\GitHub\Hysteresis
```

### 2. Create and activate a Python environment

Python 3.10 or newer is recommended.

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

If PowerShell blocks activation, allow local scripts for the current user and then run the activation command again:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
.\.venv\Scripts\Activate.ps1
```

### 3. Confirm Gurobi before optimization runs

Runs that include `--methods optimization` require `gurobipy` and a valid Gurobi license. Confirm the Python package can be imported before running those experiments:

```powershell
python -c "import gurobipy; print(gurobipy.gurobi.version())"
```

### 4. Inspect the command-line options

```powershell
python main.py --help
```

The main options are:

- `--scenario`: `analytical`, `ngsim`, or `openacc`.
- `--dataset`: input instance name for `ngsim` or `openacc`.
- `--methods`: comma-separated method list. Available methods are `analytical`, `benchmark`, `rectangle`, and `optimization`.
- `--output-dir`: folder for generated figures, CSV summaries, text summaries, and solver logs.
- `--show-plots`: display Matplotlib windows after saving figures.
- `--draw-wave-speed-heat-map`: add the analytical wave-speed heat map.

### 5. Run common experiments

Analytical smoke test without Gurobi:

```powershell
python main.py --scenario analytical --methods analytical --output-dir results/analytical_smoke
```

Analytical comparison with the optimization method:

```powershell
python main.py --scenario analytical --methods analytical,optimization --output-dir results/analytical_opt
```

NGSIM benchmark and optimization:

```powershell
python main.py --scenario ngsim --dataset Platoon_Position_3 --methods benchmark,optimization --output-dir results/ngsim/Platoon_Position_3
```

OpenACC/ZalaZONE benchmark and optimization:

```powershell
python main.py --scenario openacc --dataset dynamic_part17 --methods benchmark,optimization --output-dir results/openacc/dynamic_part17
```

Include the rectangular baseline by adding `rectangle` to `--methods`:

```powershell
python main.py --scenario analytical --methods analytical,rectangle,optimization --output-dir results/analytical_with_rectangle
```

Show figures interactively after they are saved:

```powershell
python main.py --scenario analytical --methods analytical --output-dir results/analytical_show --show-plots
```

### Adjust Analytical Scenario Parameters

The analytical single-frequency scenario is generated from the active parameter values in `CommonDefines.py`. The file `analytical/params.txt` records archived analytical settings, but `main.py` does not load it when running experiments.

To run a custom analytical case:

1. Edit the analytical parameters in `CommonDefines.py`.
2. Run the analytical scenario again:

```powershell
python main.py --scenario analytical --methods analytical --output-dir results/analytical_custom
```

3. Use a fresh `--output-dir` for each parameter set so generated figures and summaries are not mixed with earlier runs.

Key analytical parameters are:

| Parameter | Meaning |
| --- | --- |
| `a` | Leading-vehicle oscillation amplitude `A` in meters. |
| `omega` | Angular frequency in rad/s, for example `0.16 * np.pi`; `f = omega / (2 * np.pi)` is derived from it. |
| `phase_shift` | Oscillation phase shift `sigma` in radians. |
| `time_length`, `step_size` | Generated time horizon, equal to `time_length * step_size` seconds. |
| `platoon_size` | Number of vehicles in the homogeneous AV platoon. |
| `ks`, `kv` | Feedback gains in the transfer function. |
| `phi` | Actuation time lag for a vehicle to realize the desired acceleration. |
| `tau` | Minimal time gap in seconds. |
| `ve` | Equilibrium speed in m/s. |
| `standstill_s` | Standstill spacing in meters. |
| `xe` | Equilibrium spacing, normally kept as `tau * ve + standstill_s`. |
| `vehicle_point_sample_interval`, `shift_interval` | Measurement-region sampling controls, expressed as integer multiples of `step_size`. |

The paper's analytical tests vary the leading-vehicle oscillation by changing `a` and `omega`. For example, set a new amplitude and angular frequency in `CommonDefines.py`, then rerun the analytical command above to generate the corresponding trajectories, measurement regions, and flow-density outputs.

### 6. Available curated datasets

NGSIM instances listed in `Ngsim/Platoon_data/time_span.csv`:

OpenACC/ZalaZONE instances listed in `ACC/ZalaZONE/time_span.csv`:

### 7. Read the outputs

Each run writes outputs to the folder passed through `--output-dir`. If no output folder is supplied, the code writes under `results/` using scenario-specific defaults.

Typical generated files include:

- trajectory and measurement-region figures such as `Analytical.png`, `Benchmark.png`, `Rectangular.png`, or `Optimization.png`;
- `Flow-density diagram.png` and `speed standard error.png`;
- `qk_points_<Method>.csv` files;
- `HL_magnitudes.txt`, `HL_center_points.txt`, `avg_wave_speeds.txt`, and `compare_infos.txt`;
- optimization artifacts `model.lp`, `gurobi_log.txt`, and, for infeasible models, `model.ilp`.

## Data Layout

Curated input data and historical artifacts are included in the repository:

- `Ngsim/Platoon_data/` contains reconstructed NGSIM platoon trajectory CSV files and time-window metadata.
- `ACC/ZalaZONE/` contains curated OpenACC/ZalaZONE commercial ACC platoon instances and time-window metadata.
- `analytical/` contains analytical scenario parameters and archived analytical materials.
- `outputs/` and selected historical result folders contain retained artifacts from earlier analyses.

See [docs/data.md](docs/data.md) for more detail on input data, generated artifacts, and file naming.

## Output Files

The project can generate figures, CSV files, text summaries, and solver artifacts, including:

- `Analytical.png`, `Benchmark.png`, `Rectangular.png`, or `Optimization.png`: trajectory plots with measurement regions.
- `Flow-density diagram.png`: hysteresis loop comparison.
- `speed standard error.png`: within-region speed variation comparison.
- `qk_points_<Method>.csv`: density-flow points for each method.
- `HL_magnitudes.txt`: hysteresis loop magnitudes.
- `HL_center_points.txt`: loop center points.
- `avg_wave_speeds.txt`: average wave-speed summaries.
- `compare_infos.txt`: compact comparison table.
- `gurobi_log.txt` and `model.lp`: solver artifacts when optimization is selected.

New generated outputs should stay under `results/` or another explicit output folder. These generated files are intentionally ignored by Git.

## Citation

If you use this repository, please cite the paper. Replace the metadata below with the final publication details when available.

```bibtex
@article{pu2026optimal,
  title = {Optimal Measurement of Traffic Hysteresis under Traffic Oscillations: A Binary Integer Programming Approach},
  author = {Pu, Fan and Zhou, Yang and Ahn, Soyoung and Li, Sixu and Kontar, Wissam and Wang, Xiubin},
  journal = {Transportation Research Part B: Methodological},
  year = {2026}
}
```
