# Data and Artifact Layout

This repository contains source code and curated input data. Future generated files should be written to `results/` or another explicit output directory so the source tree stays readable.

## Curated Inputs

### NGSIM

`Ngsim/Platoon_data/` contains reconstructed platoon trajectory CSV files and `time_span.csv`, which records the start slot for each curated instance.

Curated NGSIM instances include:

- `Platoon_Position_3`
- `Platoon_Position_4`
- `Platoon_Position_7`
- `Platoon_Position_9`
- `Platoon_Position_12`
- `Platoon_Position_17`

### OpenACC / ZalaZONE

`ACC/ZalaZONE/` contains curated commercial ACC platoon instances and `time_span.csv`, which records the start slot and number of points for each instance.

Curated OpenACC/ZalaZONE instances include:

- `dynamic_part17`
- `dynamic_part1_ins1`
- `dynamic_part1_ins2`
- `dynamic_part1_ins3`
- `dynamic_part1_ins4`

### Analytical Scenario

The analytical scenario generates trajectories from parameters in `CommonDefines.py` and `AnalyticalMethod.py`. It does not require external trajectory data.

## Generated Artifacts

Generated artifacts include:

- comparison summaries: `compare_infos.txt`, `HL_magnitudes.txt`, `HL_center_points.txt`, `avg_wave_speeds.txt`;
- generated q-k point CSV files;
- generated trajectory, flow-density, and speed-variation figures;

New artifacts should be written under `results/` or another explicit output directory. The `.gitignore` file ignores those generated outputs by default.

## File Naming

Method names in generated output use:

- `Analytical` for analytical validation;
- `Benchmark` for the parallelogram method;
- `Rectangular` for the traditional rectangular baseline;
- `Optimization` for the proposed binary integer programming method.
