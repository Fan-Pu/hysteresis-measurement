import argparse
import csv
import math
import time
from pathlib import Path


DEFAULT_DATASETS = {
    "analytical": "single_frequency",
    "ngsim": "Platoon_Position_3",
    "openacc": "dynamic_part17",
}

METHOD_ALIASES = {
    "optimization": "optimization",
    "solver": "optimization",
    "analytical": "analytical",
    "benchmark": "benchmark",
    "rectangle": "rectangle",
    "rectangular": "rectangle",
}

DEFAULT_METHODS = {
    "analytical": ["analytical", "optimization"],
    "ngsim": ["benchmark", "optimization"],
    "openacc": ["benchmark", "optimization"],
}


def load_runtime_modules():
    global plt, np
    global CommonDefines, CommonHelper, ReadInput
    global AnalyticalMethod, Benchmark, RectangleMethod
    global PlotCrossProduct, PlotMeasurmentLines, PlotMeasurmentLines4RectangleMethod
    global PlotMeasurmentLines4Solver, PlotQKAlltogether, PlotReactionTimeInfo
    global PlotSpeedsInfo, PlotVehicles, DrawWaveHeatmap

    import matplotlib.pyplot as plt
    import numpy as np

    import CommonDefines
    import CommonHelper
    import ReadInput
    from AnalyticalMethod import AnalyticalMethod
    from Benchmark import Benchmark
    from PlotHelper import (
        DrawWaveHeatmap,
        PlotCrossProduct,
        PlotMeasurmentLines,
        PlotMeasurmentLines4RectangleMethod,
        PlotMeasurmentLines4Solver,
        PlotQKAlltogether,
        PlotReactionTimeInfo,
        PlotSpeedsInfo,
        PlotVehicles,
    )
    from RectangleMethod import RectangleMethod


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run the traffic hysteresis measurement experiments from "
            "'Optimal Measurement of Traffic Hysteresis under Traffic Oscillations'."
        )
    )
    parser.add_argument(
        "--scenario",
        choices=["analytical", "ngsim", "openacc"],
        default="analytical",
        help="Experiment family to run.",
    )
    parser.add_argument(
        "--dataset",
        help=(
            "Dataset or instance name. Examples: Platoon_Position_3 for NGSIM, "
            "dynamic_part17 for OpenACC. Ignored by the analytical scenario."
        ),
    )
    parser.add_argument(
        "--methods",
        help=(
            "Comma-separated methods to run: analytical, benchmark, rectangle, "
            "optimization. Defaults depend on the scenario."
        ),
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for generated figures, CSV files, text summaries, and solver logs.",
    )
    parser.add_argument(
        "--show-plots",
        dest="show_plots",
        action="store_true",
        help="Display Matplotlib windows after saving generated figures.",
    )
    parser.add_argument(
        "--no-show-plots",
        dest="show_plots",
        action="store_false",
        help="Save figures without opening Matplotlib windows.",
    )
    parser.add_argument(
        "--draw-wave-speed-heat-map",
        action="store_true",
        help="Draw the analytical wave-speed heat map in addition to standard plots.",
    )
    parser.set_defaults(show_plots=False)
    return parser.parse_args()


def normalize_methods(raw_methods, scenario):
    if raw_methods is None:
        return DEFAULT_METHODS[scenario]
    methods = []
    for item in raw_methods.split(","):
        method = item.strip().lower()
        if method == "":
            continue
        if method not in METHOD_ALIASES:
            raise ValueError(f"Unknown method: {item}")
        methods.append(METHOD_ALIASES[method])
    if len(methods) == 0:
        raise ValueError("At least one method must be selected.")
    if "analytical" in methods and scenario != "analytical":
        raise ValueError("The analytical method is only available for --scenario analytical.")
    return methods


def make_output_dir(scenario, dataset, output_dir):
    if output_dir is not None:
        path = Path(output_dir)
    elif scenario == "analytical":
        path = Path("results") / "analytical" / str(int(time.time()))
    else:
        path = Path("results") / scenario / dataset
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_time_span(csv_path, dataset):
    with open(csv_path, mode="r", encoding="utf-8") as file:
        csv_reader = csv.reader(file)
        next(csv_reader, None)
        for row in csv_reader:
            if row[0] == dataset:
                return row
    raise ValueError(f"Dataset {dataset} was not found in {csv_path}.")


def load_input_data(scenario, dataset):
    lane_id = 1
    if scenario == "analytical":
        analytical_method = AnalyticalMethod()
        return analytical_method.vehicles, analytical_method.vehicles4lane, analytical_method, lane_id

    if scenario == "ngsim":
        CommonDefines.reaction_time_ub = 3
        root_path = Path("Ngsim") / "Platoon_data"
        row = read_time_span(root_path / "time_span.csv", dataset)
        ReadInput.slot_id_lb = int(row[1])
        ReadInput.slot_id_ub = ReadInput.slot_id_lb + CommonDefines.read_point_num
        vehicles, vehicles4lane = ReadInput.readNgsim(root_path / f"{dataset}.csv")
        return vehicles, vehicles4lane, None, lane_id

    root_path = Path("ACC") / "ZalaZONE"
    row = read_time_span(root_path / "time_span.csv", dataset)
    ReadInput.slot_id_lb = int(row[1])
    ReadInput.slot_id_ub = ReadInput.slot_id_lb + int(row[2])
    vehicles, vehicles4lane = ReadInput.readACC(root_path / f"{dataset}.csv")
    CommonDefines.speed_limit = 27
    return vehicles, vehicles4lane, None, lane_id


def append_method_outputs(
    group_name,
    method,
    speeds_info,
    speed_color,
    qk_color,
    info_list,
    speed_color_list,
    qk_points_list,
    mean_speed_std,
    speed_std_infos,
    wave_speed_infos,
    hl_magnitudes,
    hl_center_points,
):
    center_point, qk_points = method.center_point, method.qk_points
    mean_speed_std.append((group_name, np.mean(speeds_info[2])))
    speed_std_infos.append((group_name, speeds_info[2]))
    if len(speeds_info) > 3:
        wave_speed_infos.append((group_name, speeds_info[3]))
    info_list.append((group_name, speeds_info))
    speed_color_list.append(speed_color)
    qk_points_list.append((group_name, center_point, qk_points, qk_color))
    method.magnitude = PlotCrossProduct(center_point, qk_points, group_name)
    hl_magnitudes.append((group_name, method.magnitude))
    hl_center_points.append((group_name, method.center_point))


def write_comparison_outputs(
    result_path,
    baseline_name,
    mean_speed_std,
    qk_points_list,
    hl_magnitudes,
    hl_center_points,
    avg_wave_speeds,
    speed_std_infos,
):
    compare_infos = {
        "avg wave speed": {},
        "cp distance": {},
        "loop magnitude": {},
        "standard error of speeds": {},
    }

    with open(result_path / "mean speed std.txt", "w", encoding="utf-8") as file:
        infos = compare_infos["standard error of speeds"]
        base_value = -1
        for group_name, std_val in mean_speed_std:
            file.write(f"{group_name},{std_val}\n")
            if group_name == baseline_name:
                base_value = std_val
                info = f"{group_name}, {std_val:.2f}, -"
            else:
                info = f"{group_name}, {std_val:.2f}, {(std_val - base_value) / base_value * 100:.2f}"
            infos[group_name] = info

    for group_name, center_point, qk_points, _ in qk_points_list:
        with open(result_path / f"qk_points_{group_name}.csv", "w", encoding="utf-8") as file:
            file.write("x,y\n")
            file.write(f"{center_point[0]},{center_point[-1]}\n")
            for qk_point in qk_points:
                file.write(f"{qk_point[-1] * 1000},{qk_point[0] * 3600}\n")

    with open(result_path / "HL_magnitudes.txt", "w", encoding="utf-8") as file:
        infos = compare_infos["loop magnitude"]
        base_value = -1
        for group_name, val in hl_magnitudes:
            file.write(f"{group_name},{val}\n")
            if group_name == baseline_name:
                base_value = val
                info = f"{group_name}, {val:.2f}, -"
            else:
                info = f"{group_name}, {val:.2f}, {(val - base_value) / base_value * 100:.2f}"
            infos[group_name] = info

    with open(result_path / "HL_center_points.txt", "w", encoding="utf-8") as file:
        infos = compare_infos["cp distance"]
        base_value = -1
        for group_name, val in hl_center_points:
            file.write(f"{group_name},{val}\n")
            if group_name == baseline_name:
                base_value = val
                info = f"{group_name}, {val[0]:.2f}, {val[1]:.2f}, -"
            else:
                dist = math.sqrt((val[0] - base_value[0]) ** 2 + (val[1] - base_value[1]) ** 2)
                info = f"{group_name}, {val[0]:.2f}, {val[1]:.2f}, {dist:.2f}"
            infos[group_name] = info

    with open(result_path / "avg_wave_speeds.txt", "w", encoding="utf-8") as file:
        infos = compare_infos["avg wave speed"]
        base_value = -1
        for group_name, val in avg_wave_speeds:
            file.write(f"{group_name},{val}\n")
            if group_name == baseline_name:
                base_value = val
                info = f"{group_name}, {val:.2f}, -"
            else:
                info = f"{group_name}, {val:.2f}, {(val - base_value) / base_value * 100:.2f}"
            infos[group_name] = info

    for group_name, std_speeds in speed_std_infos:
        with open(result_path / f"std_speeds_points_{group_name}.csv", "w", encoding="utf-8") as file:
            file.write("id,value\n")
            for idx, std_speed in enumerate(std_speeds):
                file.write(f"{idx},{std_speed}\n")

    with open(result_path / "compare_infos.txt", "w", encoding="utf-8") as file:
        for row_name, infos in compare_infos.items():
            file.write(f"{row_name}:\n")
            for _, text in infos.items():
                file.write(f"{text}\n")
            file.write("\n")


def save_figures(result_path):
    for fig, title_name in CommonDefines.fig_list:
        ax = fig.axes[0]
        ax.set_title("")
        fig.savefig(result_path / f"{title_name}.png", dpi=300)


def run_experiment(args):
    load_runtime_modules()

    dataset = args.dataset or DEFAULT_DATASETS[args.scenario]
    methods = normalize_methods(args.methods, args.scenario)
    result_path = make_output_dir(args.scenario, dataset, args.output_dir)
    CommonDefines.fig_list.clear()

    vehicles, vehicles4lane, analytical_method, lane_id = load_input_data(args.scenario, dataset)

    info_list = []
    speed_color_list = []
    reaction_time_color_list = []
    mean_speed_std = []
    speed_std_infos = []
    reaction_time_list = []
    wave_speed_infos = []
    qk_points_list = []
    hl_magnitudes = []
    hl_center_points = []

    if "analytical" in methods:
        group_name = "Analytical"
        analytical_reaction_time = analytical_method.reaction_time
        PlotVehicles(analytical_method.vehicles, CommonDefines.speed_limit, group_name)
        PlotMeasurmentLines(analytical_method)
        if args.draw_wave_speed_heat_map:
            DrawWaveHeatmap(analytical_method, group_name)
        analytical_speeds_info = analytical_method.gen_qk_points()
        append_method_outputs(
            group_name,
            analytical_method,
            analytical_speeds_info,
            "#A2A2A2",
            "#4169E1",
            info_list,
            speed_color_list,
            qk_points_list,
            mean_speed_std,
            speed_std_infos,
            wave_speed_infos,
            hl_magnitudes,
            hl_center_points,
        )
    else:
        analytical_reaction_time = None

    if "benchmark" in methods:
        group_name = "Benchmark"
        if args.scenario == "openacc":
            best_w, _ = CommonHelper.find_best_wave_speed(-30, 0, 0.1, vehicles, vehicles4lane[lane_id])
            print(f"Best estimated wave speed: {best_w} km/h")
            CommonDefines.general_w = best_w * CommonDefines.kmph2mps
        benchmark = Benchmark(vehicles, vehicles4lane[lane_id])
        PlotVehicles(vehicles, CommonDefines.speed_limit, group_name)
        PlotMeasurmentLines(benchmark)
        benchmark_speeds_info = benchmark.gen_qk_points()
        append_method_outputs(
            group_name,
            benchmark,
            benchmark_speeds_info,
            "#A4CEE1",
            "#2E8B57",
            info_list,
            speed_color_list,
            qk_points_list,
            mean_speed_std,
            speed_std_infos,
            wave_speed_infos,
            hl_magnitudes,
            hl_center_points,
        )
        reaction_time_color_list.append("#DC143C")
        if analytical_reaction_time is not None:
            reaction_time_list.append((group_name, benchmark.reaction_times))

    if "rectangle" in methods:
        group_name = "Rectangular"
        rectangle_method = RectangleMethod(vehicles, vehicles4lane[lane_id])
        PlotVehicles(vehicles, CommonDefines.speed_limit, group_name)
        PlotMeasurmentLines4RectangleMethod(rectangle_method)
        rectangle_speeds_info = rectangle_method.gen_qk_points()
        append_method_outputs(
            group_name,
            rectangle_method,
            rectangle_speeds_info,
            "#DADADA",
            "#808080",
            info_list,
            speed_color_list,
            qk_points_list,
            mean_speed_std,
            speed_std_infos,
            wave_speed_infos,
            hl_magnitudes,
            hl_center_points,
        )
        reaction_time_color_list.append("#FFD700")
        if analytical_reaction_time is not None:
            reaction_time_list.append((group_name, rectangle_method.reaction_times))

    if "optimization" in methods:
        from Solver import Solver

        group_name = "Optimization"
        solver = Solver(vehicles, vehicles4lane[lane_id], str(result_path) + "/")
        PlotVehicles(vehicles, CommonDefines.speed_limit, group_name)
        PlotMeasurmentLines4Solver(vehicles, solver)
        optimization_speeds_info = solver.gen_qk_points()
        append_method_outputs(
            group_name,
            solver,
            optimization_speeds_info,
            "#FD9999",
            "#DC143C",
            info_list,
            speed_color_list,
            qk_points_list,
            mean_speed_std,
            speed_std_infos,
            wave_speed_infos,
            hl_magnitudes,
            hl_center_points,
        )
        reaction_time_color_list.append("#2E8B57")
        if analytical_reaction_time is not None:
            reaction_time_list.append((group_name, solver.reaction_times))

    if len(info_list) > 0:
        PlotSpeedsInfo(info_list, speed_color_list)

    if len(qk_points_list) > 0:
        avg_wave_speeds = PlotQKAlltogether(qk_points_list)
    else:
        avg_wave_speeds = []

    if analytical_reaction_time is not None and len(reaction_time_list) > 0:
        PlotReactionTimeInfo(reaction_time_list, reaction_time_color_list, analytical_reaction_time)

    baseline_name = "Analytical" if "analytical" in methods else qk_points_list[0][0]
    write_comparison_outputs(
        result_path,
        baseline_name,
        mean_speed_std,
        qk_points_list,
        hl_magnitudes,
        hl_center_points,
        avg_wave_speeds,
        speed_std_infos,
    )
    save_figures(result_path)

    print(f"Results written to: {result_path}")
    if args.show_plots:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    run_experiment(parse_args())
