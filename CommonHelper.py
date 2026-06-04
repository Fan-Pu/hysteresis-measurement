import numpy as np

import CommonDefines
from Benchmark import Benchmark


def find_best_wave_speed(wave_lb, wave_ub, step_size, vehicles, vehicles4lane):
    best_w = best_std_error = np.inf
    for w in np.arange(wave_lb, wave_ub, step_size):
        CommonDefines.general_w = w * CommonDefines.kmph2mps
        temp_benchmark = Benchmark(vehicles, vehicles4lane)
        if not temp_benchmark.check_valid():
            continue
        temp_speeds_info = temp_benchmark.gen_qk_points()
        temp_std_error = np.mean(temp_speeds_info[2])
        if temp_std_error < best_std_error:
            best_w = w  # km/h
            best_std_error = temp_std_error
    return best_w, best_std_error
