from collections import defaultdict

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoMinorLocator

from CommonDefines import *

plt.rcParams['font.family'] = 'Times New Roman'

measure_line_width = 1
qk_line_width = 1
measure_vertical_alpha = 1

title_font_size = 16
axis_font_size = 14

dpi = 300
traj_fig_size = [3000, 1650]
qk_fig_size = [2000, 1700]

scalar_ms2kmh = 3.6

# for q-k plot
qk_krange = (40, 120)
qk_qrange = (350, 1000)


def PlotVehicles(vehicles, max_speed, title, enable_color=True):
    # plot trajectories
    temp_fig = plt.figure()
    plt.title(title)

    # colormap
    hsv = plt.cm.get_cmap('hsv')
    new_colors = hsv(np.linspace(0.0, 1 / 3, 256))
    new_cmap = mcolors.LinearSegmentedColormap.from_list('Red_to_Green', new_colors)

    k = 0
    norm = Normalize(0, max_speed * scalar_ms2kmh)
    for vehicle_id, vehicle in vehicles.items():
        k += 1
        x = []
        y = []
        z = []
        for i in range(len(vehicle.Frame_IDs)):
            x_value = vehicle.Frame_IDs[i]
            x.append(x_value)
            y.append(vehicle.Local_Y[i])
            z.append(vehicle.Mean_Speed[i] * scalar_ms2kmh)
        if enable_color:
            sc = plt.scatter(x, y, c=z, cmap=new_cmap, norm=norm, s=5)
        else:
            sc = plt.scatter(x, y, color='0.8', s=5)
    if enable_color:
        cb = plt.colorbar(sc)
        cb.ax.tick_params(labelsize=axis_font_size)
        cb.set_label('Speed (km/h)', fontsize=title_font_size)
    plt.xlabel("Time (s)", fontsize=title_font_size)
    plt.ylabel("Space (m)", fontsize=title_font_size)
    plt.tick_params(axis='both', which='major', labelsize=axis_font_size)
    plt.tick_params(axis='both', which='major', length=8)
    plt.tick_params(axis='both', which='minor', length=4)
    plt.tick_params(axis='both', which='both', direction='out')
    ax = plt.gca()
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    plt.tight_layout()

    fig_list.append((temp_fig, title))

    return norm


def PlotMeasurmentLines4Solver(vehicles, solver=None):
    # plot measurement lines
    measurement_lines = defaultdict(list)  # key: (k1,k2), value:[(i,j)]
    if solver is not None:
        for (k1, i, k2, j), value in solver.var_vals.items():
            if abs(value - 1) <= 0.001:
                x1 = solver.vehicle_slots[k1][i]
                x2 = solver.vehicle_slots[k2][j]
                y1 = solver.vehicle_pos[k1][i]
                y2 = solver.vehicle_pos[k2][j]
                plt.plot([x1, x2], [y1, y2], color='0.7', linewidth=measure_line_width)
                measurement_lines[(k1, k2)].append((i, j))
        for key in measurement_lines:
            measurement_lines[key].sort(key=lambda x: x[0])
    # fill the region
    for k1 in solver.vehicle_id_sequence[:-1]:
        k2 = k1 + 1
        edges = measurement_lines[(k1, k2)]
        i = 0
        while i <= len(edges) - 2:
            (i1, j1) = edges[i]
            (i2, j2) = edges[i + 1]
            plt.fill([solver.vehicle_slots[k1][i1], solver.vehicle_slots[k2][j1],
                      solver.vehicle_slots[k2][j2], solver.vehicle_slots[k1][i2]],
                     [solver.vehicle_pos[k1][i1], solver.vehicle_pos[k2][j1],
                      solver.vehicle_pos[k2][j2], solver.vehicle_pos[k1][i2]], color='lightgrey', alpha=0.4,
                     edgecolor='none')
            i += 2
    # plot vertical lines
    for (k1, k2), lines in measurement_lines.items():
        if k1 == solver.vehicle_id_sequence[0]:
            for m in range(0, len(lines) - 1, 2):
                (i1, j1) = lines[m]
                (i2, j2) = lines[m + 1]
                # line from i1 to i2
                x1 = solver.vehicle_slots[k1][i1]
                x2 = solver.vehicle_slots[k1][i2]
                y1 = solver.vehicle_pos[k1][i1]
                y2 = solver.vehicle_pos[k1][i2]
                plt.plot([x1, x2], [y1, y2], color='0.7', linewidth=measure_line_width, alpha=measure_vertical_alpha)
        elif k2 == solver.vehicle_id_sequence[-1]:
            for m in range(0, len(lines) - 1, 2):
                (i1, j1) = lines[m]
                (i2, j2) = lines[m + 1]
                # line from j1 to j2
                x1 = solver.vehicle_slots[k2][j1]
                x2 = solver.vehicle_slots[k2][j2]
                y1 = solver.vehicle_pos[k2][j1]
                y2 = solver.vehicle_pos[k2][j2]
                plt.plot([x1, x2], [y1, y2], color='0.7', linewidth=measure_line_width, alpha=measure_vertical_alpha)
    solver.measurement_lines = measurement_lines


def distance(x1, y1, x2, y2):
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5


def triangle_area(x1, y1, x2, y2, x3, y3):
    a = distance(x1, y1, x2, y2)
    b = distance(x2, y2, x3, y3)
    c = distance(x3, y3, x1, y1)
    s = (a + b + c) / 2
    return (s * (s - a) * (s - b) * (s - c)) ** 0.5


def quadrilateral_area(x1, y1, x2, y2, x3, y3, x4, y4):
    '''clock-wise from p1'''
    area1 = triangle_area(x1, y1, x2, y2, x3, y3)
    area2 = triangle_area(x1, y1, x3, y3, x4, y4)
    return area1 + area2


def PlotMeasurmentLines(method=None, enable_color=False):
    if method is not None:
        objective_value = 0
        min_region_count = min(len(edge) for edge in method.region_edges.values())
        for key, region_edge_list in method.region_edges.items():
            if len(region_edge_list) > min_region_count:
                method.region_edges[key] = region_edge_list[:min_region_count]

        l_width = measure_line_width if not enable_color else 2 * measure_line_width
        if enable_color:
            wave_speeds = GetWaveSpeed(method)
            wave_speeds_abs = np.abs(wave_speeds)
            hsv = plt.cm.get_cmap('hsv')
            new_colors = hsv(np.linspace(0.0, 1 / 3, 256))
            new_cmap = mcolors.LinearSegmentedColormap.from_list('Red_to_Green', new_colors)
            norm_abs = Normalize(min(wave_speeds_abs), max(wave_speeds_abs))
            norm_neg = Normalize(min(wave_speeds), max(wave_speeds))

        last_v_id = method.vehicle_id_sequence[-1]
        for region_id in range(len(method.region_edges[last_v_id])):
            for k1 in method.vehicle_id_sequence[::-1][:-1]:
                k2 = k1 - 1
                t1, y1, t2, y2 = method.region_edges[k1][region_id]
                t3, y3, t4, y4 = method.region_edges[k2][region_id]
                if enable_color:
                    wave_speed_normalized = norm_abs(np.abs((y1 - y3) / (t1 - t3)))
                    plt.plot([t1, t3], [y1, y3], color=new_cmap(wave_speed_normalized), linewidth=l_width)
                    wave_speed_normalized = norm_abs(np.abs((y2 - y4) / (t2 - t4)))
                    plt.plot([t2, t4], [y2, y4], color=new_cmap(wave_speed_normalized), linewidth=l_width)
                else:
                    plt.plot([t1, t3], [y1, y3], color='0.7', linewidth=l_width)
                    plt.plot([t2, t4], [y2, y4], color='0.7', linewidth=l_width)
                plt.fill([t1, t3, t4, t2], [y1, y3, y4, y2], color='lightgrey', alpha=0.4, edgecolor='none')
                objective_value += 0
                # plot vertical lines
                v_color = '0.7' if not enable_color else '0.5'
                if not enable_color:
                    if k1 == last_v_id:
                        plt.plot([t1, t2], [y1, y2], color=v_color, linewidth=l_width,
                                 alpha=measure_vertical_alpha)
                    elif k2 == method.vehicle_id_sequence[0]:
                        plt.plot([t3, t4], [y3, y4], color=v_color, linewidth=l_width,
                                 alpha=measure_vertical_alpha)

        if enable_color:
            cb = plt.colorbar(mappable=plt.cm.ScalarMappable(norm=norm_abs, cmap=new_cmap), ax=plt.gca())
            cb.ax.tick_params(labelsize=axis_font_size)
            cb.set_label('Wave speed (m/s)', fontsize=title_font_size)


def PlotMeasurmentLines4RectangleMethod(rectangle_method=None):
    if rectangle_method is not None:
        min_region_count = min(map(len, rectangle_method.region_edges.values()))
        for key, region_edge_list in rectangle_method.region_edges.items():
            if len(region_edge_list) > min_region_count:
                rectangle_method.region_edges[key] = region_edge_list[:min_region_count]
        for region_id in range(len(rectangle_method.region_edges[0])):
            last_veh = rectangle_method.vehicle_id_sequence[-1]
            t1, y1, t2, y2 = rectangle_method.region_edges[0][region_id]
            t3, y3, t4, y4 = rectangle_method.region_edges[last_veh][region_id]
            # horizontal
            plt.plot([t1, t2], [y2, y2], color='0.7', linewidth=measure_line_width)
            plt.plot([t3, t4], [y3, y3], color='0.7', linewidth=measure_line_width)
            # vertical
            plt.plot([t1, t3], [y2, y3], color='0.7', linewidth=measure_line_width,
                     alpha=measure_vertical_alpha)
            plt.plot([t2, t4], [y2, y3], color='0.7', linewidth=measure_line_width,
                     alpha=measure_vertical_alpha)
            plt.fill([t1, t3, t4, t2], [y2, y3, y3, y2], color='lightgrey', alpha=0.4, edgecolor='none')


def PlotSpeedsInfo(info_list, color_list):
    '''measure standard errors of speeds in each mesurement region'''
    min_size = 9999
    # determin min_size
    for group_name, info in info_list:
        min_size = min(min_size, len(info[0]))
    # plot figue
    width = 8
    height = width / 1.8
    temp_fig = plt.figure(figsize=(width, height))
    title_name = 'speed standard error'
    ax = plt.gca()

    groups = None
    n_groups = min_size
    bar_width = 5
    n_bars = len(info_list)
    group_spacing = 1.5 * bar_width

    # x coordinates for groups
    index = np.arange(n_groups) * (n_bars * bar_width + group_spacing)

    for k in range(n_bars):
        group_name, info = info_list[k]
        # resize the info list
        for i in range(len(info)):
            size_diff = len(info[i]) - min_size
            if size_diff > 0:
                if i == 0:
                    info[i] = info[i][:-size_diff]
                else:
                    info[i] = info[i][size_diff:]
        groups = info[0]
        current_offset = k * bar_width
        ax.bar(index + current_offset, info[2], bar_width, label=group_name, color=color_list[k])

    # Adding labels, title, and customizations for the first axsis
    ax.set_xlabel('Index of measurement region', fontsize=title_font_size)
    ax.set_ylabel('Standard error of speeds (m/s)', fontsize=title_font_size)
    ax.set_xticks(index + (bar_width * n_bars) / 2 - bar_width / 2)
    ax.set_xticklabels(groups)
    ax.legend(fontsize=axis_font_size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax.tick_params(axis='both', which='major', length=8)
    ax.tick_params(axis='both', which='minor', length=4)
    plt.tick_params(axis='both', which='both', direction='out')
    plt.xticks(rotation=45)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    plt.tight_layout()
    fig_list.append((temp_fig, title_name))


def PlotReactionTimeInfo(info_list, color_list, analytical_reaction_time):
    '''plot cross product based on q-k points using right-hand rule'''
    x_size = np.min([len(info_list[i][-1]) for i in range(len(info_list))])
    plt.figure()
    for i in range(len(info_list)):
        group_name, reaction_times = info_list[i]
        plt.plot(range(x_size), reaction_times[0:x_size], label=group_name, color=color_list[i])

    # plot analytical reaction time
    plt.plot([0, x_size - 1], [analytical_reaction_time, analytical_reaction_time], color='black', linestyle='--')
    # plot figue
    ax = plt.gca()
    plt.title("Reaction time differences")
    plt.legend()
    plt.ylabel('Reaction time (s)')
    plt.xlabel('Index of reaction time measurement')
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    plt.tick_params(axis='both', which='both', direction='out')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def PlotCrossProduct(center_point, qk_points, group_name):
    '''plot cross product based on q-k points using right-hand rule'''
    y_values = []
    center_point = np.array(center_point)
    magnitude = 0
    for i in range(len(qk_points) - 1):
        point_i = np.array([qk_points[i][1] * 1000, qk_points[i][0] * 3600])
        point_i1 = np.array([qk_points[i + 1][1] * 1000, qk_points[i + 1][0] * 3600])
        product = np.cross(point_i - center_point, point_i1 - center_point)
        y_values.append(product)
        magnitude += abs(product / 2)
    return magnitude


def PlotQKAlltogether(qk_points_list):
    temp_fig = plt.figure()
    title_name = 'Flow-density diagram'
    plt.title(title_name)
    avg_wave_speeds = []
    for group_name, center_point, qk_points, color in qk_points_list:
        q_list = [point[0] * 3600 for point in qk_points]
        k_list = [point[-1] * 1000 for point in qk_points]
        q_min, q_max = np.min(q_list), np.max(q_list)
        k_min, k_max = np.min(k_list), np.max(k_list)
        avg_wave_speeds.append((group_name, (q_min - q_max) / (k_max - k_min)))
        plt.scatter(k_list[1:-1], q_list[1:-1], color=color, s=5)
        # start pint
        plt.scatter(k_list[0], q_list[0], color=color, facecolors='none', s=50, marker='P')
        # end point
        plt.scatter(k_list[-1], q_list[-1], color=color, facecolors='none', s=50, marker='s')
        # center point
        plt.scatter(center_point[0], center_point[1], color=color, facecolors='none', s=50)
        plt.plot(k_list, q_list, linewidth=qk_line_width, color=color, label=group_name)
    plt.xlabel("Density (veh/km)", fontsize=title_font_size)
    plt.ylabel("Flow (veh/hr)", fontsize=title_font_size)
    plt.tick_params(axis='both', which='major', labelsize=axis_font_size)
    plt.tick_params(axis='both', which='major', length=8)
    plt.tick_params(axis='both', which='minor', length=4)
    plt.tick_params(axis='both', which='both', direction='out')
    plt.legend(fontsize=axis_font_size)
    ax = plt.gca()
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig_list.append((temp_fig, title_name))
    return avg_wave_speeds


def GetWaveSpeed(method):
    wave_speeds = []
    last_v_id = method.vehicle_id_sequence[-1]
    for region_id in range(len(method.region_edges[last_v_id])):
        for k1 in method.vehicle_id_sequence[::-1][:-1]:
            k2 = k1 - 1
            t1, y1, t2, y2 = method.region_edges[k1][region_id]
            t3, y3, t4, y4 = method.region_edges[k2][region_id]
            wave_speeds.append((y1 - y3) / (t1 - t3))
            wave_speeds.append((y2 - y4) / (t2 - t4))
    return wave_speeds


def DrawWaveHeatmap(method, group_name):
    PlotVehicles(method.vehicles, np.Infinity, f"wave speed {group_name}", False)
    PlotMeasurmentLines(method, True)
