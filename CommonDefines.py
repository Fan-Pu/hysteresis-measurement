import numpy as np

# general info
kmph2mps = 0.2778
general_w = -20 * kmph2mps  # m/s
speed_limit = 7.82  # m/s
max_speed = 0  # m/s
min_speed = 0  # m/s
w_tolerance = 0.001
small_value = 0.00001
t_tolerance = 0.01
read_point_num = 350
scale = 0.3048  # ft to m
step_size = 0.1
vf = 72 / 3.6  # km/h to m/s

# for optimization solver
reaction_time_lb = 0.7
reaction_time_ub = 2.0
outflow_sample_size = 5
vehicle_point_sample_interval = 20  # = \Delta\tau_1 * 10
shift_interval = 19  # = \Delta\tau_2 * 10
obj_speed_diff_penalty = 1

# for analytical method
time_length = 200  # how many slots
platoon_size = 4
# the oscillation part of the first vehicle
a = 15
omega = 0.16 * np.pi
phase_shift = 0 * np.pi
phi = 0.1  # time-lag for vehicle k to realize the desired acceleration
f = omega / (2 * np.pi)  # frequency (hz)
# other parameters
ks = 1  # feedback gain
kv = 1  # feedback gain
tau = 1.2  # s minimal time gap
ve = 10  # m/s equilibrium speed
standstill_s = 5  # 5m standstill spacing
xe = tau * ve + standstill_s  # equilibrium spacing

delta_xe = 0

fig_list = []


def get_slot_ranges_by_w(vehicles, vehicle_id_sequence):
    vehicle_speeds = []
    vehicle_slots = []
    vehicle_pos = []
    for vehicle_id in vehicle_id_sequence:
        vehicle = vehicles[vehicle_id]
        temp_speeds = []
        temp_slots = []
        temp_pos = []
        for idx in range(len(vehicle.Frame_IDs)):
            temp_speeds.append(vehicle.Mean_Speed[idx])
            temp_slots.append(vehicle.Frame_IDs[idx])
            temp_pos.append(vehicle.Local_Y[idx])
        vehicle_speeds.append(temp_speeds)
        vehicle_slots.append(temp_slots)
        vehicle_pos.append(temp_pos)

    # exclude some points
    temp_w_line_bottom = {str(round(vehicle_slots[0][0], 1)): vehicle_pos[0][0]}  # key: t; value: y
    temp_w_line_top = {str(round(vehicle_slots[-1][-1], 1)): vehicle_pos[-1][-1]}
    # for temp_w_line_bottom
    for t in np.arange(vehicle_slots[0][0] + step_size, vehicle_slots[-1][-1] + step_size, step_size):
        pre_t = t - step_size
        pre_pos = temp_w_line_bottom[str(round(pre_t, 1))]
        temp_w_line_bottom[str(round(t, 1))] = pre_pos + general_w * step_size
    # for temp_w_line_top
    for t in np.arange(vehicle_slots[-1][-1] - step_size, vehicle_slots[0][0], -step_size):
        pre_t = t + step_size
        pre_pos = temp_w_line_top[str(round(pre_t, 1))]
        key = str(round(t, 1))
        if key == '-0.0':
            key = '0.0'
        temp_w_line_top[key] = pre_pos - general_w * step_size

    # set the start and end points indices for each vehicle
    slots_ranges = []  # [start_i, end_i]
    for k in range(len(vehicle_id_sequence)):
        start_i = 0
        end_i = len(vehicle_slots[-1]) - 1
        if k != 0:
            for i in range(len(vehicle_slots[k])):
                key = str(round(vehicle_slots[k][i], 1))
                if vehicle_pos[k][i] >= temp_w_line_bottom[key]:
                    start_i = i
                    break

        slots_ranges.append((start_i, end_i))

    return slots_ranges, vehicle_speeds, vehicle_slots, vehicle_pos


def get_slot_ranges_by_reaction_time(vehicles, vehicle_id_sequence):
    vehicle_speeds = []
    vehicle_accelerations = []
    vehicle_slots = []
    vehicle_pos = []
    for vehicle_id in vehicle_id_sequence:
        vehicle = vehicles[vehicle_id]
        temp_speeds = []
        temp_accs = []
        temp_slots = []
        temp_pos = []
        for idx in range(len(vehicle.Frame_IDs)):
            temp_speeds.append(vehicle.Mean_Speed[idx])
            temp_slots.append(vehicle.Frame_IDs[idx])
            temp_pos.append(vehicle.Local_Y[idx])
            temp_accs.append(vehicle.Accelerations[idx])
        vehicle_speeds.append(temp_speeds)
        vehicle_slots.append(temp_slots)
        vehicle_pos.append(temp_pos)
        vehicle_accelerations.append(temp_accs)
    # set the start and end points indices for each vehicle
    slots_ranges = []  # [(start_i, end_i)]
    # for the start points
    start_indices = {0: 0}
    for k in range(1, len(vehicle_id_sequence)):
        last_vehicle_t = vehicle_slots[k - 1][start_indices[k - 1]]
        for i in range(len(vehicle_slots[k])):
            temp_t = vehicle_slots[k][i]
            if temp_t - last_vehicle_t >= reaction_time_lb:
                start_indices[k] = i
                break
    # for the end points
    end_indices = {vehicle_id_sequence[-1]: len(vehicle_slots[vehicle_id_sequence[-1]]) - 1}
    for k in range(len(vehicle_id_sequence) - 2, -1, -1):
        last_vehicle_t = vehicle_slots[k + 1][end_indices[k + 1]]
        for i in range(len(vehicle_slots[k]) - 1, 0, -1):
            temp_t = vehicle_slots[k][i]
            if last_vehicle_t - temp_t >= reaction_time_ub:
                end_indices[k] = i
                break
    for k in vehicle_id_sequence:
        slots_ranges.append((start_indices[k], end_indices[k]))

    return slots_ranges, vehicle_speeds, vehicle_slots, vehicle_pos
