import cmath
from collections import defaultdict
import CommonDefines
from CommonDefines import *
from PlotHelper import quadrilateral_area
from Vehicle import *

class AnalyticalMethod:
    def __init__(self):
        self.vehicles = {}  # key: vehicle id; value: vehicle
        self.vehicles4lane = defaultdict(list)  # key: lane id; value: vehicle id
        self.objective_value = 0

        self.GS = complex(ks, kv * omega) / complex(-omega ** 2 + ks, -omega ** 3 * phi + omega * (kv + ks * tau))
        self.GS_magnitude = abs(self.GS)
        self.GS_phase = cmath.phase(self.GS)
        self.reaction_time = -self.GS_phase / omega
        # update newell's general wave speed (m/s)
        CommonDefines.general_w = -(xe / self.reaction_time) * 0.2778
        CommonDefines.speed_limit = 24.5872  # m/s
        print(f"analytical reaction time: {self.reaction_time} s")
        self.gen_plantoon_points()
        self.gen_measurement_regions()

    def gen_plantoon_points(self):
        t_list = np.arange(0, time_length * step_size, step_size)
        self.max_speed = 0
        # generate points for the first vehicle
        vehicle_id = 0
        self.lane_id = lane_id = 1
        if vehicle_id == 0:
            sin_wave = np.sin(omega * t_list + phase_shift)
            cos_wave = np.cos(omega * t_list + phase_shift)
            vehicle = Vehicle(vehicle_id, lane_id)
            positions_nomi = []  # m
            positions_osci = []
            positions = []
            speeds = []  # m/s
            accelerations = []  # m/s^2
            for i in range(len(t_list)):
                t = t_list[i]
                positions_nomi.append(ve * t)
                positions_osci.append(a * sin_wave[i])
                positions.append(positions_nomi[-1] + positions_osci[-1])
                speeds.append(ve + a * omega * cos_wave[i])
                self.max_speed = max(self.max_speed, speeds[-1])
                accelerations.append(-a * omega ** 2 * sin_wave[i])
            vehicle.Frame_IDs = t_list
            vehicle.Local_Y_nominal = positions_nomi
            vehicle.Local_Y_oscillation = positions_osci
            vehicle.Local_Y = positions
            vehicle.Mean_Speed = speeds
            vehicle.Accelerations = accelerations
            self.vehicles[vehicle_id] = vehicle
            self.vehicles4lane[lane_id].append(vehicle_id)
            vehicle_id += 1
        # add following vehicles
        while (vehicle_id < platoon_size):
            sin_wave = np.sin(omega * t_list + phase_shift + vehicle_id * self.GS_phase)
            cos_wave = np.cos(omega * t_list + phase_shift + vehicle_id * self.GS_phase)

            vehicle = Vehicle(vehicle_id, lane_id)
            positions_nomi = []  # m
            positions_osci = []
            positions = []
            speeds = []  # m/s
            accelerations = []  # m/s^2
            for i in range(len(t_list)):
                t = t_list[i]
                positions_nomi.append(ve * t - vehicle_id * xe)
                positions_osci.append(a * self.GS_magnitude ** vehicle_id * sin_wave[i])
                positions.append(positions_nomi[-1] + positions_osci[-1])
                speeds.append(ve + a * omega * self.GS_magnitude ** vehicle_id * cos_wave[i])
                self.max_speed = max(self.max_speed, speeds[-1])
                accelerations.append(-a * omega ** 2 * self.GS_magnitude ** vehicle_id * sin_wave[i])
            vehicle.Frame_IDs = t_list
            vehicle.Local_Y_nominal = positions_nomi
            vehicle.Local_Y_oscillation = positions_osci
            vehicle.Local_Y = positions
            vehicle.Mean_Speed = speeds
            vehicle.Accelerations = accelerations
            self.vehicles[vehicle_id] = vehicle
            self.vehicles4lane[lane_id].append(vehicle_id)
            vehicle_id += 1
        self.vehicle_id_sequence = self.vehicles4lane[lane_id]

    def gen_measurement_regions(self):
        vehicles = self.vehicles
        vehicle_id_sequence = self.vehicles4lane[self.lane_id]
        first_v_id = vehicle_id_sequence[0]
        # set the start and end points indices for each vehicle
        slots_ranges, self.vehicle_speeds, self.vehicle_slots, self.vehicle_pos = (
            get_slot_ranges_by_w(vehicles, vehicle_id_sequence))  # [start_i, end_i]
        r_slots_ranges, _, _, _ = get_slot_ranges_by_reaction_time(vehicles, vehicle_id_sequence)
        slots_ranges[first_v_id] = (max(slots_ranges[first_v_id][0], r_slots_ranges[first_v_id][0]),
                                    min(slots_ranges[first_v_id][1], r_slots_ranges[first_v_id][1]))

        i1 = slots_ranges[first_v_id][0]
        i2 = i1 + shift_interval
        self.region_edges = defaultdict(list)  # key: vehicle_id; value: [(t1,y1,t2,y2)] the ith edge on this vehicle
        while i2 <= slots_ranges[first_v_id][1]:
            self.region_edges[first_v_id].append((self.vehicle_slots[first_v_id][i1], self.vehicle_pos[first_v_id][i1],
                                                  self.vehicle_slots[first_v_id][i2], self.vehicle_pos[first_v_id][i2]))
            i1 += vehicle_point_sample_interval
            i2 = i1 + shift_interval

        for k1 in self.vehicle_id_sequence[:-1]:
            k2 = k1 + 1
            i3 = 0
            for idx in range(len(self.region_edges[k1])):
                t1, y1, t2, y2 = self.region_edges[k1][idx]
                i3, t3, y3, speed3 = self.find_reaction_time_match_infos(t1, y1, k2, i3)
                i4, t4, y4, speed4 = self.find_reaction_time_match_infos(t2, y2, k2, i3 + 1)
                if t3 == -1 or t4 == -1:
                    break
                self.region_edges[k2].append((t3, y3, t4, y4))

    def find_reaction_time_match_infos(self, t1, y1, k2, start_t3_idx):
        i3 = t3 = y3 = speed = -1
        for temp_i in range(start_t3_idx, len(self.vehicle_slots[k2]) - 1):
            t31, t32 = self.vehicle_slots[k2][temp_i], self.vehicle_slots[k2][temp_i + 1]
            y31, y32 = self.vehicle_pos[k2][temp_i], self.vehicle_pos[k2][temp_i + 1]
            react_t31 = t31 - t1
            react_t32 = t32 - t1
            # (t31, y31) is the point we founded
            if abs(react_t31 - self.reaction_time) < t_tolerance:
                i3, t3, y3 = temp_i, t31, y31
                speed = self.vehicle_speeds[k2][i3]
                break
            # (t32, y32) is the point we founded
            elif abs(react_t32 - self.reaction_time) < t_tolerance:
                i3, t3, y3 = temp_i + 1, t32, y32
                speed = self.vehicle_speeds[k2][i3]
                break
            # interpolation
            elif react_t31 < self.reaction_time < react_t32:
                ratio = (self.reaction_time - react_t31) / (react_t32 - react_t31)
                i3, t3, y3 = temp_i, t31 + ratio * (t32 - t31), y31 + ratio * (y32 - y31)
                speed_lb = self.vehicle_speeds[k2][temp_i]
                speed_ub = self.vehicle_speeds[k2][temp_i + 1]
                speed = speed_lb + ratio * (speed_ub - speed_lb)
                break
        return i3, t3, y3, speed

    def gen_qk_points(self):
        self.qk_points = []  # value: (q,k)
        self.q_list = []
        self.k_list = []
        self.region_speeds = []  # id: region_id, value: [speeds]
        self.wave_speeds = defaultdict(list)  # id: v_i - v_i+1
        self.reaction_times = []
        # each region_id has a q-k point
        for region_id in range(len(self.region_edges[0])):
            region_speeds = []
            t_list = []
            y_list = []
            area_list = []
            for k1 in range(len(self.vehicle_id_sequence) - 1):
                k2 = k1 + 1
                t1, y1, t2, y2 = self.region_edges[k1][region_id]
                p1_idx = int(round((t1 - self.vehicle_slots[k1][0]) / step_size))
                p2_idx = int(round((t2 - self.vehicle_slots[k1][0]) / step_size))

                t3, y3, t4, y4 = self.region_edges[k2][region_id]
                p3_idx = int(round((t3 - self.vehicle_slots[k2][0]) / step_size))
                p4_idx = int(round((t4 - self.vehicle_slots[k2][0]) / step_size))
                area_list.append(quadrilateral_area(t1, y1, t3, y3, t4, y4, t2, y2))
                t_list.append(t4 - t3)
                y_list.append(y4 - y3)

                self.reaction_times.append(0.5 * (t3 - t1 + t4 - t2))

                self.wave_speeds[(k1, k2)].append([(y1 - y3) / (t1 - t3), (y2 - y4) / (t2 - t4)])

                if k1 == 0:
                    t_list.append(t2 - t1)
                    y_list.append(y2 - y1)
                    for i in range(p1_idx, p2_idx + 1):
                        region_speeds.append(self.vehicle_speeds[k1][i])

                for i in range(p3_idx, p4_idx + 1):
                    region_speeds.append(self.vehicle_speeds[k2][i])

            plantoon_size = len(self.vehicle_slots)
            total_area = sum(area_list)
            total_area *= ((plantoon_size) / (plantoon_size - 1))
            q = sum(y_list) / total_area  # veh/s
            k = sum(t_list) / total_area  # veh/m
            self.q_list.append(q * 3600)
            self.k_list.append(k * 1000)
            self.qk_points.append((q, k))
            self.region_speeds.append(region_speeds)

        self.center_point = (np.mean(self.k_list), np.mean(self.q_list))
        groups = [i for i in range(len(self.region_speeds))]
        means = [np.mean(speeds) for speeds in self.region_speeds]
        std_dev = [np.std(speeds) for speeds in self.region_speeds]
        # return speed info
        return [groups, means, std_dev, self.wave_speeds]
