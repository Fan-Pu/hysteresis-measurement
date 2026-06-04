import math
from collections import defaultdict

import numpy as np

import CommonDefines
from PlotHelper import quadrilateral_area


class Benchmark:
    def __init__(self, vehicles, vehicle_id_sequence):
        self.vehicles = vehicles
        self.vehicle_id_sequence = vehicle_id_sequence
        self.selected_links = []  # value: (k,i,j)
        temp_infos = {}
        first_v_id = vehicle_id_sequence[0]
        # set the start and end points indices for each vehicle
        slots_ranges, self.vehicle_speeds, self.vehicle_slots, self.vehicle_pos = (
            CommonDefines.get_slot_ranges_by_w(vehicles, vehicle_id_sequence))  # [start_i, end_i]
        r_slots_ranges, _, _, _ = (CommonDefines.get_slot_ranges_by_reaction_time(vehicles, vehicle_id_sequence))
        slots_ranges[first_v_id] = (max(slots_ranges[first_v_id][0], r_slots_ranges[first_v_id][0]),
                                    min(slots_ranges[first_v_id][1], r_slots_ranges[first_v_id][1]))

        i1 = slots_ranges[first_v_id][0]
        i2 = i1 + CommonDefines.shift_interval
        self.region_edges = defaultdict(list)  # key: vehicle_id; value: [(t1,y1,t2,y2)] the ith edge on this vehicle
        while i2 <= slots_ranges[first_v_id][1]:
            self.region_edges[first_v_id].append((self.vehicle_slots[first_v_id][i1], self.vehicle_pos[first_v_id][i1],
                                                  self.vehicle_slots[first_v_id][i2], self.vehicle_pos[first_v_id][i2]))
            temp_infos[(first_v_id, len(self.region_edges[first_v_id]) - 1)] = (i1, i2)
            i1 += CommonDefines.vehicle_point_sample_interval
            i2 = i1 + CommonDefines.shift_interval

        for k1 in self.vehicle_id_sequence[:-1]:
            k2 = k1 + 1
            i3 = 0
            for idx in range(len(self.region_edges[k1])):
                t1, y1, t2, y2 = self.region_edges[k1][idx]
                (i1, i2) = temp_infos[(k1, idx)]
                i3, t3, y3 = self.find_w_match_infos(t1, y1, k2, i3)
                i4, t4, y4 = self.find_w_match_infos(t2, y2, k2, i3)
                if t3 == -1 or t4 == -1:
                    break
                self.region_edges[k2].append((t3, y3, t4, y4))
                temp_infos[(k2, len(self.region_edges[k2]) - 1)] = (i3, i4)
                self.selected_links.extend([(k1, i1, i3), (k1, i2, i4)])

    def find_w_match_infos(self, t1, y1, k2, start_t3_idx):
        i3 = t3 = y3 = -1
        for temp_i in range(start_t3_idx, len(self.vehicle_slots[k2]) - 1):
            t31, t32 = self.vehicle_slots[k2][temp_i], self.vehicle_slots[k2][temp_i + 1]
            y31, y32 = self.vehicle_pos[k2][temp_i], self.vehicle_pos[k2][temp_i + 1]
            A = CommonDefines.general_w
            B = -1
            C = y1 - A * t1
            d_l = A * t31 + B * y31 + C
            d_r = A * t32 + B * y32 + C
            dist_l = abs(A * t31 + B * y31 + C) / math.sqrt(A ** 2 + B ** 2)
            dist_r = abs(A * t32 + B * y32 + C) / math.sqrt(A ** 2 + B ** 2)
            # (t31, y31) is the point we founded
            if abs(dist_l) <= 0.001:
                i3, t3, y3 = temp_i, t31, y31
                break
            # (t32, y32) is the point we founded
            elif abs(dist_r) <= 0.001:
                i3, t3, y3 = temp_i + 1, t32, y32
                break
            # in the middle
            elif d_l * d_r < 0:
                ratio = dist_l / (dist_l + dist_r)
                i3 = temp_i + 1 if ratio > 0.5 else temp_i
                t3, y3 = t31 + ratio * (t32 - t31), y31 + ratio * (y32 - y31)
                break
        return i3, t3, y3

    def gen_qk_points(self):
        self.qk_points = []  # value: (q,k)
        self.q_list = []
        self.k_list = []
        self.region_speeds = []  # id: region_id, value: [speeds]
        self.wave_speeds = {}  # id: v_i - v_i+1
        self.reaction_times = []
        # each region_id has a q-k point
        min_length = min(len(edge) for edge in self.region_edges.values())
        self.region_edges = {key: edge[:min_length] for key, edge in self.region_edges.items()}
        for region_id in range(len(self.region_edges[0])):
            region_speeds = []
            t_list = []
            y_list = []
            area_list = []
            for k1 in range(len(self.vehicle_id_sequence) - 1):
                k2 = k1 + 1
                if (k1, k2) not in self.wave_speeds.keys():
                    self.wave_speeds[(k1, k2)] = []
                t1, y1, t2, y2 = self.region_edges[k1][region_id]
                p1_idx = int(round((t1 - self.vehicle_slots[k1][0]) / CommonDefines.step_size))
                p2_idx = int(round((t2 - self.vehicle_slots[k1][0]) / CommonDefines.step_size))

                t3, y3, t4, y4 = self.region_edges[k2][region_id]
                p3_idx = int(round((t3 - self.vehicle_slots[k2][0]) / CommonDefines.step_size))
                p4_idx = int(round((t4 - self.vehicle_slots[k2][0]) / CommonDefines.step_size))
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

    def check_valid(self):
        valid = True
        for k in range(0, len(self.region_edges)):
            if len(self.region_edges[k]) == 0:
                valid = False
                break
        return valid
