from collections import defaultdict

import CommonDefines
from CommonDefines import *
from PlotHelper import quadrilateral_area


class RectangleMethod:
    def __init__(self, vehicles, vehicle_id_sequence):
        self.vehicles = vehicles
        self.vehicle_id_sequence = vehicle_id_sequence
        first_v_id = vehicle_id_sequence[0]
        # set the start and end points indices for each vehicle
        slots_ranges, self.vehicle_speeds, self.vehicle_slots, self.vehicle_pos = (
            CommonDefines.get_slot_ranges_by_w(vehicles, vehicle_id_sequence))  # [start_i, end_i]
        r_slots_ranges, _, _, _ = (CommonDefines.get_slot_ranges_by_reaction_time(vehicles, vehicle_id_sequence))
        slots_ranges[first_v_id] = (max(slots_ranges[first_v_id][0], r_slots_ranges[first_v_id][0]),
                                    min(slots_ranges[first_v_id][1], r_slots_ranges[first_v_id][1]))

        temp_infos = {}

        i1 = slots_ranges[first_v_id][0]
        i2 = i1 + shift_interval
        self.region_edges = defaultdict(list)  # key: vehicle_id; value: [(t1,y1,t2,y2)] the ith edge on this vehicle
        while i2 <= slots_ranges[first_v_id][1]:
            self.region_edges[first_v_id].append((self.vehicle_slots[first_v_id][i1], self.vehicle_pos[first_v_id][i1],
                                                  self.vehicle_slots[first_v_id][i2], self.vehicle_pos[first_v_id][i2]))

            temp_infos[(first_v_id, len(self.region_edges[first_v_id]) - 1)] = (i1, i2)

            i1 += CommonDefines.vehicle_point_sample_interval
            i2 = i1 + shift_interval

        # for following vehicles
        for k1 in self.vehicle_id_sequence[:-1]:
            k2 = k1 + 1
            for idx in range(len(self.region_edges[k1])):
                t1, y1, t2, y2 = self.region_edges[k1][idx]
                (i1, i2) = temp_infos[(k1, idx)]
                i3 = i1
                t3, y3 = self.vehicle_slots[k2][i3], self.vehicle_pos[k2][i3]
                i4 = i2
                t4, y4 = self.vehicle_slots[k2][i4], self.vehicle_pos[k2][i4]
                temp_infos[(k2, idx)] = (i3, i4)
                self.region_edges[k2].append((t3, y3, t4, y4))

    def gen_qk_points(self):
        self.qk_points = []  # value: (q,k)
        self.q_list = []
        self.k_list = []
        self.region_speeds = []  # id: region_id, value: [speeds]
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
                p1_idx = int(round((t1 - self.vehicle_slots[k1][0]) / CommonDefines.step_size))
                p2_idx = int(round((t2 - self.vehicle_slots[k1][0]) / CommonDefines.step_size))

                t3, y3, t4, y4 = self.region_edges[k2][region_id]
                p3_idx = int(round((t3 - self.vehicle_slots[k2][0]) / CommonDefines.step_size))
                p4_idx = int(round((t4 - self.vehicle_slots[k2][0]) / CommonDefines.step_size))
                if k1 == self.vehicle_id_sequence[0]:
                    area_list.append(quadrilateral_area(t1, y2, t3, y3, t4, y4, t2, y2))
                elif k1 == self.vehicle_id_sequence[-2]:
                    area_list.append(quadrilateral_area(t1, y1, t3, y3, t4, y3, t2, y2))
                else:
                    area_list.append(quadrilateral_area(t1, y1, t3, y3, t4, y4, t2, y2))

                t_list.append(t4 - t3)
                y_list.append(y4 - y3)

                self.reaction_times.append(0.5 * (t3 - t1 + t4 - t2))

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
        return [groups, means, std_dev]
