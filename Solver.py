from gurobipy import *

import CommonDefines
from CommonDefines import *
from PlotHelper import quadrilateral_area
from collections import defaultdict


class Solver:
    def __init__(self, vehicles, vehicle_id_sequence, output_folder, benchmark=None):
        reaction_time_ub = CommonDefines.reaction_time_ub
        reaction_time_lb = CommonDefines.reaction_time_lb

        self.vehicles = vehicles
        self.vehicle_id_sequence = vehicle_id_sequence
        first_v_id = vehicle_id_sequence[0]
        slots_ranges, self.vehicle_speeds, self.vehicle_slots, self.vehicle_pos = (
            get_slot_ranges_by_reaction_time(vehicles, vehicle_id_sequence))
        w_slots_ranges, _, _, _ = (
            get_slot_ranges_by_w(vehicles, vehicle_id_sequence))
        slots_ranges[first_v_id] = (max(slots_ranges[first_v_id][0], w_slots_ranges[first_v_id][0]),
                                    min(slots_ranges[first_v_id][1], w_slots_ranges[first_v_id][1]))

        first_vehicle_slot_set = []
        i = slots_ranges[first_v_id][0]
        i_shift = i + shift_interval
        while i_shift <= slots_ranges[first_v_id][1]:
            first_vehicle_slot_set += [i, i_shift]
            i += vehicle_point_sample_interval
            i_shift = i + shift_interval

        self.slots_ranges = slots_ranges

        # construct link set
        self.inflows = defaultdict(list)  # key: (k,i), value:[j]
        self.outflows = defaultdict(list)  # key: (k,i), value:[j]
        # positive direction
        for k in self.vehicle_id_sequence[:-1]:
            # for k and k+1
            i_range = range(slots_ranges[k][0], slots_ranges[k][1] + 1) if k != first_v_id \
                else first_vehicle_slot_set
            for i in i_range:
                temp_outflows = []
                # this is the point at k+1
                pivot_j = max(0, slots_ranges[k + 1][0], int(i + reaction_time_lb / 0.1))
                for j in range(pivot_j, slots_ranges[k + 1][1] + 1):
                    if j < slots_ranges[k + 1][0]:
                        continue
                    i_t = self.vehicle_slots[k][i]
                    j_t = self.vehicle_slots[k + 1][j]

                    approved = False
                    # equals the lower bound
                    if math.isclose(j_t - i_t, reaction_time_lb):
                        approved = True
                    # equals the upper bound
                    elif math.isclose(j_t - i_t, reaction_time_ub):
                        approved = True
                    elif j_t - i_t > reaction_time_ub:
                        break
                    elif j_t - i_t < reaction_time_lb:
                        continue
                    else:
                        approved = True

                    if approved:
                        if k == first_v_id:
                            temp_outflows.append(j)
                        # only if the node has inflow links, it can have outflow links
                        elif len(self.inflows[(k, i)]) > 0:
                            temp_outflows.append(j)
                # add links
                for j in temp_outflows:
                    if j not in self.outflows[(k, i)]:
                        self.outflows[(k, i)].append(j)
                    if i not in self.inflows[(k + 1, j)]:
                        self.inflows[(k + 1, j)].append(i)

        # sample some points to reduce computational burden (except last the vehicle)
        for (k, i), links in self.outflows.items():
            if len(links) > outflow_sample_size and k != first_v_id:
                interval = len(links) / float(outflow_sample_size)
                sampled_indices = set([int(round(i * interval)) for i in range(outflow_sample_size)])
                remove_indices = sorted(
                    set(range(len(links))).difference(sampled_indices), reverse=True)
                for idx in remove_indices:
                    j = links[idx]
                    # remove outflow links
                    del links[idx]
                    # remove from inflow links
                    temp_id = self.inflows[(k + 1, j)].index(i)
                    del self.inflows[(k + 1, j)][temp_id]

        self.img_links = defaultdict()  # key: (k,i,j)  value: (i1,j1)
        self.img_inverse_links = defaultdict()  # key: (k,i,j)  value: (i1,j1)

        self.gen_img_links()

        self.delete_trivial_links()

        self.update_slot_ranges()

        # get conflict outflow links---key: (k,i,j), value:[(i2,j2)]. link (k,i,k+1,j) conflict with link (k,i2,k+1,j2)
        self.conflict_links = defaultdict(list)
        for k in self.vehicle_id_sequence[:-1]:
            k_point_size = len(self.vehicle_slots[k])
            for i in range(slots_ranges[k][0], slots_ranges[k][1]):
                for j in self.outflows[(k, i)]:
                    for i2 in range(i + 1, k_point_size):
                        conflict_elems = [(i2, j2) for j2 in self.outflows[(k, i2)] if j2 <= j]
                        if len(conflict_elems) == 0 and len(self.outflows[(k, i2)]) != 0:
                            break
                        else:
                            self.conflict_links[(k, i, j)] += conflict_elems

        # set model
        self.model = Model("main")
        self.vars = {}
        self.var_list = list()
        self.var_names = {}  # key:(k,i,k+1,j). value:index in var_list
        self.z_vars = {}
        # auxiliary variables
        self.vars_l = {}
        self.var_names_l = {}  # key:(k,i). value:index in var_list_l
        self.var_list_l = list()
        # generate variables
        for k in self.vehicle_id_sequence[:-1]:
            for i in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                for j in self.outflows[(k, i)]:
                    key = (k, i, k + 1, j)
                    x = self.model.addVar(vtype=GRB.BINARY, name=f"x[{k},{i},{k + 1},{j}]", lb=0, ub=1)
                    self.vars[key] = x
                    self.var_list.append(x)
                    self.var_names[key] = len(self.var_list) - 1
        for k in self.vehicle_id_sequence[1:]:
            for i in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                # auxiliary variables
                key = (k, i)
                x_l = self.model.addVar(vtype=GRB.BINARY, name=f"xl[{k},{i}]", lb=0, ub=1)
                self.vars_l[key] = x_l
                self.var_list_l.append(x_l)
                self.var_names_l[key] = len(self.var_list_l) - 1

        # set objective
        objective = 0
        # obj 1
        for (k, i, m, j), var in self.vars.items():
            objective += var * abs(self.vehicle_speeds[k][i] - self.vehicle_speeds[m][j]) * obj_speed_diff_penalty

        # add constraints---group 1
        outflow_count = 0
        cons_num_list = []
        for i in first_vehicle_slot_set:
            if len(self.outflows[(first_v_id, i)]) == 0:
                continue
            var_names = [(first_v_id, i, first_v_id + 1, j) for j in self.outflows[(first_v_id, i)]]
            lhs = quicksum([self.vars[var_name] for var_name in var_names])
            self.model.addConstr(lhs == 1, f"cons_1_{i}")
            outflow_count += 1
        # group 2 (only match once)
        id = 0
        for k in self.vehicle_id_sequence[1:]:
            for i in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                var_names = [(k - 1, j, k, i) for j in self.inflows[(k, i)]]
                if len(var_names) == 0:
                    continue
                lhs = quicksum([self.vars[var_name] for var_name in var_names])
                self.model.addConstr(lhs <= 1, f"cons_2_{id}")
                id += 1
        cons_num_list.append(id)
        # group 3 (flow conservation)
        id = 0
        for k in vehicle_id_sequence[1:-1]:
            for i in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                if len(self.outflows[(k, i)]) == 0 and len(self.inflows[(k, i)]) == 0:
                    continue
                var_names = [(k, i, k + 1, j) for j in self.outflows[(k, i)]]
                sum_outflow = quicksum([self.vars[var_name] for var_name in var_names])
                var_names = [(k - 1, j, k, i) for j in self.inflows[(k, i)]]
                sum_inflow = quicksum([self.vars[var_name] for var_name in var_names])
                self.model.addConstr(sum_outflow - sum_inflow == 0, f"cons_3_{id}")
                id += 1
        cons_num_list.append(id)
        # group 4 (no overlap)
        id = 0
        for k in self.vehicle_id_sequence[:-1]:
            for i in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                for j in self.outflows[(k, i)]:
                    M = len(self.conflict_links[(k, i, j)])
                    if M == 0:
                        continue
                    lhs = 0
                    for (i2, j2) in self.conflict_links[(k, i, j)]:
                        lhs += self.vars[(k, i2, k + 1, j2)]
                    cons = self.model.addConstr(lhs <= (1 - self.vars[(k, i, k + 1, j)]) * M, f"cons_4_{id}")
                    id += 1
        cons_num_list.append(id)
        # group 5 (for auxiliary variables)
        id = 0
        # part 1
        k = first_v_id + 1
        for j in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
            rhs = 0
            for i in self.inflows[(k, j)]:
                if math.isclose(first_vehicle_slot_set.index(i) % 2, 0):
                    rhs += self.vars[(k - 1, i, k, j)]
            self.model.addConstr(self.vars_l[(k, j)] == rhs, f"cons_5_{id}")
            id += 1
        # part 2
        for k in vehicle_id_sequence[2:]:
            for j in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                temp_sum_inflow = temp_sum_l = 0
                for i in self.inflows[(k, j)]:
                    temp_sum_inflow += self.vars[(k - 1, i, k, j)]
                    temp_sum_l += self.vars_l[(k - 1, i)]
                    self.model.addConstr(
                        self.vars_l[(k, j)] >= self.vars[(k - 1, i, k, j)] + self.vars_l[(k - 1, i)] - 1,
                        f"cons_5_{id}")
                    id += 1
                    self.model.addConstr(
                        self.vars[(k - 1, i, k, j)] >= self.vars_l[(k, j)] + self.vars_l[(k - 1, i)] - 1,
                        f"cons_5_{id}")
                    id += 1
                    self.model.addConstr(
                        self.vars_l[(k - 1, i)] >= self.vars_l[(k, j)] + self.vars[(k - 1, i, k, j)] - 1,
                        f"cons_5_{id}")
                    id += 1
                # new one
                self.model.addConstr(
                    self.vars_l[(k, j)] <= temp_sum_inflow, f"cons_5_{id}")
                id += 1
                self.model.addConstr(
                    self.vars_l[(k, j)] <= temp_sum_l, f"cons_5_{id}")
                id += 1
        # group 6 (replica links)
        id = 0
        for (k, i, j), (i1, j1) in self.img_links.items():
            if k == first_v_id:
                self.model.addConstr(self.vars[(k, i, k + 1, j)] == self.vars[(k, i1, k + 1, j1)], f"cons_6_{id}")
            else:
                self.model.addConstr(
                    self.vars[(k, i1, k + 1, j1)] >= self.vars[(k, i, k + 1, j)] + self.vars_l[(k, i)] - 1,
                    f"cons_6_{id}")
            id += 1
        # group 7 (cut set inequalities)
        id = 0
        for k in vehicle_id_sequence[1:]:
            sum_inflow = 0
            for i in range(self.slots_ranges[k][0], slots_ranges[k][1] + 1):
                for pre_i in self.inflows[(k, i)]:
                    sum_inflow += self.vars[(k - 1, pre_i, k, i)]
            self.model.addConstr(sum_inflow == outflow_count, f"cons_7_{id}")
            id += 1
        # group 8 (regulate margin of regions)
        id = 0
        region_num = outflow_count / 2
        for k in self.vehicle_id_sequence[-1:]:
            window_size = vehicle_point_sample_interval
            for j in range(slots_ranges[k][0], slots_ranges[k][1] + 1):
                sum_lhs_names = []
                for temp_j in range(slots_ranges[k][0], j + 1):
                    sum_lhs_names.append((k, temp_j))
                sum_rhs_names = []
                last_j = min(j + shift_interval + 1 + window_size + 1, slots_ranges[k][1] + 1)
                for temp_j in range(j + 1, last_j):
                    sum_rhs_names.append((k, temp_j))
                if len(sum_rhs_names) > 0:
                    rhs = quicksum([self.vars_l[var_name] for var_name in sum_rhs_names])
                    lhs = quicksum([self.vars_l[var_name] for var_name in sum_lhs_names])
                    self.model.addConstr(rhs >= self.vars_l[(k, j)] + 1 - lhs / region_num - 1, f"cons8_{id}")
                    id += 1

        # initial solution
        if benchmark is not None:
            self.set_initial_solution(benchmark)

        self.model.setObjective(objective, GRB.MINIMIZE)
        self.model.setParam(GRB.Param.TimeLimit, 300)
        self.model.update()
        self.model.write(output_folder + "model.lp")
        self.model.setParam(GRB.Param.LogFile, output_folder + "gurobi_log.txt")
        self.model.optimize()
        if (self.model.status == GRB.Status.OPTIMAL or self.model.status == GRB.Status.INTERRUPTED or
                self.model.status == GRB.Status.TIME_LIMIT):
            self.obj_val = self.model.objVal
            self.var_vals = {}
            for var_name, var in self.vars.items():
                self.var_vals[var_name] = round(var.X)
        elif self.model.status == GRB.Status.INFEASIBLE:
            self.model.computeIIS()
            self.model.write(output_folder + "model.ilp")
            raise Exception('model is infeasible!')


    def delete_trivial_links(self):
        need_continue = True
        while need_continue:
            # delete the links that only has inflow links
            for k in self.vehicle_id_sequence[::-1][1:-1]:
                k_point_size = len(self.vehicle_slots[k])
                for i in range(k_point_size):
                    if len(self.inflows[(k, i)]) > 0 and len(self.outflows[(k, i)]) == 0:
                        # delete all inflow links k,i,j
                        for j in self.inflows[(k, i)]:
                            # remove from outflows
                            if i in self.outflows[(k - 1, j)]:
                                self.outflows[(k - 1, j)].remove(i)
                            # update imaginary links
                            if (k, i, j) in self.img_links.keys():
                                self.recursive_delete_img_link(k, i, j)
                            if (k, i, j) in self.img_inverse_links.keys():
                                (ori_k, ori_i, ori_j) = self.img_inverse_links[(k, i, j)]
                                self.img_links.pop((ori_k, ori_i, ori_j))
                                self.img_inverse_links.pop((k, i, j))
                        self.inflows[(k, i)].clear()
            # delete the links that only has outflow links
            for k in self.vehicle_id_sequence[1:-1]:
                k_point_size = len(self.vehicle_slots[k])
                for i in range(k_point_size):
                    if len(self.outflows[(k, i)]) > 0 and len(self.inflows[(k, i)]) == 0:
                        # delete all outflow links (k,i,j)
                        for j in self.outflows[(k, i)]:
                            if i in self.inflows[(k + 1, j)]:
                                self.inflows[(k + 1, j)].remove(i)
                            # update imaginary links
                            if (k, i, j) in self.img_links.keys():
                                self.recursive_delete_img_link(k, i, j)
                            if (k, i, j) in self.img_inverse_links.keys():
                                (ori_i, ori_j) = self.img_inverse_links[(k, i, j)]
                                self.img_links.pop((k, ori_i, ori_j))
                                self.img_inverse_links.pop((k, i, j))
                        self.outflows[(k, i)].clear()
            need_continue = False
            # check the first vehicle
            first_v_id = self.vehicle_id_sequence[0]
            slot_list = []
            for i in range(self.slots_ranges[first_v_id][0], self.slots_ranges[first_v_id][1] + 1,
                           vehicle_point_sample_interval):
                if len(self.outflows[(first_v_id, i)]) > 0:
                    slot_list.append(i)
            # odd number
            for i in range(self.slots_ranges[first_v_id][0], slot_list[0]):
                if len(self.outflows[(first_v_id, i)]) > 0:
                    # remove links (first_v_id,i,j)
                    for j in self.outflows[(first_v_id, i)]:
                        if i in self.inflows[(first_v_id + 1, j)]:
                            self.inflows[(first_v_id + 1, j)].remove(i)
                        # update imaginary links
                        if (first_v_id, i, j) in self.img_links.keys():
                            self.recursive_delete_img_link(first_v_id, i, j)
                        if (first_v_id, i, j) in self.img_inverse_links.keys():
                            (ori_k, ori_i, ori_j) = self.img_inverse_links[(first_v_id, i, j)]
                            self.img_links.pop((ori_k, ori_i, ori_j))
                            self.img_inverse_links.pop((first_v_id, i, j))
                    self.outflows[(first_v_id, i)].clear()
                    need_continue = True

    def gen_img_links(self):
        first_v_id = self.vehicle_id_sequence[0]
        """generate imaginary (replica) links"""
        for k in self.vehicle_id_sequence[:-1]:
            i = self.slots_ranges[k][0]
            while i + shift_interval <= self.slots_ranges[k][1]:
                i1 = i + shift_interval
                # examine each outflow link of i
                for j in self.outflows[(k, i)]:
                    j1 = j + shift_interval
                    # add link (k,i1,j1)
                    if j1 <= self.slots_ranges[k + 1][1]:
                        self.img_links[(k, i, j)] = (i1, j1)
                        self.img_inverse_links[(k, i1, j1)] = (i, j)
                        if j1 not in self.outflows[(k, i1)]:
                            self.outflows[(k, i1)].append(j1)
                        if i1 not in self.inflows[(k + 1, j1)]:
                            self.inflows[(k + 1, j1)].append(i1)

                if k == first_v_id:
                    i += vehicle_point_sample_interval
                else:
                    i += 1

    def gen_qk_points(self):
        self.qk_points = []  # value: (q,k)
        # key: (k1,k2), value:[(x1,y1,x2,y2,x3,y3,x4,y4)], each row is for a small quadrangle
        self.regions = defaultdict(list)
        self.region_speeds = []  # id: region_id, value: [speeds]
        self.reaction_times = []
        self.wave_speeds = {}  # id: v_i - v_i+1
        # obtain regions
        for (k1, k2), lines in self.measurement_lines.items():
            self.wave_speeds[(k1, k2)] = []
            for m in range(0, len(lines) - 1, 2):
                (i_l, j_l) = lines[m]  # left i and j
                (i_r, j_r) = lines[m + 1]  # right i and j
                x_bl, y_bl = self.vehicle_slots[k2][j_l], self.vehicle_pos[k2][j_l]
                x_br, y_br = self.vehicle_slots[k2][j_r], self.vehicle_pos[k2][j_r]
                x_ul, y_ul = self.vehicle_slots[k1][i_l], self.vehicle_pos[k1][i_l]
                x_ur, y_ur = self.vehicle_slots[k1][i_r], self.vehicle_pos[k1][i_r]
                self.wave_speeds[(k1, k2)].append([(y_ul - y_bl) / (x_ul - x_bl), (y_ur - y_br) / (x_ur - x_br)])
                self.regions[(k1, k2)].append((x_ul, y_ul, x_bl, y_bl, x_br, y_br, x_ur, y_ur))  # anti-clockwise
        # get qk_points
        self.q_list = []
        self.k_list = []
        first_v_id = self.vehicle_id_sequence[0]
        for region_id in range(len(self.regions[(first_v_id, first_v_id + 1)])):
            k1 = first_v_id
            k2 = first_v_id + 1
            region_speeds = []
            q = k = 0
            total_area = 0
            while k2 <= self.vehicle_id_sequence[-1]:
                (x_ul, y_ul, x_bl, y_bl, x_br, y_br, x_ur, y_ur) = self.regions[(k1, k2)][region_id]
                p_ul_idx = int(round((x_ul - self.vehicle_slots[k1][0]) / step_size))
                p_ur_idx = int(round((x_ur - self.vehicle_slots[k1][0]) / step_size))
                p_bl_idx = int(round((x_bl - self.vehicle_slots[k2][0]) / step_size))
                p_br_idx = int(round((x_br - self.vehicle_slots[k2][0]) / step_size))

                self.reaction_times.append(0.5 * (x_bl - x_ul + x_br - x_ur))

                if k1 == first_v_id:
                    for i in range(p_ul_idx, p_ur_idx + 1):
                        region_speeds.append(self.vehicle_speeds[k1][i])

                for i in range(p_bl_idx, p_br_idx + 1):
                    region_speeds.append(self.vehicle_speeds[k2][i])

                temp_area = quadrilateral_area(x_ul, y_ul, x_bl, y_bl, x_br, y_br, x_ur, y_ur)
                total_area += temp_area

                # the last vehicle
                if k1 == self.vehicle_id_sequence[0]:
                    k += x_ur - x_ul
                    q += y_ur - y_ul

                k += x_br - x_bl
                q += y_br - y_bl

                # move to the next (k1,k2)
                k1 += 1
                k2 += 1
            plantoon_size = len(self.vehicle_slots)
            total_area *= ((plantoon_size) / (plantoon_size - 1))
            q = q / total_area
            k = k / total_area
            self.q_list.append(q * 3600)  # veh/s to veh/hr
            self.k_list.append(k * 1000)  # veh/m to veh/km
            self.qk_points.append((q, k))
            self.region_speeds.append(region_speeds)

        self.center_point = (np.mean(self.k_list), np.mean(self.q_list))
        groups = [i for i in range(len(self.region_speeds))]
        means = [np.mean(speeds) for speeds in self.region_speeds]
        std_dev = [np.std(speeds) for speeds in self.region_speeds]
        # return speed info
        return [groups, means, std_dev, self.wave_speeds]

    def update_slot_ranges(self):
        for k in self.vehicle_id_sequence:
            # for start index
            start_idx = i = 0
            while i < len(self.vehicle_slots[k]):
                # the first vehicle
                if k == self.vehicle_id_sequence[0] and len(self.outflows[(k, i)]) > 0:
                    start_idx = i
                    break
                elif k != self.vehicle_id_sequence[0] and len(self.inflows[(k, i)]) > 0:
                    start_idx = i
                    break
                i += 1
            # for end index
            end_idx = j = len(self.vehicle_slots[k]) - 1
            while j > 0:
                # the first vehicle
                if k == self.vehicle_id_sequence[0] and len(self.outflows[(k, j)]) > 0:
                    end_idx = j
                    break
                elif k != self.vehicle_id_sequence[0] and len(self.inflows[(k, j)]) > 0:
                    end_idx = j
                    break
                j -= 1
            self.slots_ranges[k] = [start_idx, end_idx]

    def recursive_delete_img_link(self, k, i, j):
        while (k, i, j) in self.img_links.keys():
            (next_i, next_j) = self.img_links[(k, i, j)]
            # delete link
            self.img_links.pop((k, i, j))
            self.img_inverse_links.pop((k, next_i, next_j))

    def set_initial_solution(self, benchmark):
        for key in self.vars.keys():
            self.vars[key].start = 0
        for (k, i, j) in benchmark.selected_links:
            key = (k, i, k + 1, j)
            if key in self.vars.keys():
                self.vars[key].start = 1
