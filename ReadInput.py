import csv
import math
from collections import defaultdict

import pandas as pd

import CommonDefines
from CommonDefines import *
from Vehicle import *

slot_id_lb = 90
slot_id_ub = slot_id_lb + read_point_num


def readNgsim(file_path):
    vehicles = {}  # key: vehicle id; value: vehicle
    vehicles4lane = defaultdict(list)  # key: lane id; value: vehicle id
    max_speed = -math.inf

    with open(file_path, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        # each row contains a vehicle
        vehicle_id = 0
        lane_id = 1
        skip = True

        for row in csvreader:
            vehicle = Vehicle(vehicle_id, lane_id)
            slots = [0.0]
            positions = [float(row[slot_id_lb]) * scale]  # ft to m
            speeds = [0.0]  # m/s
            accelerations = [0.0]  # m/s^2
            for i in range(slot_id_lb, len(row) - 1):
                if i > slot_id_ub:
                    break
                slots.append(slots[-1] + step_size)
                positions.append(float(row[i + 1]) * scale)
                speeds.append((positions[-1] - positions[-2]) / step_size)
                max_speed = max(max_speed, speeds[-1])

            speeds[0] = speeds[1]
            for i in range(len(speeds) - 1):
                accelerations.append((speeds[i + 1] - speeds[i]) / step_size)
            vehicle.Frame_IDs = slots
            vehicle.Local_Y = positions
            vehicle.Mean_Speed = speeds
            vehicle.Accelerations = accelerations

            vehicles[vehicle_id] = vehicle
            vehicles4lane[lane_id].append(vehicle_id)

            vehicle_id += 1
    CommonDefines.max_speed = max_speed
    return vehicles, vehicles4lane


def readACC(file_path):
    '''read OpenAcc dataset'''
    vehicles = {}  # key: vehicle id; value: vehicle
    vehicles4lane = defaultdict(list)  # key: lane id; value: vehicle id
    max_speed = -math.inf
    min_speed = math.inf

    df = pd.read_csv(file_path)
    vehicle_num = 6
    lane_id = 1
    for vehicle_id in range(vehicle_num):
        vehicle = Vehicle(vehicle_id, lane_id)
        vehicles[vehicle_id] = vehicle
        vehicles4lane[lane_id].append(vehicle_id)
    # for each row
    origin_lat = origin_lon = -1
    base_t = -1
    for index, row in df.iterrows():
        if index < slot_id_lb or index > slot_id_ub:
            continue
        slot = float(row['Time'])
        if index == slot_id_lb:
            base_t = slot
        slot -= base_t
        # determine the origin point in the coordinates
        if index == slot_id_lb:
            origin_lat = math.radians(float(row[f'Lat{vehicle_num}']))
            origin_lon = math.radians(float(row[f'Lon{vehicle_num}']))

        for vehicle_id in range(vehicle_num):
            vehicle = vehicles[vehicle_id]
            vehicle.Frame_IDs.append(slot)
            vehicle.Mean_Speed.append(float(row[f'Speed{vehicle_id + 1}']))
            lat, lon = math.radians(float(row[f'Lat{vehicle_id + 1}'])), math.radians(float(
                row[f'Lon{vehicle_id + 1}']))
            vehicle.Latitudes.append(lat)
            vehicle.Longitudes.append(lon)
            # Haversine Formula to get the distance
            dlat = lat - origin_lat
            dlon = lon - origin_lon
            a = math.sin(dlat / 2) ** 2 + math.cos(origin_lat) * math.cos(lat) * math.sin(dlon / 2) ** 2
            c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
            # Radius of earth in meters. Use 6371 for kilometers
            R = 6371000
            distance = R * c
            vehicle.Local_Y.append(distance)

            max_speed = max(max_speed, vehicle.Mean_Speed[-1])
            min_speed = min(min_speed, vehicle.Mean_Speed[-1])
            if index == slot_id_lb:
                vehicle.Accelerations.append(0.0)
            else:
                vehicle.Accelerations.append((vehicle.Mean_Speed[-1] - vehicle.Mean_Speed[-2]) / step_size)
    CommonDefines.max_speed = max_speed
    CommonDefines.min_speed = min_speed

    return vehicles, vehicles4lane
