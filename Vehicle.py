class Vehicle:
    def __init__(self, v_id, lane_id):
        self.Vehicle_ID = v_id
        self.Frame_IDs = []
        self.Lane_ID = lane_id
        self.Local_Y = []
        self.Mean_Speed = []
        self.Follower_IDs = []
        self.Leader_IDs = []
        self.Follower_ID = -1
        self.Leader_ID = -1
        self.Accelerations = []
        # for analytical method
        self.Local_Y_nominal = []
        self.Local_Y_oscillation = []
        # for OpenAcc
        self.Latitudes = []  # degree
        self.Longitudes = []
