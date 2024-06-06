import json
import random
import math
import numpy as np
import scipy.stats as stats
import utm
from geopy.distance import geodesic
from scipy.optimize import fsolve
from shapely.geometry import Polygon, MultiPolygon, mapping
from shapely.ops import unary_union
import joblib


class MCBackend:
    """
    This class contains the backend functions for the Monte Carlo simulation. 
    Ready to be used from the FrontEnd interface or the BackendOffline script.
    It generates lightning strikes, calculates the TDOA for the given network of stations, and saves the results.
    """

    def __init__(self, input_stations_file, curr_dist_parameters_file, svm_parameters_negLS, svm_parameters_posLS):
        self.mu, self.sigma = 0, 1e-6               # mean and standard deviation for the error of detection in the TOA of a station
        with open(curr_dist_parameters_file, 'r') as json_file:
            parameters = json.load(json_file)
            self.pos_params = parameters[0]         # positive current distribution parameters are stored at 1st element
            self.neg_params = parameters[1]         # negative current distribution parameters are stored at 2nd element
        self.min_curr = 2E3                         # min current allowed for lightning strike
        self.max_curr = 200E3                       # max current allowed for lightning strike
        self.svm={}
        self.read_pkl(svm_parameters_posLS, 'pos')  # read the SVM model for positive lightning strikes
        self.read_pkl(svm_parameters_negLS, 'neg')  # read the SVM model for negative lightning strikes
        with open(input_stations_file, 'r') as json_file:
            self.stations = json.load(json_file)
        self.input_stations = {}
        self.mode = True                            # True: implies all detectors record all lightning strikes / False: uses the svm as detector pattern method (True by default)

    def read_pkl(self, filename, polarity):
        """
        This function reads the SVM model from a .pkl file and stores the parameters in the backend object.
        The important parameters of the model are the ones that define the separating line: w1, w2, and b.
        In case of another model is generated, ensure to have separate behavior for positive and negative lightning strikes.
        """
        svm_model = joblib.load(filename)
        w1 = svm_model.coef_[0][0]
        w2 = svm_model.coef_[0][1]
        b = svm_model.intercept_[0]
        self.svm[polarity] = {'w1': w1, 'w2': w2, 'b': b}
    
    def get_regions(self):
        """Returns the list of regions in the station network"""
        regions_set = set()
        for station in self.stations:
            regions_set.add(station['reg'])
        return list(regions_set)

    def filter_stations_by_region(self, region):
        """Returns the list of stations in the given region"""
        return [station for station in self.stations if station['reg'] == region]

    def add_input_stations(self,input_stations):
        """Add the input stations to the backend object for the solution"""
        self.input_stations = input_stations

    def calculate_average_coordinates(self, selected_stations):
        """Calculate the average coordinates of the selected stations"""
        if not selected_stations:
            return None
        avg_lat = sum(station['lat'] for station in selected_stations) / len(selected_stations)
        avg_lon = sum(station['lon'] for station in selected_stations) / len(selected_stations)
        return [avg_lat, avg_lon]

    def generate_random_curr(self):
        """Generate a random current value from the lognormal distribution of the lightning strikes"""
        random_number = random.random()
        if random_number < 0.9:         # Important: 90% of generated currents are negative
            shape, loc, scale = self.neg_params['shape'], self.neg_params['loc'], self.neg_params['scale']
            while True:
                random_value = stats.lognorm.rvs(shape, loc=loc, scale=scale)
                if self.min_curr <= random_value <= self.max_curr:
                    break
            return -random_value
        else:
            shape, loc, scale = self.pos_params['shape'], self.pos_params['loc'], self.pos_params['scale']
            while True:
                random_value = stats.lognorm.rvs(shape, loc=loc, scale=scale)
                if self.min_curr <= random_value <= self.max_curr:
                    break
            return random_value

    def create_boundary(self, offset=0.2):      # offset is the degrees of latitude and longitude to create the square
        """Create a boundary around the input stations to generate random lightning strikes within the boundary"""
        all_squares = []
        for coord in self.input_stations:        # Create squares around each coordinate
            lat = coord['lat']
            lon = coord['lon']
            square = Polygon([
                (lon - offset, lat - offset),
                (lon + offset, lat - offset),
                (lon + offset, lat + offset),
                (lon - offset, lat + offset)
            ])
            all_squares.append(square)
        boundary = unary_union(all_squares)             # Remove overlapping parts and aggregate into a single polygon
        if isinstance(boundary, MultiPolygon):          # Convert MultiPolygon to Polygon if necessary
            boundary = boundary.convex_hull
        return mapping(boundary)['coordinates'][0]

    def generate_random_coordinate_within_boundaries(self):
        """Generate a random coordinate within the boundary created around the input stations"""
        boundary_flat = self.create_boundary()
        min_lat = min([point[1] for point in boundary_flat])
        max_lat = max([point[1] for point in boundary_flat])
        min_lon = min([point[0] for point in boundary_flat])
        max_lon = max([point[0] for point in boundary_flat])
        while True:
            lat = random.uniform(min_lat, max_lat)
            lon = random.uniform(min_lon, max_lon)
            point_to_check = (lon, lat)
            if self.is_point_within_boundary(point_to_check, boundary_flat):
                return (lat, lon)

    def is_point_within_boundary(self, point, boundary):
        """Check if the point is within the boundary polygon"""
        x, y = point
        n = len(boundary)
        inside = False
        p1x, p1y = boundary[0]
        for i in range(n + 1):
            p2x, p2y = boundary[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    def convert_xy(self, point):
        x = utm.from_latlon(point[0], point[1])[0]
        y = utm.from_latlon(point[0], point[1])[1]
        return {'x': round(x, 3), 'y': round(y, 3)}

    def convert_sta_xy(self, stations):
        stations_xy = []
        for station in stations:
            point = (station['lat'], station['lon'])
            station_xy = self.convert_xy(point)
            station_xy['id'] = station.get('id', None)
            stations_xy.append(station_xy)
        return stations_xy

    def random_ls_with_coords(self, id):
        """Generate a random lightning strike with random coordinates and current value within the boundaries of the input stations"""
        ls_coord = self.generate_random_coordinate_within_boundaries()
        ls_coord_round = (round(ls_coord[0], 6), round(ls_coord[1], 6))
        xy = self.convert_xy(ls_coord_round)
        curr = round(self.generate_random_curr(), 0)
        random_ls = {"id": id, "timestamp": 0, "curr": curr, "x": xy['x'], "y": xy['y'],
                     "lat": ls_coord_round[0], "lon": ls_coord_round[1]}
        return random_ls

    def distcoord(self, coord1, coord2):
        return geodesic(coord1, coord2).kilometers  # calculate the distance between two coordinates in kilometers

    def distcart(self, p1, p2):
        """Calculate the distance between two cartesian points in meters"""
        x1, y1 = p1
        x2, y2 = p2
        dist_m = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        return dist_m

    def curr_x_dist_lookup(self, curr, dist_m):
        """Check if the current and distance are within the SVM model boundaries for the lightning strike detection"""
        dist_km = dist_m / 1e3
        if curr >= 0:
            predicted_dist = self.svm['pos']['w1'] * curr + self.svm['pos']['w2'] * dist_km + self.svm['pos']['b']
        else:
            predicted_dist = self.svm['neg']['w1'] * curr + self.svm['neg']['w2'] * dist_km + self.svm['neg']['b']
        return predicted_dist > 0 or self.mode

    def check_sta(self, new_ls, station_xy, id):
        """Check if the station detects the lightning strike and calculate the time of arrival (TOA) if detected"""
        p1 = (station_xy['x'], station_xy['y'])
        p2 = (new_ls['x'], new_ls['y'])
        dist_m = self.distcart(p1, p2)
        curr = new_ls['curr']
        s = round(np.random.normal(self.mu, self.sigma), 9)
        if self.curr_x_dist_lookup(curr, dist_m):
            toa = round(dist_m / 3e8 + s, 9)
            return {"sta": station_xy.get('id', id), "time": toa, "x": station_xy['x'], "y": station_xy['y'], 'err': s}
        else:
            return None

    def init_ls(self, stations_xy, num_ls):
        """Initialize the lightning strikes with random coordinates and current values within the boundaries of the input stations"""
        lstrikes = []
        for idx in range(0, num_ls):
            new_ls = self.random_ls_with_coords(idx)
            stations_that_arrived = []
            for id, sta in enumerate(stations_xy):
                station = self.check_sta(new_ls, sta, id)
                if station is not None:
                    stations_that_arrived.append(station)
            new_ls['nsta'] = len(stations_that_arrived)
            new_ls['sig'] = stations_that_arrived
            lstrikes.append(new_ls)
        return lstrikes

    def func(self, variables, *args):
        """Function to solve the TDOA equations using fsolve (geometrical solution with 3 stations)"""
        (x, y) = variables
        (x0, x1, x2, y0, y1, y2, d10, d20) = args
        eqn_1 = (((x - x2) ** 2 + (y - y2) ** 2) ** 0.5 - ((x - x0) ** 2 + (y - y0) ** 2) ** 0.5 - d20)
        eqn_2 = (((x - x1) ** 2 + (y - y1) ** 2) ** 0.5 - ((x - x0) ** 2 + (y - y0) ** 2) ** 0.5 - d10)
        return [eqn_1, eqn_2]
    
    def tdoa_dont_solve(self, ls1):
        """If the number of stations that detected the lightning strike is less than 3, the solution is not performed"""
        statlist = ls1['sig']
        ls1['x_calc'] = None
        ls1['y_calc'] = None
        ls1['tdoa_estim_error'] = None
        if len(statlist) > 0:
            J = np.array([d['err'] for d in statlist])
            ls1['avg_sta_error'] = round(np.average(J), 9)
        else:
            ls1['avg_sta_error'] = None
        return ls1

    def centroid(self,x,y):
        """Calculate the centroid of the triangle formed by the 3 stations for initial guess of geometrical solution"""
        if len(x) != 3 or len(y) != 3:
            raise ValueError("Please provide exactly three vertices.")
        x_centroid = sum(x) / 3
        y_centroid = sum(y)  / 3
        return x_centroid, y_centroid

    def tdoa_solve_geo(self, ls1):
        """Solve the TDOA equations using geometrical solution with 3 stations"""
        statlist = ls1['sig']
        xs1 = np.array([d['x'] for d in statlist])
        ys1 = np.array([e['y'] for e in statlist])
        ts1 = [f['time'] for f in statlist]
        xs = [x for _, x in sorted(zip(ts1, xs1))]
        ys = [x for _, x in sorted(zip(ts1, ys1))]
        ts = sorted(ts1)
        x0 = xs[0]
        y0 = ys[0]
        t0 = ts[0]
        tds = np.delete(ts, 0) - t0
        rtoa = np.array(tds) * 3e8
        initial_sol = np.array([self.centroid(xs,ys)]) 
        res = fsolve(self.func, initial_sol, args=(xs[0], xs[1], xs[2], ys[0], ys[1], ys[2], rtoa[0], rtoa[1]))
        ls1['x_calc'] = round(res[0], 3)
        ls1['y_calc'] = round(res[1], 3)
        tdoa_estim_error = ((ls1['x'] - ls1['x_calc']) ** 2 + (ls1['y'] - ls1['y_calc']) ** 2) ** 0.5
        ls1['tdoa_estim_error'] = round(tdoa_estim_error, 3)
        J = np.array([d['err'] for d in statlist])
        ls1['avg_sta_error'] = round(np.average(J), 9)
        return ls1
    
    def tdoa_solve_mat(self, ls1):
        """Solve the TDOA equations using matrix solution with more than 3 stations
        
        The method tdoa_solve can't be disclosed, but it should act in the ls dictionary to solve the TDOA equations when more than 3 stations are present
        
        ls1 has the format: (example)
        {   "id": 0,"timestamp": 0,"curr": 17827.0,"x": 689477.073,"y": 6085648.072,"lat": -35.354494,"lon": 149.085304,"nsta": 5,
            "sig": [
                {"sta": 44,"time": 7.4993e-05,"x": 692400.408,"y": 6107982.76,"err": -9.1e-08},
                {"sta": 45, "time": 4.4022e-05, "x": 677634.547, "y": 6089619.274,"err": 2.387e-06},
                {"sta": 46, "time": 2.7301e-05, "x": 697381.107, "y": 6082732.569, "err": -7.81e-07},
                {"sta": 47, "time": 8.4065e-05, "x": 666050.89, "y": 6076376.349, "err": 8.4e-08},
                {"sta": 49, "time": 6.0279e-05, "x": 687418.052, "y": 6068488.102, "err": 2.669e-06}
                ]
        }
                
        and the method should provide:
        ls1['x_calc'] X cartesian coordinate of the solution 
        ls1['y_calc'] Y cartesian coordinate of the solution
        ls1['tdoa_estim_error'] error of the solution, can be calculated as the euclidean distance between the calculated and the real position (as in the geometrical solution)
        ls1['avg_sta_error'] average error of the stations, can be calculated as the average of the error of each station        
        
        """
        from tdoa_solve import tdoa_solve    # a method for solving stations with more than 3 detectors is needed to be imported, see method tdoa_solve_mat for more detail
        
        return tdoa_solve(ls1)

    def commander(self, num_ls):
        """Commander function to generate lightning strikes and calculate the TDOA for the given network of stations"""
        stations_xy = self.convert_sta_xy(self.input_stations)
        ls = self.init_ls(stations_xy, num_ls)
        lightning_dic = ls
        return lightning_dic

    def save(self, lstrikes, filename):
        """Save the lightning data to a JSON file"""
        print('------------------------------------------------------------------------------------')
        json_string = json.dumps(lstrikes)
        with open(filename, "w") as json_file:
            json_file.write(json_string)
            json_file.close()
        print("Lightning Data was saved at ", filename)
        print('------------------------------------------------------------------------------------')

    def solver(self, strike_list):
        """Solve the TDOA equations for the lightning strikes in the list of strikes"""
        lstrikes = []
        for i in range(0, len(strike_list)):
            single_ls = strike_list[i]
            if len(single_ls['sig']) > 3:
                try:
                    ls1 = self.tdoa_solve_mat(single_ls)
                    tdoa_warning = False
                except:
                    tdoa_warning = True
                    single_ls['sig'] = single_ls['sig'][:3]
                    ls1 = self.tdoa_solve_geo(single_ls)
            elif len(single_ls['sig']) == 3:
                ls1 = self.tdoa_solve_geo(single_ls)
            else:
                ls1 = self.tdoa_dont_solve(single_ls)
            lstrikes.append(ls1)
        return lstrikes, tdoa_warning 