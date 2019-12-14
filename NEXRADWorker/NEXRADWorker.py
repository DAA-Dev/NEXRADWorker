# Imports for proper program functionality
import nexradaws, config, math, logging, os
from txtparsing import DataWorker

LOC_FOLS = config.LOC_FOLS
TAG = 'NEXRADWorker - '

# Class to manage all the files downloaded, deleting and updating as necessary
class NEXRADStationManager():
    # Get the metadata for radar stations and store as a list of RadarStations, sorted from lowest
    # longitudinal value to greatest
    def __init__(self, left_bot_cor, right_top_cor, time):
        self.__radar_stations = [] # All the radar stations, generated from metadata
        self.__sim_time = time     # Current time for the simulation

        worker = DataWorker(LOC_FOLS['meta']+'nexrad-stations-template.txt')
        worker.quicksort_lg(LOC_FOLS['meta']+'nexrad-stations.txt',
                            LOC_FOLS['meta']+'nexrad-stations-sorted.txt',
                            'longitude')
        worker.replace(LOC_FOLS['meta']+'nexrad-stations.txt',
                       LOC_FOLS['meta']+'nexrad-stations-sorted.txt')
        meta_data = worker.get_vals(LOC_FOLS['meta']+'nexrad-stations.txt', 
                                    ['icao', 'state', 'elevation', 'latitude', 'longitude'])
        for station_list in meta_data:
            self.__radar_stations.append(RadarStation(station_list[0], station_list[1],
                                                      station_list[2], station_list[3],
                                                      station_list[4]))
        self.update_area(left_bot_cor, right_top_cor)

    # Updates the area for which relevant radar stations are found
    def update_area(self, left_bot_cor, right_top_cor, pull_new_data=False):
        self.__relevant_stations = []
        right_top_bounds = self.get_point(right_top_cor, 45, 231.5)
        left_bottom_bounds = self.get_point(left_bot_cor, 225, 231.5)
        logging.info(TAG+'updating the relevant radar stations based on area parameters')
        logging.info(TAG+'getting all the stations in the bounds: ' + str(left_bottom_bounds) + ', ' + str(right_top_bounds))

        lon_stations = self.bin_search_longitude(left_bottom_bounds[1], right_top_bounds[1])
        for station in lon_stations:
           if station.latitude > left_bottom_bounds[0] and station.latitude < right_top_bounds[0]:
               self.__relevant_stations.append(station)

        if pull_new_data:
            self.pull_new_data()

    # Updates the time for the NEXRAD data with a new time passed as an argument
    def update_time(self, time, pull_new_data=False):
        self.__sim_time = time
        if pull_new_data:
            self.pull_new_data()
    
    # Steps the time of the object by the passed timedelta argument
    def step_time(self, timedelta, pull_new_data=False):
        self.__sim_time = self.__sim_time + timedelta
        if pull_new_data:
            self.pull_new_data()

    # Find the stations within a range of longitudes 
    def bin_search_longitude(self, lower_bound, upper_bound):
        logging.info(TAG+'searching NEXRAD stations based on latitude')
        # Finds the closest value to the target from RadarStations in array
        # Modified version of the classic recursive binary search algorithm
        def bin_search_stations(array, left, right, target):
            middle = int((left + right) / 2)
            if middle is left:
                if abs(array[right].longitude - target) < abs(array[left].longitude - target):
                    return right
                else:
                    return left
            if right >= 0:
                if array[middle].longitude > target:
                    return bin_search_stations(array, left, middle - 1, target)
                if array[middle].longitude < target:
                    return bin_search_stations(array, middle + 1, right, target)
            else:
                return 0
        
        # Runs recursive searches to find closest points to lower and upper bounds, and then uses the
        # indecies returned from the binary search to return a list of RadarStation objects
        left = bin_search_stations(self.__radar_stations, 0, len(self.__radar_stations) - 1, lower_bound)
        right = bin_search_stations(self.__radar_stations, 0, len(self.__radar_stations) - 1, upper_bound)
        results = []
        if left is right:
            results.append(self.__radar_stations[left])
            return results
        else:
            for x in range(left, right):
                results.append(self.__radar_stations[x])
            return results

    # Function to get a new gps point based on an origin gps point, a bearing from that point
    # as well as the distance travelled along that bearing, used to calculate max bounds
    # Bearing - degrees clockwise from true north
    def get_point(self, origin_point, bearing, distance):
        earth_radius = 6378.1
        bearing = math.radians(bearing)

        lat1 = math.radians(origin_point[0]) 
        lon1 = math.radians(origin_point[1]) 

        lat2 = math.asin(math.sin(lat1)*math.cos(distance/earth_radius) +
                math.cos(lat1)*math.sin(distance/earth_radius)*math.cos(bearing))
        lon2 = lon1 + math.atan2(math.sin(bearing)*math.sin(distance/earth_radius)*math.cos(lat1),
                        math.cos(distance/earth_radius)-math.sin(lat1)*math.sin(lat2))

        lat2 = math.degrees(lat2)
        lon2 = math.degrees(lon2)
        return [lat2, lon2]


    # Method to pull data for each station from AWS
    # First, queries AWS and logs the relevant queries, then initiates downloads of relevant files
    def pull_new_data(self):
        for station in self.__relevant_stations:
            print()
        return None

    # Clears the downloaded files from local directory for NEXRAD data
    @staticmethod
    def cl_wd():
        logging.info(TAG+'clearing downloaded files in the NEXRAD directory')
        for file in os.listdir(LOC_FOLS['nexrad']):
            try:
                os.unlink(LOC_FOLS['nexrad']+file)
            except:
                logging.error(TAG+'error in deleting nexrad files, check for memory leak')

class RadarStation():
    def __init__(self, icao, state, elevation, latitude, longitude):
        self.icao = icao
        self.state = state
        self.elevation = elevation
        self.latitude = latitude
        self.longitude = longitude

    @property
    def elevation(self):
        return self._elevation

    @property
    def latitude(self):
        return self._latitude

    @property
    def longitude(self):
        return self._longitude

    @elevation.setter
    def elevation(self, elevation):
        self._elevation = int(elevation)

    @latitude.setter
    def latitude(self, latitude):
        self._latitude = float(latitude)

    @longitude.setter
    def longitude(self, longitude):
        self._longitude = float(longitude)

    def create_map(self):
        return None

    def __str__(self):
        printout  = '*****' + self.icao  + '*****' + ' \n'
        printout += 'Longitude: ' + str(self.longitude) + '\n'
        printout += 'Latitude: ' + str(self.latitude) + '\n'
        return printout

# Important metadata:
#   - NEXRAD stations have a range of 231.5 km
# GPS coordinates are in the following format (lat, lon)
# When looking at a Mercator projection that correlates to (y, x)