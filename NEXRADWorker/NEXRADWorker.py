# Imports for proper program functionality
import nexradaws, config, math, logging, os, pyart
import matplotlib.pyplot as plt
from txtparsing import DataWorker

LOC_FOLS = config.LOC_FOLS
TAG = 'NEXRADWorker - '

# Simple sign-extension method
def s_ext(str, length):
    while len(str) != length:
        str = '0'+str
    return str

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
        year = self.__sim_time.year
        month = self.__sim_time.month
        day = self.__sim_time.day
        self.cl_wd()
        new_relevant_stations = []

        logging.info(TAG+'starting a pull of new data')
        aws_interface = nexradaws.NexradAwsInterface()
        radar_list = aws_interface.get_avail_radars(year, month, day)
        # Updates the scans for each station to ones closest to the target time
        for station in self.__relevant_stations:
            name = station.icao
            if name in radar_list:
                scans = aws_interface.get_avail_scans(s_ext(str(year), 4), s_ext(str(month), 2), s_ext(str(day), 2), name)
                closest_time_delta = abs((self.__sim_time - scans[0].scan_time).total_seconds()) 
                closest_scan  = scans[0]
                for index in range(1, len(scans)-1):
                    delta = abs((self.__sim_time - scans[index].scan_time).total_seconds()) 
                    if delta > closest_time_delta:
                        break
                    else:
                        closest_time_delta = delta
                        closest_scan = scans[index]
                logging.info(TAG+name)
                logging.info(TAG+'target time: ' + str(self.__sim_time))
                logging.info(TAG+'scan time: ' + str(closest_scan.scan_time))
                station.scan = closest_scan
                new_relevant_stations.append(station)

        # Have each RadarStation update itself, downloading its respective scan
        self.__relevant_stations = new_relevant_stations
        for station in self.__relevant_stations:
            station.pull_data()
        logging.info(TAG+'pulled data for all the stations in the station manager')
        del aws_interface

    # Method that combines all the scans in the 'nexrad' folder and creates an overlay for our simulator
    def combine_scans_into_overlay(self):
        local_scans = []
        for station in self.__relevant_stations:
            scan = station.create_map()
            if scan is not None:
                local_scans.append(scan)
        local_scans = tuple(local_scans)
        # None of this works. Figure out how to make it work. 


        #grid = pyart.map.grid_from_radars(local_scans,
        #                                  grid_shape=(1, 201, 201),
        #                                  grid_limits=((100000, 100000), (-100000, 100000), (-100000, 100000)))
        #print(type(grid))
        #print(str(grid.fields))

        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.imshow(grid.fields['reflectivity']['data'][0], origin='lower', extent=(-60, 40, -50, 40), vmin=0, vmax=48)
        #plt.show()


    # Method to print metadata for all the relevant stations
    def print_scan_data(self):
        for station in self.__relevant_stations:
            station.print_scan_data()

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
        self.scan = None

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

    def pull_data(self):
        aws_interface = nexradaws.NexradAwsInterface()
        try:
            self.local_data = aws_interface.download([self.scan], LOC_FOLS['nexrad'])._successfiles[0]
        except:
            self.local_data = None
        del aws_interface

    # Work here to create subplots, and then later to return proper object and metadata in order to data
    # from multiple stations
    def create_map(self):
        if self.local_data is not None:
            radar = self.local_data.open_pyart()
            return radar
            # Code to display a plot from one individual station
            #display = pyart.graph.RadarDisplay(radar)
            #fig = plt.figure(figsize=(6,5))
            #ax = fig.add_subplot(111)
            #display.plot('reflectivity', 0, title='', vmin=-32, vmax=64, colorbar_label='', ax=ax)
            #plt.show()
        else:
            return None

    def print_scan_data(self):
         if self.local_data is not None:
             logging.info(TAG+'Radar ' + self.icao + ' has data stored locally at: ' + self.local_data.filepath)
             logging.info(TAG+'Scan time is: ' + str(self.local_data.scan_time))

    def __str__(self):
        printout  = '*****' + self.icao  + '*****' + ' \n'
        printout += 'Longitude: ' + str(self.longitude) + '\n'
        printout += 'Latitude: ' + str(self.latitude) + '\n'
        return printout

# Important metadata:
#   - NEXRAD stations have a range of 231.5 km
# GPS coordinates are in the following format (lat, lon)
# When looking at a Mercator projection that correlates to (y, x)

# To improve the program runtime:
#           - Look at possibly pulling metadata from AWS, not a local required txt file