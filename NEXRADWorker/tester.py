# Imports for proper program functionality
import nexradaws, config, pytz, pyart, tempfile
import matplotlib.pyplot as plt
from NEXRADWorker import NEXRADStationManager
from txtparsing import DataWorker
from datetime import datetime, timezone

templocation = tempfile.mkdtemp()

# Initializing variables and local folder environment
config.init_environment()
LOC_FOLS = config.LOC_FOLS
BING_KEY = config.BING_MAPS_API_KEY

worker = DataWorker(LOC_FOLS['meta']+'nexrad-stations-template.txt')

time = datetime(2016, 3, 23, 12, 34, tzinfo=timezone.utc)
manager = NEXRADStationManager([38.149284, -108.755224], [41.951239, -102.351951], time)
list = manager.bin_search_longitude(-108.755224, -102.351951)
NEXRADStationManager.cl_wd()
print(len(list))


# Example code from the nexradaws examples page
#conn = nexradaws.NexradAwsInterface()
#central_timezone = pytz.timezone('US/Central')
#radar_id = 'KTLX'
#start = central_timezone.localize(datetime(2013,5,31,17,0))
#end = central_timezone.localize (datetime(2013,5,31,19,0))
#scans = conn.get_avail_scans_in_range(start, end, radar_id)
#print("There are {} scans available between {} and {}\n".format(len(scans), start, end))
#print(scans[0:4])


#results = conn.download(scans[0:4], LOC_FOLS['nexrad'])
#print(results.success)

#fig = plt.figure(figsize=(16,12))
#for i,scan in enumerate(results.iter_success(),start=1):
#    ax = fig.add_subplot(2,2,i)
#    radar = scan.open_pyart()
#    display = pyart.graph.RadarDisplay(radar)
#    display.plot('reflectivity',0,ax=ax,title="{} {}".format(scan.radar_id,scan.scan_time))
#    display.set_limits((-250, 250), (-250, 250), ax=ax)
#plt.show()

# Takes around 30 seconds to generate the 4 plot display in this example