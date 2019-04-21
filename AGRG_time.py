import datetime

def time_datetime2GPS(dt):
	#reverse of time_gps2datetime
	currentLeapSeconds = 18
	
	gps0 = datetime.datetime(1980, 1, 6, 0, 0, 0)
	gpsWeek = int((dt - gps0).days/7)
	gps1 = datetime.timedelta(weeks=gpsWeek)
	gps2 = gps0+gps1
	gpsSecond = int((dt - gps2).total_seconds()-currentLeapSeconds)
	
	return [gpsWeek, gpsSecond]

def time_gps2datetime(gpsWeek,gpsSeconds):
	#WORKS!
	#gps time starts January 6, 1980 @ 00:00:00
	#differs from UTC by leap seconds https://www.nist.gov/pml/time-and-frequency-division/atomic-standards/leap-second-and-ut1-utc-information
	#leap seconds must be maintained as the rotational speed of the earth varies 
	#currently 18 seconds (as of july 2017)
	
	currentLeapSeconds = 18
	
	gps0 = datetime.datetime(1980, 1, 6, 0, 0, 0)
	gps1 = datetime.timedelta(weeks=gpsWeek, seconds=gpsSeconds+currentLeapSeconds)
	gps2 = gps0+gps1
	# print gps2
	return gps2

def time_tstring2datetime(tstring):
	#bit of a WIP
	#format:
	#YYYYMMDDHHMMSSIIIIII(where I is microseconds [10e-6 seconds])
	
	#input is padded to full 20 chars (using 0's)
	if len(tstring)<4:
		tstring = tstring + '0'*(4-len(tstring))
	# print tstring
	if len(tstring)<8:
		tstring = tstring + '01'*(6-len(tstring))
	# print tstring
	tstringPad = tstring + '0'*(20-len(tstring))
	
	# print 'YYYYMMDDHHMMSSIIIIII'
	# print tstring
	# print tstringPad
	
	dt = datetime.datetime.strptime(tstringPad, '%Y%m%d%H%M%S%f')
	# print dt
	return dt