########################################################################################################################################################
# Last updated: Aug 1, 2022 by Xinyu Liu
# This file specifies the input data to all other files required for trajectory-based simulation, with:
# Configuration -- Double-Sided 0-degree with short spots.
########################################################################################################################################################

import math
from math import ceil
import random
import numpy as np

dirname = ''

suffix = ''

debug = False
VEHICLE_IDX = None

########################################################################################################################################################
''' The following section specifies a default setup but user may make changes '''

''' number of vehicles types '''
N_vt = 2

''' probability of replacement vehicle type '''
frac_v = [0.2, 1]		

''' preferred alternative for spot placement '''
# place = 'entr'
place = 'exit'
# place = 'user'

''' admission policy '''
# pol = 'P' # primary only
# pol = 'F' # feasible only
# pol = 1 # set an integer for k

''' initialize w_0 '''
w_0 = 0
# w_0 = 1

''' single vs double sided '''
# side = 'single'
side = 'double'

''' 0-degree vs 90-degree '''
angle = 0
# angle = 90

''' only valid for 0-degree '''
''' long vs short spots '''
# mode = 'long'
mode = 'short'

if angle == 0:

	''' number of spot types '''
	N_st = 2 

	''' infeasible (0) / primary (1) / secondary (2) spot types for vehicle types '''
	vs_mat = [[1, 0],
			  [2, 1]] 

else:
	assert angle == 90

	''' number of spot types '''
	N_st = 3 

	''' infeasible (0) / primary (1) / secondary (2) spot types for vehicle types '''
	vs_mat = [[1, 0, 0],
			  [2, 1, 0]] 

''' mean service time '''
meanSERV = 60. # sec
# meanSERV = 90. # sec
# meanSERV = 300. # sec

''' average (desired) speed '''
avgSPEED = 44./3. # ft/sec

''' types of input distributions '''
# simType = 'det'
# simType = 'cav'
simType = 'exp'
# simType = 'unif'

''' simulation length and warmup period length (if any) '''
WARMUP_HOUR = 2 # hr
SIM_HOUR = 20 # hr
SIM_ITER = 20

''' max number of spots on each side to be evaluated '''
if side == 'single':
	LOT_TOTAL = 40
else:
	assert side == 'double'
	LOT_TOTAL = 20

if mode == 'singlelane':
	LOT_TOTAL = 60
########################################################################################################################################################
# The following section of the file should not be changed unless you know exactly what you are trying to do.
# WARNING: Unintended changes to this section can have unexpected effects to the simulation results.
if angle == 0:
	filename = '%s_%s_%s_%s_%s_%s'%(angle, mode, side, meanSERV, simType, suffix)
elif angle == 90:
	filename ='%s_%s_%s_%s_%s'%(angle, side, meanSERV, simType, suffix)

''' full vs partial access control '''
# control = 'partial' 
control = 'full'

ptype = 0
# ptype = 1

if control == 'partial' and (not (angle == 0 and mode == 'long')):
	pcontrol = control + str(ptype)
else:
	pcontrol = control

''' delay start? '''
delay = False

''' free curb? '''
# free_curb = True
free_curb = False
if free_curb:
	suffix = 'free_curb'

''' sample path matching with MC simulation? '''
spmatch = False
if spmatch:
	suffix += 'spmatch'
	assert simType == 'cav'
	WARMUP_HOUR = 0 # hr
	SIM_HOUR = 100 # hr
	SIM_ITER = 1

WARMUP_UNIT = 3600 * WARMUP_HOUR # sec
SIM_UNIT = 3600 * SIM_HOUR # sec

########################################################################################################################################################
# The following section of the file should not be changed unless you know exactly what you are trying to do.
# WARNING: Unintended changes to this section can have unexpected effects to the simulation results.

''' configuration specific parameters '''
if angle == 90:
	# 90-degree configurations
	nUnit = 3
	dgap = 2
	meanPOUT = 4.6 # sec
	meanPLIN = 9.7 # sec
	LOT_LENGTH = 10. # ft
	LOT_WIDTH = 17.75 # ft
	LANE_WIDTH = 23. # ft

if angle == 0:
	# 0-degree configurations
	nUnit = 1
	dgap = 1
	meanPOUT = 5.6 # 5.6 sec
	LOT_WIDTH = 10. # ft
	LANE_WIDTH = 12. # ft
	if mode == 'long':
		meanPLIN = 4.8 # sec
		LOT_LENGTH = 32.5 # ft
	if mode == 'short':
		meanPLIN = 16.5 # sec
		LOT_LENGTH = 25.0 # ft
	if mode == 'singlelane':
		meanPLIN = 1.
		LOT_LENGTH = 23.0 # ft

if angle == 45:
	# 45-degree configurations
	nUnit = 2 # ?
	dgap = 2 # ?
	meanPOUT = 10.3 # sec
	meanPLIN = 4.2 # sec
	LOT_LENGTH = 13.0 # ft
	LOT_WIDTH = 17.5 # ft
	LANE_WIDTH = 12.0 # ft

rateDRIV = avgSPEED # ft/sec
CAR_LENGTH = LOT_LENGTH * nUnit # ft
meanDRIV = CAR_LENGTH / rateDRIV # time to drive through one boarding spot in sec
# meanPOUT += meanDRIV
if (angle == 0 and mode == 'short') or (angle == 90 and spmatch):
	meanPLIN += meanDRIV
M_out = int( round(meanPOUT / meanDRIV) )
M_in = int( round(meanPLIN / meanDRIV) ) 
G_out = max(M_out, 3)
# meanPOUT = M_out * meanDRIV
# meanPLIN = M_in * meanDRIV
# rateSERV = 1. / meanSERV
# ratePOUT = 1. / meanPOUT
# ratePLIN = 1. / meanPLIN

if angle == 0:
	''' (vehicle type, spot type)-dependent mean service times '''
	mu_mat = [[meanSERV + 10., meanSERV + 10.],
			  [meanSERV - 5.,  meanSERV]]

	''' (vehicle type, spot type)-dependent mean enter maneuver times'''
	# in_mat = [[M_in,     M_in],
	# 		  [max(2, M_in - 1), M_in]]
	in_mat = [[M_in,     M_in],
			  [M_in - 1, M_in]]

	''' (vehicle type, spot type)-dependent mean exit maneuver times '''
	# out_mat = [[M_out,     M_out],
	# 		   [max(2, M_out - 1), M_out]]
	out_mat = [[M_out,     M_out],
			   [M_out - 1, M_out]]

else: 
	assert angle == 90
	''' (vehicle type, spot type)-dependent mean service times '''
	mu_mat = [[meanSERV + 10., meanSERV + 10., np.nan],
			  [meanSERV - 5.,  meanSERV,       np.nan]]

	''' (vehicle type, spot type)-dependent mean enter maneuver steps '''
	# in_mat = [[M_in,     M_in],
	# 		  [max(2, M_in - 1), M_in]]
	in_mat = [[M_in,     M_in, np.nan],
			  [M_in - 1, M_in, np.nan]]

	''' (vehicle type, spot type)-dependent mean exit maneuver steps '''
	# out_mat = [[M_out,     M_out],
	# 		   [max(2, M_out - 1), M_out]]
	out_mat = [[M_out,     M_out, np.nan],
			   [M_out - 1, M_out, np.nan]]

meanPOUT = None
meanPLIN = None
meanSERV = None 

SEED_SERV = random.Random(2020).sample(range( int(1e12) ), LOT_TOTAL)
SEED_POUT = random.Random(9031).sample(range( int(1e12) ), LOT_TOTAL)
SEED_PLIN = random.Random(4654).sample(range( int(1e12) ), LOT_TOTAL)
SEED_DRIV = random.Random(1203).sample(range( int(1e12) ), LOT_TOTAL)
SEED_ENTR = random.Random(7482).sample(range( int(1e12) ), LOT_TOTAL)

SMALL_INTERVAL = 1e-7

event_priority = {'leave_system':        1,
				  'finish_pulling_out':  2,
				  'prepare_pulling_out': 3,
				  'start_service':       4,
				  'start_second_enter':  5,
				  'prepare_first_exit':  6,
				  'start_pulling_in':    7,
				  'enter_system':        8,
				  'add_stop_idx':        7.5
				  }

#############################################################################
#############################################################################

# the speed limit within parking area is 15 miles per hour i.e. 22 ft/sec
# we assume an approximate value of 23 ft/sec 
# so that for 0-degree (short/normal) configurations it takes 1 sec to drive through the length of one boarding spot

# the data source is the handbook unless noted otherwise
# in general it assumes the design vehicle has a width of 6'7'' and a length of 17'1''

# the data for 0-degree configurations is obtained from the parallel parking stall design
# the parallel stall length is 23'0'' 
# the width projection with minimum level of comfort is 8'3'' 
# with an aisle width of 12' for one-way traffic; or with an aisle width of 22'4'' for two-way traffic
# the width projection with generous level of comfort is 9'0'' 
# with an aisle width of 15'0'' for one-way traffic and 25'4'' for two-way traffic
# our model does not model two-way traffic within the pickup and dropoff facility

# now if we increase the length of each boarding spot,
# it is likely that time for enter and exit maneuvers are likely to decrease 
# thus we also define the parameters for 0-degree (long) configurations
# the data source is from Peter's writeup
# the long parallel stall length is 32'6''

# the data for 90-degree configurations is obtained from the 90-degree parking stall design
# the 90-degree vehicle projection is 17'9''
# and the width projection with minimum level of comfort is 8'3'' with an aisle width of 23' for one-way and two-way traffic
# and the width projection with generous level of comfort is 9'0'' with an aisle width of 26'0'' for one-way and two-way traffic