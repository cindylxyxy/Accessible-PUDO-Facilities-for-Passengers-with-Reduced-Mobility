########################################################################################################################################################
# Last updated: Jul 24, 2022 by Xinyu Liu
# This file defines 2 objects: the system object which stores and updates the system state in a trajectory-based simulation,
# the vehicle object which stores and updates the information related to each vehicle in the system;
# and several helper functions.
# This file assumes for 0-degree configurations with long spots. 
# The other configurations also use the system object in this file as the parent class.
########################################################################################################################################################

import sys
import json
from math import ceil, floor, sqrt
import numpy as np
from heapq import *
from utils import *
from params import *


# a helper function which turns the index of a spot into an equivalent in single-sided configuration
def idx2spot(j):
	if side == 'double':
		return ceil(j / 2)
	else:
		assert side == 'single'
		return j

# i(j) for spot j in {1, ..., N} which returns the lane block next to j
def spot2blk(j):
	if nUnit == 1:
		return idx2spot(j)
	else:
		return ceil( idx2spot(j) / nUnit)

# i_in(j) for spot j in {1, ..., N} which returns the location from which a vehicle starts an enter maneuver into j
def i_in(j):
	if angle == 0 and mode == 'long':
		# i.e. if the vehicle enters a boarding spot by going forward
		return (spot2blk(j) - 1) * CAR_LENGTH
	elif (angle == 0 and mode == 'short') or (angle == 90 and spmatch): 
		# i.e. if the vehicle enters a boarding spot by going backward from 1 lane block ahead
		return spot2blk(j) * CAR_LENGTH
	else:
		assert (angle == 90 and not spmatch)
		# i.e. if the vehicle enters a boarding spot by going backward
		return (idx2spot(j) + dgap) * LOT_LENGTH

# i_out(j) for spot j in {1, ..., N} which returns the location at which a vehicle ends an exit maneuver from j
def i_out(j):
	if (angle == 90 and not spmatch):
		return i_in(j)
	else:
		return spot2blk(j) * CAR_LENGTH

def in_range(j, N):

	if (side == 'double'):
		
		if (angle == 0) and (mode == 'long'):
			if (j % 2 == 1):
				return range( min(N, j + 2), min(N, j + 3) + 1 )
			return range( min(N, j + 1), min(N, j + 2) + 1 )
		
		if (angle == 0) and (mode == 'short'):
			if (j % 2 == 1):
				return range( max(1, j - 2), j + 2 )
			return range( max(1, j - 3), j )

		if (angle == 90):
			if (j % 2 == 1):
				return range( max(1, j - 4), min(N, j + 5) + 1 )
			return range( max(1, j - 5), min(N, j + 4) + 1 )

	assert (side == 'single')
	if (angle == 0) and (mode == 'long'):
		return range( min(N, j + 1), min(N, j + 1) + 1)
	if (angle == 0) and (mode == 'short'):
		return range( max(1, j - 1), j )
	assert (angle == 90)
	return range( max(1, j - 2), min(N, j + 2) + 1 )

def out_range(j, N):
	
	if (side == 'double'):

		if (angle == 0) and (mode == 'long'):
			if (j % 2 == 1):
				return range( max(1, j - 2), j + 2 )
			return range( max(1, j - 3), j )

		if (angle == 0) and (mode == 'short'):
			if (j % 2 == 1):
				return range( max(1, j - 2), j + 2 )
			return range( max(1, j - 3), j )

		if (angle == 90):
			if (j % 2 == 1):
				return range( max(1, j - 4), min(N, j + 5) + 1 )
			return range( max(1, j - 5), min(N, j + 4) + 1 )

	assert side == 'single'

	if (angle == 0):
		return range( max(1, j - 1), j )

	assert (angle == 90)
	return range( max(1, j - 2), min(N, j + 2) + 1 )

class system():

	def __init__(self, N, N_s, pol, seedSERV = None, seedPOUT = None, seedPLIN = None, seedDRIV = None):

		self.make_input(N, N_s, pol, print_input = False)
		# import pdb; pdb.set_trace()

		if simType in ['det', 'cav', 'unif2' ]:
			self.timeDRIV = ParamGen(Cons(rateDRIV, seed = seedDRIV))

		elif simType in ['exp', 'exp1', 'unif', 'unif1']:
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))

		elif simType in ['triag0', 'triag0', 'triag']:
			self.timeDRIV = ParamGen(Disc(rateDRIV, seed = seedDRIV))

		self.timeSERV = [[] for _ in range(N_vt)]
		self.timePOUT = [[] for _ in range(N_vt)]
		self.timePLIN = [[] for _ in range(N_vt)]

		for v in range(N_vt):
			for j in range(len(self.spot_types)):
				if simType in ['cav', 'exp']:
					self.timeSERV[v].append(ParamGen(Expo(self.mu[v][j], seed = seedSERV)))
				elif simType in ['det', 'exp1', 'unif1']:
					self.timeSERV[v].append(ParamGen(Cons(self.mu[v][j], seed = seedSERV)))
				elif simType in ['unif', 'unif2']:
					self.timeSERV[v].append(ParamGen(Unif2(self.mu[v][j], seed = seedSERV)))
				elif simType == 'triag1':
					self.timeSERV[v].append(ParamGen(Tri1(self.mu[v][j], seed = seedSERV)))
				elif simType == 'triag0':
					self.timeSERV[v].append(ParamGen(Tri0(self.mu[v][j], seed = seedSERV)))
				else:
					assert simType == 'triag'
					self.timeSERV[v].append(ParamGen(Tria(self.mu[v][j], seed = seedSERV)))

				if simType in ['exp', 'exp1']:
					self.timePOUT[v].append(ParamGen(Erla(meanDRIV, self.m_out[v][j], seed = seedPOUT)))
					self.timePLIN[v].append(ParamGen(Erla(meanDRIV, self.m_in[v][j], seed = seedPLIN)))
				elif simType in ['det', 'cav', 'unif2']:
					self.timePOUT[v].append(ParamGen(Cons(self.m_out[v][j] * meanDRIV, seed = seedPOUT)))
					self.timePLIN[v].append(ParamGen(Cons(self.m_in[v][j] * meanDRIV, seed = seedPLIN)))
				elif simType in ['unif', 'unif1']:
					self.timePOUT[v].append(ParamGen(Unif2(self.m_out[v][j] * meanDRIV, seed = seedPOUT)))
					self.timePLIN[v].append(ParamGen(Unif2(self.m_in[v][j] * meanDRIV, seed = seedPLIN)))
				else:
					self.timePOUT[v].append(ParamGen(Tria(self.m_out[v][j] * meanDRIV, seed = seedPOUT)))
					self.timePLIN[v].append(ParamGen(Tria(self.m_in[v][j] * meanDRIV, seed = seedPLIN)))

		self.randomVT = ParamGen(Unif(.5))

		self.N = len(self.spot_types)

		if side == 'double':
			self.half_N = N
			if angle == 90:
				self.half_N = int(self.N/2)
				assert self.half_N == N + ceil(N_s/2)

		if angle == 0 and side == 'single':
			assert self.N == N
		elif angle == 0:
			assert side == 'double'
			assert self.N == N*2
		else:
			assert angle == 90

		if angle == 90 and spmatch:
			self.n = spot2blk(self.N) * CAR_LENGTH
		else:
			self.n = idx2spot(self.N) * LOT_LENGTH

		# print ('length of facility:', self.n)

		self.w_0 = w_0
		self.curr = 0.0
		self.start_time = 0.0
		self.eventheap = []
		self.waiting = [[] for _ in range(N_st)]
		self.N_w = 0
		self.inservice = [None for _ in range(self.N)]
		self.head = None
		self.inCount = 0
		self.outCount = 0	
		self.entry_blocked = self.curr
		self.entry_cleared = self.curr
		self.debug = debug
		self.debug_times = []
		self.debug_speed = []
		self.debug_unit = SIM_UNIT

		self.max_idle = 0.0
		self.max_idle_idx = None
		self.max_pct_idle = 0.0
		self.max_pct_idle_idx = None
		self.prod_time = []
		self.idle_time = []
		self.pct_idle = []	
		self.first_service = None
		self.last_departure = - meanDRIV

		self.wait_out = [{'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'oppo_out': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} for _ in range(self.N)]
		
	def make_input(self, N, N_s, pol, print_input = True):

		'''
		N: number of spots on each side
		N_vt: number of vehicle types
		frac_v: fractions of each vehicle type [vector]
		N_st: number of spot types  
		N_s: number of spots per spot type [vector]
		vs_mat: 0 for infeasible, 1 for primary, and 2 for secondary [N_v by N_s matrix]
		mu_mat: (vehicle type, spot type) dependent mean service times
		in_mat: (vehicle type, spot type) dependent mean enter maneuver steps
		out_mat: (vehicle type, spot type) dependent mean exit maneuver steps
		'''
		assert len(frac_v) == N_vt
		assert len(vs_mat) == N_vt
		assert len(vs_mat[0]) == N_st

		if print_input:

			print ('Number of spots on each side:', N)
			print ('Number of vehicle types:', N_vt)
			print ('Fractions of each vehicle type:', frac_v)
			print ('Number of spot types:', N_st)

			priority = {0: 'infeasible', 1: 'primary', 2: 'secondary'}
			for s in range(N_st):
				for v in range(N_vt):
					print ('\tSpot type %s is %s for vehicle type %s, with mean service time %s seconds, m_in = %s steps and m_out = %s steps.'
						%(s, priority[vs_mat[v][s]], v, mu_mat[v][s], in_mat[v][s], out_mat[v][s]))

			print ('Numbers of spots for each type:', N_s)

		''' the spot type of each spot '''
		
		try:
			assert len(N_s) == N_st
			assert np.sum(N_s) == N
		except TypeError:
			assert N_st == 2 or angle == 90
			N_s = [N_s, N - N_s]
			
		self.spot_types = []
		if angle == 0 and side == 'single':
			for s, n_s in enumerate(N_s):
				self.spot_types += [s] * n_s

		elif angle == 0 and side == 'double':
			for s, n_s in enumerate(N_s):
				self.spot_types += [s] * 2 * n_s

		elif angle == 90 and side == 'single':
			self.spot_types += [0, 2, 0] * floor(N_s[0]/2)
			if (N_s[0]%2) == 1:
				self.spot_types += [0, 2]
			self.spot_types += [1] * N_s[1]

		else:
			assert angle == 90 and side == 'double'
			self.spot_types += [0, 0, 2, 2, 0, 0] * floor(N_s[0]/2)
			if (N_s[0]%2) == 1:
				self.spot_types += [0, 0, 2, 2]
			self.spot_types += [1] * N_s[1] * 2

		if place == 'entr':
			pass
		elif place == 'exit':
			self.spot_types = np.flip(self.spot_types)
		else:
			assert place == 'user'	

		# assert len(self.spot_types) == N

		if print_input:
			print ('Side:', side)
			print ('Ascending spot types from:', place)
			print (self.spot_types)

		''' threshold values for each spot type '''
		self.threshold = [0 for _ in range(N_st)]
		for s in range(N_st):
			if pol == 'P':
				try:
					self.threshold[s] = N_s[s]
				except IndexError:
					self.threshold[s] = -1
			elif pol == 'F':
				self.threshold[s] = -1
			else:
				assert type(pol) == int
				self.threshold[s] = pol

		''' type-dependent g_out '''
		gout_mat = np.maximum(out_mat, 3)

		''' (vehicle type, spot)-dependent input '''

		self.mu = [[] for _ in range(N_vt)]
		self.m_in = [[] for _ in range(N_vt)]
		self.m_out = [[] for _ in range(N_vt)]
		self.g_out = [[] for _ in range(N_vt)]

		for v in range(N_vt):
			for s in self.spot_types:
				self.mu[v].append(mu_mat[v][s])
				self.m_in[v].append(in_mat[v][s])
				self.m_out[v].append(out_mat[v][s])
				self.g_out[v].append(gout_mat[v][s])

			# assert (side == 'single' and len(self.mu[v]) == N) or (side == 'double' and len(self.mu[v])== 2*N)
			assert (len(self.mu[v]) == len(self.spot_types))
			
		if print_input:

			for v in range(N_vt):
				print ('Mean service times for v-type %s:'%v, self.mu[v])
				print ('Enter maneuver steps for v-type %s:'%v, self.m_in[v])
				print ('Exit maneuver steps for v-type %s:'%v, self.m_out[v])

	def add_event(self, event):
		heappush( self.eventheap, event) 

	def state_print(self):

		print (self.curr, self.curr_vehicle.j, self.curr_typ)
		x = []
		z = []
		y = [0 for _ in range(spot2blk(self.N))]
		w = [0 for _ in range(spot2blk(self.N))]
		if not (angle == 0 and mode == 'long'):
			y.append(0)
			w.append(0)

		for car in self.inservice:
			if car is None:
				x.append(0)
				z.append(None)
			else:
				z.append(car.type)
				if car.status == 3:
					x.append(1)
				elif car.status in [3.5, 4]:
					x.append(2)
				elif car.status == 2:
					x.append(4)
				else:
					assert car.status == 5
					x.append(3)

		car = self.head
		while car is not None:
			if car.status in [1, 6]:
				car.update_loc()
				if np.abs( car.curr_loc % CAR_LENGTH - CAR_LENGTH ) < 1:
					lane_block = ceil(car.curr_loc / CAR_LENGTH)
				else:
					lane_block = floor( car.curr_loc / CAR_LENGTH )
				if lane_block > len(y):
					print ('one vehicle is leaving ...')
				else:
					if y[lane_block - 1] != 0:
						import pdb; pdb.set_trace()
					if car.status == 1:
						y[lane_block - 1] = car.j
						w[lane_block - 1] = car.type
					else:
						y[lane_block - 1] = self.N + 1
						w[lane_block - 1] = car.type
			car = car.nex

		print ('x:', x)
		print ('z:', z)
		print ('y:', y)
		print ('w:', w)

		if self.N <= 8:
			for idx in range(len(x)):
				if x[idx] == 1:
					print (idx + 1, self.inservice[idx].serv_end)

	def debug_print(self):

		print (self.curr, self.curr_vehicle.j, self.curr_typ)

		if side == 'single':
			print ('Boarding area ...')
			for car in self.inservice:
				if car is not None:
					print ('Vehicle %s is of type %s and assigned spot %s, with current status as %s.' %(car.idx, car.type, car.j, car.status))

		else: 
			assert side == 'double'
			print ('Left boarding area ...')
			for idx in range(self.half_N):
				if self.inservice[2 * idx] is not None:
					car = self.inservice[2 * idx]
					print ('Vehicle %s is of type %s and assigned spot %s, with current status as %s.' %(car.idx, car.type, car.j, car.status))

			print ('Right boarding area ...')
			for idx in range(self.half_N):
				if self.inservice[2 * idx + 1] is not None:
					car = self.inservice[2 * idx + 1]
					print ('Vehicle %s is of type %s and assigned spot %s, with current status as %s.' %(car.idx, car.type, car.j, car.status))

		print ('Through lane ...')
		car = self.head
		while car != None:
			car.update_loc()
			print ('Vehicle %s is of type %s and assigned spot %s, with current status as %s and current location as %s.' %(car.idx, car.type, car.j, car.status, car.curr_loc) )
			car = car.nex			

	def idle_time_calc(self, curr_vehicle):

		total = self.curr - curr_vehicle.enter_time 
		if (total < curr_vehicle.prod_time) and np.abs(total - curr_vehicle.prod_time) > 21 * SMALL_INTERVAL:
			import pdb; pdb.set_trace()
		idle = max(0.0, total - curr_vehicle.prod_time)
		if idle > self.max_idle:
			self.max_idle = idle
			self.max_idle_idx = curr_vehicle.idx
		if idle + curr_vehicle.prod_time == 0.0:
			import pdb; pdb.set_trace()
		if idle / (idle + curr_vehicle.prod_time) > self.max_pct_idle:
			self.max_pct_idle = idle / (idle + curr_vehicle.prod_time)
			self.max_pct_idle_idx = curr_vehicle.idx
		self.prod_time.append(curr_vehicle.prod_time)
		self.idle_time.append(idle)
		self.pct_idle.append( idle / (idle + curr_vehicle.prod_time) )

	def run(self):

		while self.curr - self.start_time <= SIM_UNIT: 

			# if self.curr - self.start_time >= self.debug_unit:
			# 	import pdb; pdb.set_trace()
						
			curr_event = heappop(self.eventheap)
			curr_vehicle = curr_event.vehicle
			curr_typ  = curr_event.typ
			self.curr_vehicle = curr_vehicle
			self.curr = float(curr_event.time)
			self.curr_typ = curr_typ			
			try:
				stop = curr_vehicle.stop
			except AttributeError:
				assert curr_typ == 'enter_system'

			# print (self.curr, curr_vehicle.idx, curr_typ)

			if VEHICLE_IDX is not None and curr_vehicle.idx == VEHICLE_IDX:
				import pdb; pdb.set_trace()

			################################### update system ###################################
			if curr_typ == 'leave_system':
				self.leave_system(curr_vehicle)

			elif curr_typ == 'start_pulling_in':
				self.start_pulling_in(curr_vehicle)

			elif curr_typ == 'start_service':
				self.start_service(curr_vehicle)
			
			elif curr_typ == 'prepare_pulling_out':
				self.prepare_pulling_out(curr_vehicle)

			elif curr_typ == 'finish_pulling_out':
				self.finish_pulling_out(curr_vehicle)

			elif curr_typ == 'add_stop_idx':
				self.add_stop_idx(curr_vehicle)

			else:
				assert curr_typ == 'enter_system'
				self.enter_system(curr_vehicle, debug_idx = VEHICLE_IDX)

			if VEHICLE_IDX is not None and curr_vehicle.idx == VEHICLE_IDX:
				import pdb; pdb.set_trace()

			if self.debug:
				self.state_print()
				if curr_typ == 'enter_system' and curr_vehicle.idx is not None:
					import pdb; pdb.set_trace()				
				'''
				if spmatch:
					self.state_print()
					if curr_typ == 'enter_system' and curr_vehicle.idx is not None:
						import pdb; pdb.set_trace()
				else:
					self.debug_print()
					import pdb; pdb.set_trace()
				'''

		self.start_time = self.curr
		return (self.outCount)

	def get_spots(self, v):

		candidate_spot = []

		for s in range(N_st):
			if vs_mat[v][s] == 1:
				candidate_spot += self.waiting[s]
			elif vs_mat[v][s] == 2 and self.threshold[s] < len(self.waiting[s]):
				candidate_spot += self.waiting[s]
		
		return candidate_spot

	def leave_system(self, curr_vehicle):

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'leave_system') )
			return

		assert self.head == curr_vehicle
		assert curr_vehicle.prev == None					
		self.head = curr_vehicle.nex
		if curr_vehicle.nex != None:
			curr_vehicle.nex.prev = None
			curr_vehicle.nex = None

		self.outCount += 1
		curr_vehicle.prod_time += ((curr_vehicle.dest_to_stop - i_out(curr_vehicle.j)) / curr_vehicle.driv)
		self.idle_time_calc(curr_vehicle)
		self.last_departure = self.curr

		return

	def start_pulling_in(self, curr_vehicle):

		# if self.curr >= 34212. and curr_vehicle.stop == 5:
		# 	import pdb; pdb.set_trace()

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		heap_len = len(self.eventheap)
		event_holder = []
		next_event = heappop(self.eventheap)
		while self.eventheap != [] and next_event.time - self.curr < 100 * SMALL_INTERVAL:
			assert next_event.time >= self.curr
			if next_event.typ == 'start_pulling_in' and next_event.vehicle.stop > curr_vehicle.stop:
				assert next_event.time + 1e-10 > self.curr
				curr_vehicle.end_time = next_event.time + 1e-10
				curr_vehicle.traj = DLinkedList()
				curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, self.curr, 0.0) )
				curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, curr_vehicle.end_time, 'D') )
				for event_visited in event_holder:
					self.add_event( event_visited )
				self.add_event( next_event )
				assert heap_len == len(self.eventheap)
				self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
				car = curr_vehicle.nex
				while car is not None:
					car.update_loc()
					car.update_traj()
					car = car.nex
				return
			heappush(event_holder, next_event)
			next_event = heappop(self.eventheap)

		for event_visited in event_holder:
			self.add_event( event_visited )
		self.add_event( next_event )
		assert heap_len == len(self.eventheap)

		curr_vehicle.curr_loc = curr_vehicle.dest_to_stop
		delayed = False
		req_time = self.curr

		# if the vehicle has arrived at the desired destination while its enter maneuver blocked by vehicles on the through lane
		# i.e. 1) if its immediate proceeding vehicle is having an exit maneuver at the same spot from the middle lane to the through lane;
		# if (curr_vehicle.prev is not None) and (curr_vehicle.prev.j == curr_vehicle.j) and (curr_vehicle.prev.status == 5):
		if (self.inservice[curr_vehicle.j - 1] is not None):
			assert (self.inservice[curr_vehicle.j - 1].status == 5)
			assert (curr_vehicle.prev is not None)
			assert (self.inservice[curr_vehicle.j - 1] == curr_vehicle.prev or side == 'double')
			assert self.inservice[curr_vehicle.j - 1].pout_end > self.curr
			delayed = True 
			req_time = max( req_time, self.inservice[curr_vehicle.j - 1].pout_end )

		elif (curr_vehicle.prev is not None) and (curr_vehicle.prev.stop == curr_vehicle.stop) and (curr_vehicle.prev.status == 5):
			assert side == 'double'
			pass

		# or 2) if the immediate proceeding vehicle has stopped.
		elif (curr_vehicle.prev is not None) and (curr_vehicle.prev.status != 2) and (curr_vehicle.prev.curr_loc <= curr_vehicle.dest_to_stop + 1.5 * CAR_LENGTH):
			cp = curr_vehicle.prev.traj.head
			if cp.data.t > self.curr:
				if cp.data.t <= self.curr + 2 * SMALL_INTERVAL:
					delayed = True
					req_time = max( req_time, cp.data.t )
				else:
					import pdb; pdb.set_trace()
			else:
				while (cp.nex is not None) and (cp.nex.data.t <= self.curr):
					cp = cp.nex
				if cp.nex is not None and cp.nex.data.t <= self.curr + SMALL_INTERVAL:
					cp = cp.nex
				assert cp is not None
				if (cp.data.v == 0.0):
					assert cp.nex.data.t > self.curr
					delayed = True
					req_time = max( req_time, cp.nex.data.t )
				elif (spmatch and cp.data.v != 'D' and cp.data.v < rateDRIV and cp.nex.data.t - meanDRIV > self.curr):
					delayed = True
					req_time = max( req_time, cp.nex.data.t - meanDRIV )

		if delayed:
			assert req_time > self.curr
			curr_vehicle.traj = DLinkedList()
			curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, self.curr, 0.0) )
			curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, req_time, 'D') )
			curr_vehicle.end_time = req_time
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			car = curr_vehicle.nex
			while car is not None:
				car.update_loc()
				car.update_traj()
				car = car.nex
			return

		assert self.inservice[curr_vehicle.j - 1] is None
		self.inservice[curr_vehicle.j - 1] = curr_vehicle
		curr_vehicle.start_in()

		car = curr_vehicle.nex
		while car != None:
			car.update_loc()
			car.update_traj()
			car = car.nex

		assert curr_vehicle.end_time is not None
		self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_service') )
		return

	def start_service(self, curr_vehicle):

		curr_vehicle.start_service()

		if curr_vehicle.prev != None:
			curr_vehicle.prev.nex = curr_vehicle.nex
		else:
			assert self.head == curr_vehicle
			self.head = curr_vehicle.nex

		if curr_vehicle.nex != None:
			car = curr_vehicle.nex
			car.prev = curr_vehicle.prev
			if car.after_plin:
				car.after_plin = False
			while car != None:
				car.update_loc()
				car.update_traj()
				car = car.nex

		curr_vehicle.prev = None
		curr_vehicle.nex = None

		assert curr_vehicle.serv_end > self.curr
		self.add_event( event(curr_vehicle.serv_end, curr_vehicle, 'prepare_pulling_out') )
		return

	def prepare_pulling_out(self, curr_vehicle):

		# if self.curr >= 255. and curr_vehicle.j == 1:
		# 	import pdb; pdb.set_trace()

		stop = curr_vehicle.stop
		assert self.inservice[curr_vehicle.j - 1] == curr_vehicle

		if curr_vehicle.status < 4:
			first_attempt = True
			if curr_vehicle.status == 3:
				assert self.curr == curr_vehicle.serv_end
				self.wait_out[curr_vehicle.j - 1]['veh_count'] += 1
		else:
			first_attempt = False

		if self.first_service is None:
			self.first_service = self.curr
			first_attempt = False

		if curr_vehicle.idx == VEHICLE_IDX:
			delay_reason = None
			delay_status = None
			delay_speed = None
			delayed, req_time, prev, delay_reason, delay_status, delay_speed = self.check_lane_zero_long(curr_vehicle, first_attempt, delay_reason, delay_status, delay_speed)
		else:
			delayed, req_time, prev = self.check_lane_zero_long(curr_vehicle, first_attempt)

		if delayed:
			if not req_time > self.curr:
				import pdb; pdb.set_trace()
			curr_vehicle.status = 4
			curr_vehicle.end_time = req_time
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
			self.wait_out[curr_vehicle.j - 1]['total'] += (req_time - self.curr)
			if curr_vehicle.idx == VEHICLE_IDX:
				self.debug_times.append( tuple( (delay_reason, delay_status, req_time) ) )
				if delay_reason == 5:
					assert delay_speed is not None
					self.debug_speed.append( delay_speed )
			return

		# Now, the current vehicle is ready to start the exit maneuver	
		curr_vehicle.status = 4
		curr_vehicle.prev = prev
		if prev != None:
			curr_vehicle.nex = prev.nex
			try:
				prev.nex.prev = curr_vehicle
			except:
				pass
			prev.nex = curr_vehicle

		elif self.head != None:
			self.head.prev = curr_vehicle
			curr_vehicle.nex = self.head 
			self.head = curr_vehicle

		else:
			self.head = curr_vehicle

		curr_vehicle.start_out()
		
		assert curr_vehicle.pout_end > self.curr
		self.add_event( event( curr_vehicle.pout_end, curr_vehicle, 'finish_pulling_out') )

		# and we update the trajectories of all upcoming vehicles
		if curr_vehicle.nex is not None:
			car = curr_vehicle.nex
			if car.after_plin:
				if not (prev is not None and prev.status == 2):
					import pdb; pdb.set_trace()
				car.after_plin = False
			while car != None:
				car.update_loc()
				car.update_traj()
				car = car.nex

		# lastly we schedule the replacement vehicles
		new_vehicle = vehicle(self, Vtype = self.w_0)
		r = self.randomVT.next()
		for v in range(N_vt):
			if r <= frac_v[v]:
				self.w_0 = v
				break

		self.add_event( event( self.curr, new_vehicle, 'enter_system') )
		self.waiting[self.spot_types[curr_vehicle.j-1]].append(curr_vehicle.j)
		self.N_w += 1
	
		return

	'''
	def add_stop_idx(self, curr_vehicle):
		self.waiting[self.spot_types[curr_vehicle.j-1]].append(curr_vehicle.j)
		self.N_w += 1
	'''

	def finish_pulling_out(self, curr_vehicle):

		assert curr_vehicle.status == 5
		curr_vehicle.status = 6

		assert self.inservice[curr_vehicle.j - 1] is not None
		self.inservice[curr_vehicle.j - 1] = None
		if not (curr_vehicle.end_time is not None and curr_vehicle.end_time > self.curr):
			import pdb; pdb.set_trace()
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'leave_system') )				
		return 

	def enter_system(self, curr_vehicle, debug_idx = None):

		# if self.curr >= 4556.:
		# 	import pdb; pdb.set_trace()

		assert curr_vehicle.status == 0
		if self.entry_blocked == self.curr:
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system' ) )
			return

		J_star = self.get_spots(curr_vehicle.type)

		# if there is no feasible spot for the incoming vehicle type, 
		# we postpone this event until the next event associated with a feasible spot
		if len(J_star) == 0:

			heap_len = len(self.eventheap)
			event_holder = []

			while self.eventheap != []:
				next_event = heappop(self.eventheap)
				event_holder.append(next_event)
				if next_event.typ == 'enter_system':
					continue
				next_vehicle = next_event.vehicle
				assert next_event.time >= self.curr
				if vs_mat[curr_vehicle.type][self.spot_types[next_vehicle.j-1]] > 0:
					break

			for event_visited in event_holder:
				self.add_event( event_visited )
			assert heap_len == len(self.eventheap)

			self.entry_blocked = self.curr
			try:
				self.entry_cleared = next_event.time
			except:
				import pdb; pdb.set_trace()
			self.add_event( event(self.entry_cleared, curr_vehicle, 'enter_system') )
			
			return

		# if there is no vehicle on the driving lane
		# the replacement vehicle is assigned to the spot with the largest index
		if self.head == None:
			assert not spmatch
			self.inCount += 1
			self.head = curr_vehicle
			if free_curb and self.N_w >= 3:
				print (self.N_w, self.waiting)
				import pdb; pdb.set_trace()
			
			if not free_curb:
				j_new = max(J_star)
			else:
				j_new = min(J_star)
			
			assert 1 <= j_new <= self.N
			curr_vehicle.assign_spot(j_new)
			self.N_w -= 1
			self.waiting[self.spot_types[j_new-1]].remove(j_new)

			curr_vehicle.update_traj()
			assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		# if there are vehicles on the driving lane, find the last one
		req_time = self.curr
		car = self.head
		while car.nex is not None:
			car = car.nex	
		car.update_loc()

		# with 0-degree long spots, 
		# if the last one is making an enter maneuver in to spot 2 (single) or spots 3 & 4 (double) 
		# i.e. the loc of the last vehicle is right at the CAR_LENGTH (i.e. first lane block occupied)
		# the replacement vehicle can enter and be assigned to spot 1 under both access control
		if angle == 0 and mode == 'long' and car.curr_loc == CAR_LENGTH and car.status == 2:
			
			assert car.plin_end is not None and car.plin_end >= self.curr
			if car.plin_end - self.curr < 100 * SMALL_INTERVAL:
				self.entry_blocked = self.curr
				self.entry_cleared = car.plin_end
				self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )				
				return

			j_new = None
			if side == 'double' and 2 in J_star and self.curr + meanDRIV * self.m_in[curr_vehicle.type][2-1] > car.plin_end:
				j_new = 2 
			elif 1 in J_star and self.curr + meanDRIV * self.m_in[curr_vehicle.type][1-1] > car.plin_end:
				j_new = 1
			
			if j_new is not None:

				if free_curb and self.N_w >= 3:
					print (self.N_w)
					import pdb;pdb.set_trace()
				
				self.N_w -= 1 
				self.waiting[self.spot_types[j_new-1]].remove(j_new)				
				self.inCount += 1
				curr_vehicle.prev = car
				car.nex = curr_vehicle
				curr_vehicle.assign_spot( j_new )
				assert curr_vehicle.stop == 1
				assert self.inservice[curr_vehicle.j - 1] is None
				
				curr_vehicle.update_traj()
				assert curr_vehicle.end_time != None
				self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
				return

		# if the last one is within CAR_LENGTH (i.e. first lane block occupied)
		# the replacement vehicle cannot enter under either access control
		if car.curr_loc <= CAR_LENGTH + SMALL_INTERVAL:
			car_time = car.calc_time(CAR_LENGTH + SMALL_INTERVAL)
			if car_time == self.curr:
				cp = car.traj.head
				assert cp.data.t <= self.curr
				while cp.nex is not None and cp.nex.data.t <= self.curr:
					cp = cp.nex
				if cp.nex is None:
					import pdb; pdb.set_trace()
				else:
					assert cp.nex.data.t > self.curr
					if cp.data.v > 0.0:
						pass
					else:
						import pdb; pdb.set_trace()

			elif angle == 0 and mode == 'long' and side == 'double' and control == 'full' and (not free_curb) and car.stop == 1 and car.status == 5:
				j_new = None
				if car.j == 1 and 2 in J_star:
					j_new = 2 
				elif car.j == 2 and 1 in J_star:
					j_new = 1
			
				if j_new is not None:

					if free_curb and self.N_w >= 3:
						print (self.N_w)
						import pdb; pdb.set_trace()

					self.N_w -= 1
					self.waiting[self.spot_types[j_new-1]].remove(j_new)
					self.inCount += 1
					curr_vehicle.prev = car
					car.nex = curr_vehicle
					curr_vehicle.assign_spot( j_new )
					assert curr_vehicle.stop == 1
					assert self.inservice[curr_vehicle.j - 1] is None

					curr_vehicle.update_traj()
					assert curr_vehicle.end_time != None
					self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
					return

				else:
					assert car_time > self.curr
					req_time = max( req_time, car_time)	 
			else:
				assert car_time > self.curr
				req_time = max( req_time, car_time)

		if req_time > self.curr:
			self.entry_blocked = self.curr
			self.entry_cleared = req_time
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
			return

		if debug_idx is not None and curr_vehicle.idx == debug_idx:
			import pdb; pdb.set_trace()

		assert control == 'full'
		assert curr_vehicle.curr_loc == 0.0

		# car is the last vehicle on the driving lane 
		# and also the prev for the replacement vehicle if the latter can enter
		J_new = []
		last = car

		# j_new will include the spots that the replacement vehicle can head to 
		# without being blocked or delayed by the last vehicle on the lane in expectation
		for j in sorted(J_star, reverse = (not free_curb)):
			assert j > 0
			J = idx2spot(j)

			if (J < car.stop) or (car.status == 6):
				J_new.append(j)

			elif car.status == 2:
				# K_in with K = car.j and J = idx2spot(j)
				assert car.j != j
				assert car.plin_start <= self.curr
				assert not J == car.stop == 1
				assert car.stop >= 3
				assert car.curr_loc == (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
				if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV >= max(0.0, meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start)) ):
					J_new.append(j)

			elif car.status == 5:
				# K_out with K = car.j and J = idx2spot(j)
				assert car.pout_start <= self.curr <= car.pout_end
				assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
				if side == 'double' and J == car.stop and j != car.j:
					J_new.append(j)
				elif ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= max(0.0, meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start)):
					J_new.append(j)

			else:
				assert car.status == 1
				# I_in with K = car.j and J = idx2spot(j)
				assert car.j != j
				assert car.stop >= 3
				assert car.curr_loc >= curr_vehicle.curr_loc
				if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV >= meanDRIV * self.m_in[car.type][car.j-1] - 100 * SMALL_INTERVAL ):
					J_new.append(j)

		if J_new != []:
			if not free_curb:
				assert J_new[0] == max(J_new)
			else:
				assert J_new[0] == min(J_new)
			car_time = 0.0

		else:
			j = sorted(J_star, reverse = (not free_curb))[0]
			assert (idx2spot(j) >= car.stop)

			if car.status == 2:
				assert car.j != j
				assert car.plin_start <= self.curr
				assert ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start) )
				car_time = meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start) - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV

			elif car.status == 5:
				assert car.pout_start <= self.curr
				assert ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start)  )
				car_time = meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start) - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV

			else:
				assert car.status == 1
				assert car.j != j
				assert ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < meanDRIV * self.m_in[car.type][car.j-1] )
				car_time = meanDRIV * self.m_in[car.type][car.j-1] - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV

		while (J_new != []) and (car.prev is not None):

			car = car.prev
			car.update_loc()
			# print (J_new)
			for idx in range(len(J_new)-1, -1, -1):
				j = J_new[idx]
				# print (idx, j)
				J = idx2spot(j)

				if (J < car.stop) or (car.status == 6):
					pass

				elif car.status == 2:
					# K_in with K = car.j and J = idx2spot(j)
					assert car.j != j
					assert car.plin_start <= self.curr
					assert car.stop >= 3
					assert car.curr_loc == (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
					if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < max(0.0, meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start)) ):
						print (J_new, J_new[:idx])
						J_new = J_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start) - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV)
						break

				elif car.status == 5:
					# K_out with K = car.j and J = idx2spot(j)
					assert car.pout_start <= self.curr <= car.pout_end
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
					if side == 'double' and J == car.stop and j != j:
						pass
					elif ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < max(0.0, meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start)):
						J_new = J_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start) - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
						break

				else:
					assert car.status == 1
					# I_in with K = car.j and J = idx2spot(j)
					assert car.j != j
					assert car.stop >= 3
					assert car.curr_loc >= curr_vehicle.curr_loc
					if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < meanDRIV * self.m_in[car.type][car.j-1] ):
						J_new = J_new[:idx]
						if (idx == 0):
							car_time = max(car_time, meanDRIV * self.m_in[car.type][car.j-1] - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV)
						break

		if J_new == []:
			car_time = self.curr + car_time
			if car_time == self.curr:
				car_time += SMALL_INTERVAL
			if car_time <= self.curr:
				import pdb; pdb.set_trace()
			curr_vehicle.curr_loc = 0.0
			self.add_event( event( car_time, curr_vehicle, 'enter_system') )
			return

		# car.prev is None
		# i.e. there is at least one spot where a replacement vehicle can head to 
		# without being blocked or delayed by any vehicle already on the lane
		# if there are multiple, choose the largest 
		self.inCount += 1
		assert J_new[0] == max(J_new)
		assert J_new[-1] == min(J_new)
		if not free_curb:
			j_new = J_new[0]
		else:
			j_new = J_new[-1]
			
		curr_vehicle.assign_spot( j_new )
		assert j_new in self.waiting[self.spot_types[j_new-1]]
		self.N_w -= 1
		self.waiting[self.spot_types[j_new-1]].remove( j_new )

		curr_vehicle.prev = last
		last.nex = curr_vehicle

		assert curr_vehicle.curr_loc == 0.0
		curr_vehicle.update_traj()
		assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )

		return

	'''
	def check_enter_zero_long(self, enter_time, j):

		assert spmatch and control == 'full'
		car = self.head
		while car is not None:
			if (j <= car.j) or (car.status == 6):
				pass
			elif car.status == 2 and car.stop == 1:
				assert car.plin_start <= self.curr
				enter_time = max(enter_time, car.plin_end + meanDRIV)
			elif car.status == 2:
				assert car.plin_start <= self.curr
				if (car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
					enter_time = max(enter_time, car.plin_start + meanPLIN - (car.stop - 1) * LOT_LENGTH / rateDRIV)
			elif car.status == 5 and car.stop == 1:
				assert car.pout_start <= self.curr <= car.pout_end
				enter_time = max(enter_time, car.pout_end + meanDRIV)
			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				if ((car.stop - 1) * LOT_LENGTH - CAR_LENGTH) / rateDRIV < car.pout_end - self.curr:
					enter_time = max(enter_time, car.pout_end - ((car.stop - 1) * LOT_LENGTH - CAR_LENGTH) / rateDRIV)
			else:
				assert car.status == 1
				assert car.stop >= 2 
				if enter_time < self.curr + meanPLIN - car.curr_loc / rateDRIV - 10 * SMALL_INTERVAL:
					enter_time = self.curr + meanPLIN - car.curr_loc / rateDRIV
			car = car.nex
		return enter_time
	'''

	def check_lane_zero_long(self, curr_vehicle, first_attempt, delay_reason = None, delay_status = None, delay_speed = None, curr_time = None):

		# if self.curr >= 88763. and curr_vehicle.j == 1:
		# 	import pdb; pdb.set_trace()

		assert angle == 0 and mode == 'long'

		delayed = False
		if curr_time is None:
			curr_time = self.curr
			est_pout_start = curr_time
		else:
			assert curr_time > self.curr
			est_pout_start = curr_time
		req_time = curr_time

		car = self.head
		prev = None
		stopped = False

		stop = curr_vehicle.stop
		idx = curr_vehicle.idx

		temp_delay = {'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'oppo_out': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} 
			
		while car != None:

			car.update_loc()
			if car.idx == idx:
				pass

			elif car.curr_loc >= stop * LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL:
				prev = car

			elif car.status == 2 and car.stop == stop + 1:
				if not ((side == 'single' and stop <= self.N - 1) or (side == 'double' and stop <= self.half_N - 1)):
					import pdb; pdb.set_trace()
				if not (self.inservice[car.j - 1] == car):
					import pdb; pdb.set_trace()
				assert car.plin_end >= self.curr > car.plin_start
				if car.plin_end > curr_time:
					delayed = True
					if not car.plin_end > curr_time:
						import pdb; pdb.set_trace()
					if idx == VEHICLE_IDX and car.plin_end > req_time:
						delay_reason = 2
						delay_status = 2
					req_time = max( req_time, car.plin_end)
					temp_delay['front_in'] = max(temp_delay['front_in'], car.plin_end - self.curr)
				else:
					prev = car

			elif car.status == 5: 
				stopped = False
				assert car.stop <= stop
				assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH
				if car.stop == stop and car.pout_end - curr_time > SMALL_INTERVAL:
					delayed = True
					assert car.pout_end > curr_time
					if curr_vehicle.idx == VEHICLE_IDX and car.pout_end > req_time:
						delay_reason = 1 
						delay_status = 5
					req_time = max(req_time, car.pout_end)
					temp_delay['oppo_out'] = car.pout_end - self.curr
					break

				if car.j in out_range(curr_vehicle.j, self.N) and car.pout_end - curr_time > SMALL_INTERVAL:
					delayed = True
					assert car.pout_end > curr_time
					if curr_vehicle.idx == VEHICLE_IDX and car.pout_end > req_time:
						delay_reason = 1 
						delay_status = 5
					req_time = max(req_time, car.pout_end)
					temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car.pout_end - self.curr)

				car_time = car.pout_start + meanDRIV * self.m_out[car.type][car.j-1] + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
				if car_time < curr_time + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1]:
					delayed = True
					car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
					assert car_time > curr_time
					if idx == VEHICLE_IDX and car_time > req_time:
						delay_reason = 3
						delay_status = 5
					req_time = max( req_time, car_time )
					temp_delay['back_out'] = max(temp_delay['back_out'], car_time - self.curr)

			elif car.status == 2:
				assert (car.stop < stop) or (car.stop == stop and side == 'double')
				if car.stop == stop:
					break
				stopped = True

			elif (car.stop == stop) and (car.status == 1):
				assert (car.curr_loc <= car.dest_to_stop == (stop - 1) * LOT_LENGTH)
				stopped = True			 

			elif (car.stop < stop) and (car.status == 1):
				if stop < 2:
					import pdb; pdb.set_trace()
				assert (car.curr_loc <= car.dest_to_stop <= (stop - 2) * LOT_LENGTH)
				stopped = True

			elif (car.stop == stop + 1) and (car.status == 1):
				assert car.dest_to_stop == stop * LOT_LENGTH
				car_time = self.curr + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV

				if stopped and car.prev.status == 1:
					assert car.prev.end_time is not None
					assert car.prev.stop < stop or (car.prev.stop == stop and side == 'double')
					assert (car.prev.stop - 1) * LOT_LENGTH >= car.curr_loc
					car_time = max(meanDRIV * self.m_in[car.prev.type][car.prev.j-1] + car.prev.end_time, self.curr + ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 
				elif stopped:
					assert car.prev.status == 2
					car_time = car.calc_time( (stop - 1) * LOT_LENGTH )

				if car_time < curr_time + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] - 100 * SMALL_INTERVAL:
					delayed = True
					if not car.end_time > curr_time - 3 * SMALL_INTERVAL:
						if curr_time > self.curr:
							delayed = True
							break
						import pdb; pdb.set_trace()
					if car.end_time <= curr_time:
						car_time = car.end_time + SMALL_INTERVAL
						if car_time <= curr_time:
							car_time += SMALL_INTERVAL
					else:
						car_time = car.end_time
					if idx == VEHICLE_IDX and car_time > req_time:
						delay_reason = 4
						delay_status = 1
					req_time = max( req_time, car_time )
					temp_delay['front_drive'] = max(temp_delay['front_drive'], car_time - self.curr)

				if stopped:
					break

			else:
				assert (car.status in [1, 6])
				assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH
				
				if self.curr < car.calc_time(stop * LOT_LENGTH + CAR_LENGTH) < curr_time:
					prev = car

				else:
					car_time = self.curr + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if stopped and car.prev.status == 1:
						assert car.prev.end_time is not None
						assert car.prev.stop < stop or (car.prev.stop == stop and side == 'double')
						assert (car.prev.stop - 1) * LOT_LENGTH >= car.curr_loc
						car_time = max(meanDRIV * self.m_in[car.prev.type][car.prev.j-1] + car.prev.end_time, self.curr + ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 		
					elif stopped:
						assert car.prev.status == 2
						car_time = max(car.prev.plin_end, self.curr + ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 

					if car_time < curr_time + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] - 100 * SMALL_INTERVAL:
						delayed = True
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
						if not car_time > curr_time:
							if control == 'full' and curr_time > self.curr and car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH:
								pass
							else:
								import pdb; pdb.set_trace()
						if idx == VEHICLE_IDX and car_time > req_time:
							delay_reason = 5
							delay_status = car.status
							cp = car.traj.head
							if cp.data.t > curr_time:
								import pdb; pdb.set_trace()
							while cp.nex is not None and cp.nex.data.t <= curr_time:
								cp = cp.nex
							if cp.nex is None:
								import pdb; pdb.set_trace()
							assert cp.nex.data.t > curr_time >= cp.data.t 
							delay_speed = cp.data.v
						req_time = max( req_time, car_time )	
						temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)
				
					if stopped:
						break					

			car = car.nex

		for item in temp_delay:
			self.wait_out[curr_vehicle.j - 1][item] += temp_delay[item] 
		
		if idx == VEHICLE_IDX: 
			return (tuple( (delayed, req_time, prev, delay_reason, delay_status, delay_speed) ))
		else:
			return (tuple( (delayed, req_time, prev) ))

class vehicle():
	
	def __init__(self, sys, Vtype, getin = False, j = None):

		self.sys = sys
		self.driv = self.sys.timeDRIV.next() 
		assert Vtype is not None and Vtype < N_vt
		self.type = Vtype
		self.status = 0

		self.idx = None
		self.j = None
		self.stop = None
		self.block_idx = None
		self.curr_loc = None
		self.dest_to_stop = None
		self.enter_time = None
		self.end_time = None
		self.prod_time = 0.0	

		self.after_plin = False
		self.prev = None
		self.nex = None
		self.traj = None

		self.plin_start = None
		self.plin_time = None
		self.plin_end = None

		self.pout_start = None
		self.pout_time = None
		self.pout_end = None

		self.serv_start = None
		self.serv_time = None
		self.serv_end = None

		if getin:
			assert j is not None
			self.assign_spot(j)
			self.status = 2
			self.start_service()
		else:
			assert j is None
			self.curr_loc = 0.0

	def assign_spot(self, j):
		assert vs_mat[self.type][self.sys.spot_types[j - 1]] > 0
		assert self.status == 0
		self.status = 1
		self.idx = self.sys.inCount
		self.j = j
		self.enter_time = self.sys.curr
		self.stop = idx2spot(j)
		self.block_idx = spot2blk(j)
		self.dest_to_stop = i_in(j)

	def start_in(self):
		assert self.curr_loc == self.dest_to_stop
		assert (self.status == 1)
		self.status = 2
		self.plin_time = self.sys.timePLIN[self.type][self.j-1].next()
		self.plin_end = self.plin_time + self.sys.curr
		self.plin_start = self.sys.curr
		self.prod_time += (self.dest_to_stop / self.driv + self.plin_time)
		self.update_traj()

	def start_service(self):
		assert (self.status == 2)
		self.status = 3
		self.serv_time = self.sys.timeSERV[self.type][self.j-1].next()
		self.serv_end = self.serv_time + self.sys.curr
		self.serv_start = self.sys.curr
		self.prod_time += self.serv_time
		self.traj = None

	def start_out(self):
		assert self.status == 4
		self.status = 5
		self.curr_loc = i_out(self.j)
		self.dest_to_stop = self.sys.n + CAR_LENGTH
		self.pout_time = self.sys.timePOUT[self.type][self.j-1].next()
		self.pout_end = self.pout_time + self.sys.curr
		self.pout_start = self.sys.curr
		self.prod_time += self.pout_time
		self.update_traj()

	def calc_time(self, loc):

		if self.curr_loc > loc:
			if control == 'full':
				pass
			else:
				import pdb; pdb.set_trace()

		cp = self.traj.head

		# if the curr_loc is the target loc,
		# then find the last time that it stays here
		# i.e. the first time that it starts to have nonzero speed.
		# NOTE: nonzero speed can be either positive or 'D'.
		if self.curr_loc == loc:

			if cp.data.t > self.sys.curr:
				if np.abs(cp.data.t - self.sys.curr) > SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				elif cp.nex is not None:
					import pdb; pdb.set_trace()
				else:
					return cp.data. t

			while cp.nex is not None and cp.nex.data.t <= self.sys.curr:
				cp = cp.nex
			assert cp.data.t <= self.sys.curr

			if cp.nex is None:
				return self.sys.curr
			else:
				assert cp.nex.data.t > self.sys.curr
				if cp.data.v > 0.0:
					return self.sys.curr
				while cp.data.v == 0.0:
					cp = cp.nex
				return cp.data.t

		# since all trajectories end with speed 'D' i.e. nonzero
		# we can expect that all vehicles with curr_loc == loc already considered
		# except for those starting with cp.data.v > 0.0,
		# which means the vehicle is arrived at the loc
		while cp.nex != None:
			if cp.nex.data.x <= loc:
				cp = cp.nex
			else:
				break
		# assert cp.data.v != 0.0
		if cp.data.v == 0.0:
			import pdb;pdb.set_trace()

		# now we start with the changepoint s.t.
		# cp.data.x <= loc < cp.nex.data.x
		# Note: if cp.data.v == 0.0, then cp.data.x = cp.nex.data.x
		# given that cp.data.x <= loc
		# we can expect cp.nex.data.x <= loc thus the iteration will go on 
		# which explains the strict inequality between loc and cp.nex.data.x
		if cp.nex == None:
			assert cp.data.v == 'D'
			return cp.data.t

		return cp.data.t + (loc - cp.data.x) / cp.data.v
	
	def update_loc(self):

		assert (self.status != 3) and (self.status != 4)

		if (self.status == 2):
			assert self.curr_loc == self.dest_to_stop
			return

		if (self.status == 5):
			assert self.curr_loc == i_out(self.j)
			return

		try:
			cp = self.traj.head
		except:
			import pdb; pdb.set_trace()

		if cp.nex is None:
			if not np.abs(cp.data.t - self.sys.curr) <= 1e-05:
				import pdb; pdb.set_trace()
			if not cp.data.x == self.curr_loc:
				if np.abs(cp.data.x - self.curr_loc) <= SMALL_INTERVAL:
					self.curr_loc = cp.data.x
				else:
					import pdb; pdb.set_trace()
			return 

		if cp.data.t > self.sys.curr:
			print ('Line 474: current time not accounted for in the trajectory!!!')
			import pdb; pdb.set_trace()

		if self.end_time < self.sys.curr - 2 * SMALL_INTERVAL:
			print ('Line 478: self.end_time, self.sys.curr', self.end_time, self.sys.curr)
			import pdb; pdb.set_trace()

		while cp.nex.data.t <= self.sys.curr:
			cp = cp.nex
			if cp.nex is None:
				if np.abs(self.end_time - self.sys.curr) > 3 * SMALL_INTERVAL:
					print ('Line 486: it should be small', np.abs(self.end_time - self.sys.curr))
					import pdb; pdb.set_trace()
				self.curr_loc = cp.data.x
				return

		assert cp.data.t <= self.sys.curr < cp.nex.data.t
		assert cp.data.v != 'D'
		self.curr_loc = cp.data.x + cp.data.v * (self.sys.curr - cp.data.t)
		if self.curr_loc > self.dest_to_stop:
			if not (self.curr_loc - self.dest_to_stop < 1.5e-04):
				import pdb; pdb.set_trace()
			self.curr_loc = self.dest_to_stop
		return

	def update_traj(self):

		# if self.idx == 48370 and self.sys.curr >= 633321.:
		# 	print (self.sys.curr)
		# 	if self.traj is not None:
		# 		self.traj.print()
		# 	import pdb; pdb.set_trace()

		# if self.sys.curr >= 12202. and self.idx == 1057:
		# 	import pdb; pdb.set_trace()

		traj = self.traj
		end_time = self.end_time
		assert self.status in [1, 2, 5, 6]

		self.traj = DLinkedList()

		if self.status == 2:
			# i.e. the vehicle has started pulling in but has not started service
			assert (self.curr_loc == self.dest_to_stop)
			if self.plin_end > self.sys.curr:
				self.end_time = self.plin_end
			else:
				assert self.end_time == self.plin_end
			if end_time is not None:
				assert self.end_time >= end_time
			assert self.end_time >= self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0) )
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		if self.curr_loc == self.dest_to_stop and self.status == 6:
			assert self.end_time is not None
			assert (np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL)
			if end_time is not None:
				assert self.end_time >= end_time
			assert self.end_time >= self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		if self.curr_loc == self.dest_to_stop and self.status == 1:
			if self.end_time is None:
				assert (self.stop == 1) or (angle == 90 and self.block_idx == 1)
				self.end_time = self.sys.curr
			if end_time is not None:
				assert self.end_time >= end_time
			if not self.end_time >= self.sys.curr:
				if not np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		start_t = self.sys.curr

		if self.status == 5:
			# i.e. the vehicle has started pulling out
			self.update_loc()
			assert self.traj.head == None
			if self.pout_end > self.sys.curr:
				self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0) )
				start_t = self.pout_end
			
		if self.prev is None:
			try:
				assert self.sys.head == self
			except:
				print ('Line 541: this car does not have a prev while self.sys.head is not itself!!!')
				import pdb; pdb.set_trace()

			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if not (self.end_time >= end_time - 70 * SMALL_INTERVAL):
					import pdb; pdb.set_trace()
				self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return

		if self.prev.end_time < start_t + 100 * SMALL_INTERVAL:
			# if self.prev.status = 1, starting the enter maneuver soon => self can travel up to self.prev.dest_to_stop - CAR_LENGTH
			# if self.prev.status = 2, finishing the enter maneuver soon => prev can be ignored
			# if self.prev.status = 5 or 6, leaving the system soon => prev can be ignored 
			if self.prev.status == 1 and self.prev.end_time > start_t:
				# import pdb; pdb.set_trace()
				if self.curr_loc >= self.prev.dest_to_stop - CAR_LENGTH:
					self.curr_loc = self.prev.dest_to_stop - CAR_LENGTH
					self.traj.addEnd(  changePoint(self.curr_loc, start_t, 0.0) )
					start_t = self.prev.end_time 
			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if not (self.end_time >= end_time - 100 * SMALL_INTERVAL):
					import pdb; pdb.set_trace()
				self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return

		if self.prev.status == 2 and self.dest_to_stop >= self.prev.curr_loc:

			if angle == 0 and mode == 'long':
				assert self.prev.curr_loc == (self.prev.stop - 1) * LOT_LENGTH

				if not self.after_plin:
					if not self.curr_loc <= self.prev.curr_loc - CAR_LENGTH + 1e-04:
						if self.sys.curr_typ == 'start_pulling_in' and self.sys.curr_vehicle == self.prev:
							print (self.prev.curr_loc - self.curr_loc - CAR_LENGTH)
							import pdb; pdb.set_trace()
							self.curr_loc = self.prev.curr_loc - CAR_LENGTH
						else:
							import pdb; pdb.set_trace()
					self.after_plin = True

					if start_t + (self.prev.curr_loc - self.curr_loc) / self.driv >= self.prev.plin_end:
						self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
						self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
						if end_time is not None and self.end_time < end_time:
							if not self.end_time >= end_time - 62 * SMALL_INTERVAL:
								import pdb; pdb.set_trace()
							self.end_time = end_time
						self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
					else:
						if self.curr_loc < self.prev.curr_loc - CAR_LENGTH:
							self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
							start_t += (self.prev.curr_loc - CAR_LENGTH - self.curr_loc) / self.driv
						new_speed = CAR_LENGTH / (self.prev.plin_end - start_t)
						self.traj.addEnd( changePoint(self.prev.curr_loc - CAR_LENGTH, start_t, new_speed) )
						self.traj.addEnd( changePoint(self.prev.curr_loc, self.prev.plin_end, self.driv))
						self.end_time = self.prev.plin_end + (self.dest_to_stop - self.prev.curr_loc) / self.driv
						if end_time is not None and self.end_time < end_time:
							if not self.end_time >= end_time - SMALL_INTERVAL:
								import pdb; pdb.set_trace()
							self.end_time = end_time
						self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )

				else:
					assert traj is not None
					assert end_time >= self.prev.plin_end
					cp = traj.head
					while cp.nex is not None and cp.nex.data.t <= self.sys.curr:
						cp = cp.nex
					while cp is not None:
						self.traj.addEnd( cp.data )
						cp = cp.nex
					assert self.end_time is not None
				return					

			else:
				if self.status != 5 and np.abs( self.curr_loc + CAR_LENGTH - self.prev.curr_loc ) < SMALL_INTERVAL:
					self.curr_loc = self.prev.curr_loc - CAR_LENGTH

		#####################################################################################
		cp = self.prev.traj.head
		if cp.nex is None:
			import pdb; pdb.set_trace()

		# the assertion above should be fine since cp.nex == None iff either one of the following holds:
		# i) the vehicle is pulling in and thus self.getin is True
		# ii) the vehicle is the head onlane thus self.prev is None
		# it implies that cp.data.v != 'D'
		assert (cp.data.t <= start_t)
		while (cp.nex is not None) and (cp.nex.data.t <= start_t):
			assert cp.data.v != 'D'
			cp = cp.nex

		# the while-loop cannot end with (cp.nex is None)
		# remember that self.prev should already has an accurate trajectory w. a valid end_time
		# if (cp.nex is None), then (self.prev.end_time == cp.nex.data.t)
		# since (self.prev.end_time <= start_t) is already accounted for
		# it is only possible to have (cp.nex.data.t > start_t) and also (cp.nex.data.v = 'D')
		if (cp.nex is None):
			print ('line 575: testing if an option is possible')
			import pdb; pdb.set_trace()
			assert (self.prev.end_time == cp.nex.data.t > start_t)
			assert cp.data.v == 'D'
			# ??? assert self.prev.status == 5
			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if end_time is not None:
				assert self.end_time >= end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return
		
		# some sanity check before proceeding
		assert (cp.data.v != 'D') and (cp.nex != None)
		assert (cp.data.t <= start_t < cp.nex.data.t)

		# let us define some counters to keep track of time and location (in the projected traj)
		iter_t = start_t
		if (start_t == self.sys.curr):
			self.prev.update_loc()
			prev_loc = self.prev.curr_loc
		else:
			assert self.status == 5
			prev_loc = cp.data.x + (start_t - cp.data.t) * cp.data.v
			if angle == 90 and 0 <= self.curr_loc + CAR_LENGTH - prev_loc < 1e-12:
				prev_loc = self.curr_loc + CAR_LENGTH
				if cp.data.x < cp.nex.data.x == prev_loc and cp.data.v == 0.0:
					curr = self.prev
					while curr.prev is not None and curr.prev.status != 2:
						curr = curr.prev
					while curr.idx != self.idx:
						curr.update_traj()
						curr = curr.nex
					# import pdb; pdb.set_trace()
					cp = self.prev.traj.head
					assert cp.nex is not None
					assert (cp.data.t <= start_t)
					while (cp.nex is not None) and (cp.nex.data.t <= start_t):
						assert cp.data.v != 'D'
						cp = cp.nex
					assert cp.nex is not None
					assert (cp.data.v != 'D') and (cp.nex != None)
					assert (cp.data.t <= start_t < cp.nex.data.t)	
					if cp.data.x + (start_t - cp.data.t) * cp.data.v != prev_loc:
						import pdb; pdb.set_trace()				

		# the RHS below is the location of prev at start_t
		# check that the distance between the prev and the curr at start_t is at least CAR_LENGTH
		if (self.status != 5) and 0 <= self.curr_loc + CAR_LENGTH - prev_loc < 15e-5:
			self.curr_loc = max( 0.0, prev_loc - CAR_LENGTH )
			if (prev_loc < CAR_LENGTH):
				self.traj.addEnd( changePoint(0.0, start_t, 0.0) )
				start_t = self.prev.calc_time(CAR_LENGTH + SMALL_INTERVAL)
				assert start_t > self.sys.curr
		start_x = self.curr_loc
		iter_x = self.curr_loc

		# check that the distance between the prev and the curr at start_t is at least CAR_LENGTH
		if self.curr_loc + CAR_LENGTH > prev_loc:
			import pdb; pdb.set_trace()

		if self.curr_loc == self.dest_to_stop and self.status == 1:
			if self.end_time is None:
				assert (self.stop == 1)
				self.end_time = self.sys.curr
			if end_time is not None:
				assert self.end_time >= end_time
			if not self.end_time >= self.sys.curr:
				if not np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		# if the distance == CAR_LENGTH AND cp.data.v == self.driv
		# the curr has to travel at cp.data.v == self.driv for some time
		# i.e. the curr already catches up with the its prev
		#      and the curr has a free-flow rate == cp.data.v s.t the headway distance is maintained 
		if (start_x + CAR_LENGTH == prev_loc) and (cp.data.v == self.driv):
			
			while (cp.nex is not None):
				cp = cp.nex
				if (cp.data.v == 'D') or (cp.data.x >= self.dest_to_stop + CAR_LENGTH):
					self.traj.addEnd( changePoint(start_x, start_t, self.driv) )
					self.end_time = start_t + (self.dest_to_stop - start_x) / self.driv
					if (end_time is not None) and (self.end_time < end_time):
						if not (self.end_time >= end_time - 100 * SMALL_INTERVAL):
							import pdb; pdb.set_trace()
						self.end_time = end_time
					self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
					return
				elif (cp.data.v < self.driv):
					self.traj.addEnd( changePoint(start_x, start_t, self.driv) )
					break
				elif (cp.data.v > self.driv):
					break
				else:
					assert (cp.data.v == self.driv)

			assert (cp is not None) and (cp.data.v != self.driv)
			if not (start_x + (cp.data.t - start_t) * self.driv + CAR_LENGTH == cp.data.x):
				if np.abs(start_x + (cp.data.t - start_t) * self.driv + CAR_LENGTH - cp.data.x) > 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
			iter_t = cp.data.t
			prev_loc = cp.data.x
			iter_x = prev_loc - CAR_LENGTH
		
		while True:

			# if the distance == CAR_LENGTH AND cp.data.v < self.driv
			# the curr has to travel at cp.data.v <= self.driv for some time
			# i.e. the curr already catches up with the its prev
			#      and the curr has a free-flow rate geq cp.data.v s.t the headway distance is maintained 
			# the curr starts to travel at cp.data.v (not constant, but always the same as the prev)

			if (iter_x + CAR_LENGTH == prev_loc) and (cp.data.v < self.driv):

				while (cp.nex is not None):

					assert cp.data.v != 'D' 
					if not (iter_x + CAR_LENGTH == cp.data.x + (iter_t - cp.data.t) * cp.data.v):
						import pdb; pdb.set_trace()

					if cp.data.v <= self.driv:

						if cp.data.v == 0 and np.abs(iter_x - self.dest_to_stop) < SMALL_INTERVAL:
							self.end_time = iter_t
							if (end_time is not None) and (self.end_time < end_time):
								if not (self.end_time >= end_time - 21 * SMALL_INTERVAL):
									import pdb; pdb.set_trace()
								self.end_time = end_time						
							self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
							return
						else: 
							self.traj.addEnd( changePoint(iter_x, iter_t, cp.data.v) )

						if cp.nex.data.x - CAR_LENGTH >= self.dest_to_stop:
							if not (cp.data.v > 0.0):
								if ( np.abs(cp.nex.data.x - cp.data.x) < SMALL_INTERVAL ):
									self.end_time = cp.nex.data.t
									if (end_time is not None) and (self.end_time < end_time):
										if not (self.end_time >= end_time - SMALL_INTERVAL):
											import pdb; pdb.set_trace()
										self.end_time = end_time
									self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
								else:
									import pdb; pdb.set_trace()
							else:
								self.end_time = iter_t + (self.dest_to_stop - iter_x) / cp.data.v
								if (end_time is not None) and (self.end_time < end_time):
									if not (self.end_time >= end_time - 92 * SMALL_INTERVAL):
										import pdb; pdb.set_trace()
									self.end_time = end_time
								self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
							return

						if np.abs(iter_x + (cp.nex.data.t - iter_t) * cp.data.v + CAR_LENGTH - cp.nex.data.x) > 1.4e-04:
							import pdb; pdb.set_trace()
						cp = cp.nex 
						iter_x = cp.data.x - CAR_LENGTH
						iter_t = cp.data.t
					
					else:
						assert cp.data.v > self.driv
						break	

				if cp.nex is None:
					assert cp.data.v == 'D'
					assert (iter_x + CAR_LENGTH == cp.data.x)
					self.traj.addEnd( changePoint(iter_x, iter_t, self.driv) )
					self.end_time = iter_t + (self.dest_to_stop - iter_x) / self.driv
					if (end_time is not None) and (self.end_time < end_time):
						if not (self.end_time >= end_time - 99 * SMALL_INTERVAL):
							import pdb; pdb.set_trace()
						self.end_time = end_time
					self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
					return

				assert cp.nex is not None
				assert (cp.data.v != 'D') and (cp.data.v > self.driv)
				assert (iter_x + CAR_LENGTH == cp.data.x)
				assert (iter_x < self.dest_to_stop)
				start_t = iter_t
				start_x = iter_x
				prev_loc = cp.data.x

			# if the distance is > CAR_LENGTH 
			# OR if the distance == CAR_LENGTH while self.driv < cp.data.v
			# the curr can drive at its own free-flow rate for some time
			if not ( (iter_x + CAR_LENGTH < prev_loc) | (cp.data.v > self.driv) ):
				import pdb; pdb.set_trace()
			assert (iter_x + CAR_LENGTH <= prev_loc)

			# the curr starts to drive at self.driv
			self.traj.addEnd( changePoint(start_x, start_t, self.driv) )

			while (cp.nex is not None):

				if (iter_x >= self.dest_to_stop - SMALL_INTERVAL):
					break

				# as long as self.driv <= cp.data.v,
				# the curr continues to travel at self.driv
				# the inequality is not strict 
				# because the distance between two vehicles (> CAR_LENGTH) will be maintained if both travels at the same rate
				if (self.driv <= cp.data.v):

					if (self.driv == cp.data.v) and (iter_x + CAR_LENGTH == prev_loc):
						cp = cp.nex
						prev_loc = cp.data.x
						iter_x = prev_loc - CAR_LENGTH
						iter_t = cp.data.t
						if (cp.data.v != 'D') and (cp.data.v < self.driv):
							break
					else:
						cp = cp.nex
						prev_loc = cp.data.x
						iter_x = start_x + (cp.data.t - start_t) * self.driv
						iter_t = cp.data.t
						if not (iter_x + CAR_LENGTH <= prev_loc):
							if np.abs(iter_x + CAR_LENGTH - prev_loc) < 15e-05:
								iter_x = prev_loc - CAR_LENGTH
							else:
								import pdb; pdb.set_trace()
					continue

				assert (cp.data.v < self.driv)
				if (iter_x + CAR_LENGTH == prev_loc):
					break

				if np.abs(iter_x + (cp.nex.data.t - iter_t) * self.driv + CAR_LENGTH - cp.nex.data.x) < SMALL_INTERVAL:
					cp = cp.nex
					iter_x = cp.data.x - CAR_LENGTH
					iter_t = cp.data.t
					prev_loc = cp.data.x
					if (cp.data.v == 'D') or (cp.data.v < self.driv):
						break

				elif (iter_x + (cp.nex.data.t - iter_t) * self.driv + CAR_LENGTH < cp.nex.data.x):
					cp = cp.nex
					iter_x = start_x + (cp.data.t - start_t) * self.driv
					iter_t = cp.data.t
					prev_loc = cp.data.x
					assert (iter_x + CAR_LENGTH < cp.data.x)

				else:
					change_t = (prev_loc - CAR_LENGTH - iter_x) / (self.driv - cp.data.v)
					if not (0.0 < change_t < cp.nex.data.t - cp.data.t):
						import pdb; pdb.set_trace()
					if cp.data.v == 0.0:
						iter_x = prev_loc - CAR_LENGTH
						iter_t += change_t
					else:
						iter_t += change_t
						prev_loc = cp.data.x + (iter_t - cp.data.t) * cp.data.v
						iter_x  = prev_loc - CAR_LENGTH
						
					break

			# if the prev ends before the curr catches up with the prev
			# i.e. the curr can travel at self.driv until its own destination (from the info available up to now)
			# OR if the curr catches up with the prev after the destination of the curr
			if (cp.nex is None) | (iter_x >= self.dest_to_stop - SMALL_INTERVAL): 
				self.end_time = start_t + (self.dest_to_stop - start_x) / self.driv
				if (end_time is not None) and (self.end_time < end_time):
					if not (self.end_time >= end_time - 100 * SMALL_INTERVAL):
						import pdb; pdb.set_trace()
					self.end_time = end_time
				self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
				return

			assert (cp.data.v != 'D')
			assert (cp.data.v < self.driv)
			start_x = iter_x
			start_t = iter_t
			if not (prev_loc == start_x + CAR_LENGTH):
				import pdb; pdb.set_trace()

		return
