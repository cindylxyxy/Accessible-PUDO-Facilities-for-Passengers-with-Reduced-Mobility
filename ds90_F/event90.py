########################################################################################################################################################
# Last updated: Jul 24, 2022 by Xinyu Liu
# This file defines 2 objects: the system object which stores and updates the system state in a trajectory-based simulation,
# the vehicle object which stores and updates the information related to each vehicle in the system;
# and several helper functions.
# This file deals specifically with 90-degree configurations.
########################################################################################################################################################

import sys
import json
from math import ceil, floor, sqrt
import numpy as np
from heapq import *
from utils import *
from params import *
from event import *


class system90(system):

	def __init__(self, N, N_s, pol, seedSERV = None, seedPOUT = None, seedPLIN = None, seedDRIV = None):
		
		super(system90, self).__init__(N, N_s, pol, seedSERV, seedPOUT, seedPLIN, seedDRIV)
		
		self.wait_out = [{'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'nearby_in': 0.0, 'oppo_in': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} for _ in range(self.N)]

	def start_pulling_in(self, curr_vehicle):

		# if self.curr >= 29894. and curr_vehicle.stop == 1:
		# 	import pdb; pdb.set_trace()

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		if self.eventheap != []:
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

		if not self.inservice[curr_vehicle.j - 1] is None:
			import pdb; pdb.set_trace()
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

	def prepare_pulling_out(self, curr_vehicle):

		# if self.curr >= 29894. and curr_vehicle.j == 1:
		# 	import pdb; pdb.set_trace()

		stop = curr_vehicle.stop
		block_idx = curr_vehicle.block_idx
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
			delayed, req_time, prev, delay_reason, delay_status, delay_speed = self.check_lane_90(curr_vehicle, first_attempt, delay_reason, delay_status, delay_speed)
		else:
			delayed, req_time, prev = self.check_lane_90(curr_vehicle, first_attempt)

		if delayed:
			if not req_time > self.curr:
				import pdb; pdb.set_trace()
			curr_vehicle.status = 4
			curr_vehicle.end_time = req_time
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )	
			self.wait_out[curr_vehicle.j - 1]['total'] += (req_time - self.curr)
			return

		if curr_vehicle.status == 4 and self.eventheap != []:
			broken = False
			heap_len = len(self.eventheap)
			event_holder = []
			while self.eventheap != []:
				next_event = heappop(self.eventheap)
				heappush(event_holder, next_event)
				if next_event.time > self.curr + 1e-05:
					break 
				if next_event.typ == 'prepare_pulling_out' and curr_vehicle.j < next_event.vehicle.j:
					if next_event.time <= self.curr:
						import pdb; pdb.set_trace()
					curr_vehicle.end_time = next_event.time
					self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
					self.wait_out[curr_vehicle.j - 1]['total'] += (curr_vehicle.end_time - self.curr)
					broken = True
					break

			for event_visited in event_holder:
				self.add_event( event_visited )	
			if broken:
				assert heap_len + 1 == len(self.eventheap)
				return
			assert heap_len == len(self.eventheap)

		est_pout_start = self.curr

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

	def enter_system(self, curr_vehicle, debug_idx = None):

		# if self.curr >= 93.:
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
			self.entry_cleared = next_event.time
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

			if (spot2blk(j) < car.block_idx - 1) or (car.status == 6):
				J_new.append(j)

			elif car.status == 2:
				# K_in with k = car.j and j
				assert car.j != j
				assert car.plin_start <= self.curr
				assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
				if (car.stop - 1) * LOT_LENGTH / rateDRIV >= max(0.0, meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start)):
					J_new.append(j)

			elif car.status == 5:
				# K_out with k = car.j and j
				assert car.pout_start <= self.curr <= car.pout_end
				assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
				if ( (car.stop - 1) * LOT_LENGTH / rateDRIV >= max(0.0, meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start) - 7e-06) ):
					J_new.append(j)

			else:
				assert car.status == 1
				# I_in with k = car.j and j
				assert car.j != j
				assert car.curr_loc >= CAR_LENGTH + curr_vehicle.curr_loc - SMALL_INTERVAL
				if (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= meanDRIV * self.m_in[car.type][car.j-1]:
					J_new.append(j)

		if J_new != []:
			if not free_curb:
				assert J_new[0] == max(J_new)
			else:
				assert J_new[0] == min(J_new)
			car_time = 0.0

		else:
			j = sorted(J_star, reverse = (not free_curb))[0]
			# assert (idx2spot(j) >= car.stop)

			if car.status == 2:
				assert car.j != j
				assert car.plin_start <= self.curr
				assert (car.stop - 1) * LOT_LENGTH / rateDRIV < meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start)
				car_time = meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV 

			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				assert (car.stop - 1) * LOT_LENGTH / rateDRIV < meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start)
				car_time = meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV

			else:
				assert car.status == 1
				assert car.j != j
				assert (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanDRIV * self.m_in[car.type][car.j-1]
				car_time = meanDRIV * self.m_in[car.type][car.j-1] - (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV

		while (J_new != []) and (car.prev is not None):

			car = car.prev
			car.update_loc()

			for idx in range(len(J_new)-1, -1, -1):
				j = J_new[idx]
				# print (idx, j)
				J = idx2spot(j)

				if (spot2blk(j) < car.block_idx - 1) or (car.status == 6):
					pass

				elif car.status == 2:
					# K_in with K = car.j and J = idx2spot(j)
					assert car.j != j
					assert car.plin_start <= self.curr
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
					if (car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start)):
						J_new = J_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanDRIV * self.m_in[car.type][car.j-1] - (self.curr - car.plin_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV)
						break						

				elif car.status == 5:
					# K_out with K = car.j and J = idx2spot(j)
					assert car.pout_start <= self.curr <= car.pout_end
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
					if (car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start) - 3e-06):
						J_new = J_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanDRIV * self.m_out[car.type][car.j-1] - (self.curr - car.pout_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV)
						break

				else:
					assert car.status == 1
					# I_in with K = car.j and J = idx2spot(j)
					assert car.j != j
					assert car.curr_loc >= CAR_LENGTH + curr_vehicle.curr_loc
					if (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanDRIV * self.m_in[car.type][car.j-1]:
						J_new = J_new[:idx]
						if (idx == 0):
							car_time = max(car_time, meanDRIV * self.m_in[car.type][car.j-1] - (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
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
		
		assert curr_vehicle.dest_to_stop >= curr_vehicle.curr_loc
		curr_vehicle.update_traj()
		assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
		
		return

	def check_lane_90(self, curr_vehicle, first_attempt, delay_reason = None, delay_status = None, delay_speed = None, curr_time = None):

		# if self.curr >= 774. and curr_vehicle.j == 1:
		# 	import pdb; pdb.set_trace()

		assert angle == 90

		delayed = False
		if curr_time is None:
			curr_time = self.curr
			est_pout_start = curr_time
			est_pout_end = est_pout_start + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1]
		else:
			if not curr_time > self.curr:
				import pdb; pdb.set_trace()
			est_pout_start = curr_time
			est_pout_end = est_pout_start + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1]
		req_time = curr_time

		car = self.head
		prev = None
		stopped = False
		
		stop = curr_vehicle.stop
		idx = curr_vehicle.idx
		block_idx = curr_vehicle.block_idx

		temp_delay = {'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'nearby_in': 0.0, 'oppo_in': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} 

		while car != None:

			car.update_loc()
			if car.idx == idx:
				pass

			elif car.status == 2:
				assert car.plin_start <= self.curr <= car.plin_end
				if car.stop == stop and est_pout_start < car.plin_end - 10 * SMALL_INTERVAL:
					assert (side == 'double' and car.j != curr_vehicle.j)
					stopped = False
					delayed = True
					req_time = max( req_time, car.plin_end + SMALL_INTERVAL )
					temp_delay['oppo_in'] = max(temp_delay['oppo_in'], car.plin_end + SMALL_INTERVAL - self.curr)
				elif stop + dgap >= car.stop >= stop - dgap and est_pout_start < car.plin_end - 1e-05:
					stopped = False
					delayed = True
					req_time = max( req_time, car.plin_end + SMALL_INTERVAL )
					temp_delay['nearby_in'] = max(temp_delay['nearby_in'], car.plin_end + SMALL_INTERVAL - self.curr)
				elif car.stop < stop - dgap and car.plin_end > curr_time:
					stopped = True
				elif car.block_idx > block_idx and curr_time <= car.plin_end:
					if car.block_idx == block_idx + 1 and (car.stop > stop + dgap) and curr_time < car.plin_end - meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1]:
						delayed = True
						req_time = max( req_time, car.plin_end - meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] )
						temp_delay['front_in'] = max(temp_delay['front_in'], car.plin_end - meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] - self.curr)

					elif control == 'full' or ptype == 0:
						if est_pout_end + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV < max(curr_time, car.plin_start + meanDRIV * self.m_in[car.type][car.j-1]) - 2e-05:
							assert car.plin_start + meanDRIV * self.m_in[car.type][car.j-1] > curr_time
							assert car.plin_start + meanDRIV * self.m_in[car.type][car.j-1] > est_pout_end + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV
							if not car.plin_start + meanDRIV * self.m_in[car.type][car.j-1] > est_pout_start + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV:
								import pdb; pdb.set_trace()
							else:	
								delayed = True
								car_time = meanDRIV * self.m_in[car.type][car.j-1] + car.plin_start - (meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV)
								req_time = max( req_time, car_time )
								temp_delay['front_in'] = max(temp_delay['front_in'], car_time - self.curr)
						else:	
							prev = car
					else:
						prev = car				

			elif car.j in out_range(curr_vehicle.j, self.N) and car.status == 5 and car.pout_end - curr_time > SMALL_INTERVAL:
				assert car.pout_end > curr_time
				delayed = True
				req_time = max(req_time, car.pout_end + meanDRIV)
				temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car.pout_end + meanDRIV - self.curr)
				break

			elif car.curr_loc >= (stop + dgap) * LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL:

				if (control == 'full' or ptype == 0) and car.stop > (stop + dgap) and car.status == 1:
					car_time = est_pout_end + (car.stop - 1 - stop - dgap) * LOT_LENGTH / rateDRIV 
					if car_time < car.end_time + meanDRIV * self.m_in[car.type][car.j-1] - 2e-05:
						if not car.end_time + meanDRIV * self.m_in[car.type][car.j-1] > est_pout_start + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] + (car.stop - 1 - stop - dgap) * LOT_LENGTH / rateDRIV:
							import pdb; pdb.set_trace()
						else:
							delayed = True
							car_time = car.end_time + meanDRIV * self.m_in[car.type][car.j-1] - meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] - (car.stop - 1 - stop - dgap) * LOT_LENGTH / rateDRIV
							req_time = max( req_time, car_time )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)
					else:
						prev = car
				else:
					prev = car

			elif car.status == 5:
				stopped = False
				assert car.stop < stop - dgap or car.pout_end - curr_time <= SMALL_INTERVAL
				assert car.pout_start <= self.curr
				assert (stop - 1) * LOT_LENGTH >= car.curr_loc or car.pout_end - curr_time <= SMALL_INTERVAL		
				if not car.dest_to_stop >= min( (stop + dgap) * LOT_LENGTH + CAR_LENGTH, self.n + CAR_LENGTH):
					import pdb; pdb.set_trace()
				car_time = max(self.curr, meanDRIV * self.m_out[car.type][car.j-1] + car.pout_start) + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
				if car_time < curr_time + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1]:
					car_time = car.calc_time( (stop + dgap) * LOT_LENGTH + CAR_LENGTH )
					assert car_time > self.curr
					delayed = True
					req_time = max( req_time, car_time )
					temp_delay['back_out'] = max(temp_delay['back_out'], car_time - self.curr)

			elif car.status == 1 and car.j == curr_vehicle.j:
				assert car.stop == stop
				assert (car.curr_loc < car.dest_to_stop)
				assert car.curr_loc < (stop - 1) * LOT_LENGTH + SMALL_INTERVAL
				stopped = True
		
			elif car.status == 1 and car.stop < stop - dgap:
				assert (car.curr_loc <= car.dest_to_stop <= (stop - 1) * LOT_LENGTH )
				stopped = True

			else:
				assert (car.status == 6) or (car.status == 1 and car.stop >= stop - dgap and car.j != curr_vehicle.j)
				assert car.dest_to_stop >= (stop + dgap) * LOT_LENGTH + CAR_LENGTH or car.stop <= stop + dgap
				car_time = self.curr + ((stop + dgap) * LOT_LENGTH - CAR_LENGTH - car.curr_loc) / rateDRIV
				if stopped and car.prev.status == 1 and car.prev.j == curr_vehicle.j:
					if car.curr_loc >= (stop + dgap) * LOT_LENGTH - CAR_LENGTH + SMALL_INTERVAL: 
						import pdb; pdb.set_trace()
					break
				elif stopped:
					assert car.prev.end_time is not None
					assert car.prev.stop < stop - dgap
					if not (car.prev.stop + dgap) * LOT_LENGTH - CAR_LENGTH >= car.curr_loc - 1e-04:
						import pdb; pdb.set_trace()
					car_time = (stop - car.prev.stop) * LOT_LENGTH / rateDRIV
					if car.prev.status == 1:		
						car_time += max(meanDRIV * self.m_in[car.prev.type][car.prev.j-1] + car.prev.end_time, self.curr + max(0.0, ((car.prev.stop + dgap) * LOT_LENGTH - CAR_LENGTH - car.curr_loc)) / rateDRIV)
					else:
						car_time += max(car.prev.plin_end, self.curr + max(0.0, ((car.prev.stop + dgap) * LOT_LENGTH - CAR_LENGTH - car.curr_loc)) / rateDRIV)

				if car_time < est_pout_start + meanDRIV * self.m_out[curr_vehicle.type][curr_vehicle.j-1] - 100 * SMALL_INTERVAL:
					delayed = True
					if (car.status == 1 and car.stop <= stop + dgap) or (car.status == 6 and car.stop < stop + dgap):
						assert car.j != curr_vehicle.j
						if car.end_time < self.curr:
							import pdb; pdb.set()
						elif car.end_time == self.curr:
							req_time = max( req_time, car.end_time + SMALL_INTERVAL )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car.end_time + SMALL_INTERVAL - self.curr)
						else:
							req_time = max( req_time, car.end_time )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car.end_time - self.curr)
					else:
						if not car.dest_to_stop >= (stop + dgap) * LOT_LENGTH + CAR_LENGTH:
							import pdb; pdb.set_trace()
						car_time = car.calc_time( (stop + dgap) * LOT_LENGTH + CAR_LENGTH )
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
	'''
	def check_enter_90(self, enter_time, stop):

		assert spmatch and control == 'full'
		car = self.head
		while car is not None:
			if (spot2blk(stop) < car.block_idx - 1) or (car.status == 6):
				pass
			elif car.status == 2 and car.block_idx == 1:
				assert car.plin_start <= self.curr
				enter_time = max(enter_time, car.plin_end + meanDRIV)
			elif car.status == 2:
				assert car.plin_start <= self.curr
				if (car.block_idx - 1) * CAR_LENGTH / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
					enter_time = max(enter_time, car.plin_start + meanPLIN - (car.block_idx - 1) * CAR_LENGTH / rateDRIV)
			elif car.status == 5 and car.block_idx == 1:
				assert car.pout_start <= self.curr <= car.pout_end
				enter_time = max(enter_time, car.pout_end + meanDRIV)
			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				if ((car.block_idx - 1) * CAR_LENGTH - CAR_LENGTH) / rateDRIV < car.pout_end - self.curr:
					enter_time = max(enter_time, car.pout_end - ((car.block_idx - 1) * CAR_LENGTH - CAR_LENGTH) / rateDRIV)
			else:
				assert car.status == 1
				if enter_time < car.end_time + meanPLIN - (car.block_idx - 1) * CAR_LENGTH / rateDRIV - 10 * SMALL_INTERVAL:
					enter_time = car.end_time + meanPLIN - (car.block_idx - 1) * CAR_LENGTH / rateDRIV
			car = car.nex
		return enter_time
	'''