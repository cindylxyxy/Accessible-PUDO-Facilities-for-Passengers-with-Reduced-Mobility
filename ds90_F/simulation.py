########################################################################################################################################################
# Last updated: Aug 1, 2022 by Xinyu Liu
# This file initiates and runs the trajectory-based simulation for Double-Sided 90-degree Configuration and the specified operational policies. 
########################################################################################################################################################

import sys
import csv
import pickle
from copy import deepcopy
import numpy as np

from params import *
from event import *
from event90 import *
from utils import *

''' calculate the unit area for each layout configuration '''
assert angle == 90
assert side == 'double'
A = LOT_LENGTH * (2 * LOT_WIDTH + LANE_WIDTH)

''' print the simulation input '''
print ('printing input params for %s sided layout ...'%side)
print ('angle:', angle)
print ('type of simulation:', simType)
print ('loading spot configuration:', 'Lot length is %s, lot width is %s, and lane width is %s'%(LOT_LENGTH, LOT_WIDTH, LANE_WIDTH) )
print ('number of hours in simulation:', SIM_HOUR)
print ('number of simulation iterations:', SIM_ITER) 
print ('average service time (sec):', meanSERV)
print ('average driving speed (ft/sec):', avgSPEED)

sys.stdout.flush()

''' a function that returns the maximal throughput and the corresponding policy given N and N_s '''
def eval_Ns(N, N_s):
	
	print ('Starting N_s =', N_s)

	''' calculate the facility area w.r.t. N and N_s '''
	assert angle == 90
	area = (ceil(N_s / 2) + N) * A

	thrupt_vec = []
	std_vec = []
	cutoff_vec = []
	
	''' iterate over cutoff values '''
	for k in range(0, 1 + N_s):	

		if N == N_s == k:
			print ('Optimal cutoff value for N = %s and N_s = %s is %s.'%(N, N_s, cutoff_vec[-1]))
			return (cutoff_vec[-1], thrupt_vec[cutoff_vec[-1]], std_vec[cutoff_vec[-1]], area)

		print ('Cutoff value:', k)
		sys.stdout.flush()

		try:
			''' check if the requested data has been generated before '''
			visited = pickle.load(open(dirname + '%s_%s_%s_%s.p'%(filename, N, N_s, k), 'rb'))		
		except FileNotFoundError:
			visited = None

		if visited is not None:
			thrupt_vec.append(visited['thrupt'])
			std_vec.append(visited['std_thrupt'])
			cutoff_vec.append(np.argmax(thrupt_vec))
			# print (thrupt_vec, cutoff_vec)

			if len(cutoff_vec) >= 3 and cutoff_vec[-1] == cutoff_vec[-2] == cutoff_vec[-3]:
				print ('Optimal cutoff value for N = %s and N_s = %s is %s.'%(N, N_s, cutoff_vec[-1]))
				return (cutoff_vec[-1], thrupt_vec[cutoff_vec[-1]], std_vec[cutoff_vec[-1]], area)
			continue

		''' initialize the simulation run '''
		test = system90(N, N_s, k, SEED_SERV[N-1], SEED_POUT[N-1], SEED_PLIN[N-1], SEED_DRIV[N-1])
		for j in range(1, test.N + 1):
			if test.spot_types[j-1] < 2:
				test.inCount += 1
				car = vehicle(test, Vtype = test.spot_types[j-1], getin = True, j = j)
				test.inservice[j-1] = car
				test.add_event( event(car.serv_time, car, 'prepare_pulling_out') )
		assert test.inservice.count(None) + test.inCount == test.N	

		''' record the # of vehicles leaving the facility '''
		count = []
		for i in range(SIM_ITER):
			count.append( test.run() )

		total_out = count[-1]
		total_in = total_out
		count = [count[0]] + [count[i+1] - count[i] for i in range(SIM_ITER-1)]

		
		''' calculate the system times for remaining vehicles '''
		car = test.head
		while car is not None:
			total_in += 1
			car.update_loc()
			if car.status == 1:
				if not car.prod_time == 0.0:
					import pdb; pdb.set_trace()
				if car.curr_loc == 0.0:
					total_in -= 1
					car = car.nex
					continue
				else:
					car.prod_time += (car.curr_loc / car.driv)

			elif car.status == 2:
				if not car.plin_end >= test.curr:
					import pdb; pdb.set_trace()
				car.prod_time -= (car.plin_end - test.curr)

			elif car.status == 3:
				assert mode == 'singlelane'
				total_in += 1
				if not car.serv_end >= test.curr:
					import pdb; pdb.set_trace()
				car.prod_time -= (car.serv_end - test.curr)

			elif car.status == 5:
				if not car.pout_end >= test.curr:
					import pdb; pdb.set_trace()
				car.prod_time -= (car.pout_end - test.curr)

			else:
				assert car.status == 6
				if spmatch:
					start_x = car.block_idx * CAR_LENGTH
				else:
					start_x = (car.stop + dgap) * LOT_LENGTH
				if not car.curr_loc >= start_x:
					import pdb; pdb.set_trace()
				car.prod_time += ((car.curr_loc - start_x) / car.driv)

			test.idle_time_calc(car)
			car = car.nex

		for j in range(1, test.N + 1):
			if test.inservice[j - 1] is not None and test.inservice[j - 1].status in [3, 4]:
				car = test.inservice[j - 1]
				total_in += 1
				if car.status == 3:
					if not car.serv_end >= test.curr:
						import pdb; pdb.set_trace()
					car.prod_time -= (car.serv_end - test.curr)
				test.idle_time_calc(car)

		total = total_in

		visited = dict()
		visited['max_idle'] = test.max_idle
		visited['max_pct_idle'] = test.max_pct_idle
		visited['tot_pct_prod'] = np.sum(test.prod_time) / total
		visited['tot_pct_idle'] = np.sum(test.idle_time) / total
		visited['avg_pct_idle'] = np.sum(test.pct_idle) / total
		visited['thrupt'] = np.mean(count) / SIM_HOUR
		visited['std_thrupt'] = np.std(count) / sqrt(SIM_ITER) / SIM_HOUR

		pickle.dump(visited, open(dirname + '%s_%s_%s_%s.p'%(filename, N, N_s, k), 'wb'))

		thrupt_vec.append(visited['thrupt'])
		std_vec.append(visited['std_thrupt'])
		cutoff_vec.append(np.argmax(thrupt_vec))
		# print (thrupt_vec, cutoff_vec)

		if len(cutoff_vec) >= 3 and cutoff_vec[-1] == cutoff_vec[-2] == cutoff_vec[-3]:
			print ('Optimal cutoff value for N = %s and N_s = %s is %s.'%(N, N_s, cutoff_vec[-1]))
			return (cutoff_vec[-1], thrupt_vec[cutoff_vec[-1]], std_vec[cutoff_vec[-1]], area)

	print ('Optimal cutoff value for N = %s and N_s = %s is %s.'%(N, N_s, cutoff_vec[-1]))
	return (cutoff_vec[-1], thrupt_vec[cutoff_vec[-1]], std_vec[cutoff_vec[-1]], area)


if (not debug):
	outfile = open(dirname + 'res_traj_%s.csv'%filename, 'w')
	writer = csv.writer(outfile, delimiter = ',')
	writer.writerow( ['half_N', 'N', 'N_s', 'area', 'cutoff', 'thrupt', 'std thrupt'] )

opt_thrupt = dict()

''' main loop - iterate over N '''
for N in range(1, 1 + LOT_TOTAL):

	print ('\nStarting N =', N)

	if N not in opt_thrupt:
		opt_thrupt[N] = dict()

	curr_max = 1

	''' iterate over N_s '''
	for N_s in range(1, 1 + N):
		assert N_s not in opt_thrupt[N]
		opt_thrupt[N][N_s] = eval_Ns(N, N_s)

		if opt_thrupt[N][N_s][1] > opt_thrupt[N][curr_max][1]:
			curr_max = N_s

	if (not debug):
		writer.writerow([N, 2 * N, N_s, opt_thrupt[N][N_s][3], opt_thrupt[N][N_s][0], opt_thrupt[N][N_s][1], opt_thrupt[N][N_s][2]])			
