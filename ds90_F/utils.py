from math import ceil, floor
from collections import *
from params import *
import random 

changePoint = namedtuple( 'changePoint', ['x', 't', 'v'])

class time():

	def __init__(self, val):
		self.val = round(val, 3)

	def __str__(self):
		return str(self.val)

	def __eq__(self, other):
		return self.val == other.val

	def __gt__(self, other):
		return self.val > other.val

	def __lt__(self, other):
		return self.val < other.val

class ParamGen():

	def __init__(self, gen):
		self.gen = gen
		self.count = 0

	def next(self):
		self.count += 1
		val = - 0.1
		while val <= 0.0:
			# val = round( next(self.gen), 2 )
			val = next(self.gen)
		return val
		# return max(floor( val / meanDRIV) * meanDRIV, SMALL_INTERVAL)

def Expo(mean, seed = None):
	# returning time
	curr = random.Random(seed)
	while True:
		yield curr.expovariate(1./mean)

def Erla(mean, m, seed = None):
	# returning time
	curr = random.Random(seed)
	while True:
		yield curr.gammavariate(alpha = m, beta = mean)

def Unif(mean, seed = None):
	# returning time
	curr = random.Random(seed)
	while True:
		yield curr.uniform(0, 2. * mean)

def Unif2(mean, seed = None):
	# returning rate
	curr = random.Random(seed)
	while True:
		yield curr.uniform(0.5 * mean, 1.5 * mean)

def Norm(mean, sigma = None, seed = None):
	# returning time
	curr = random.Random(seed)
	if sigma == None:
		sigma = mean
	while True:
		yield curr.gauss(mean, sigma)	

def Tri1(mode, low = None, high = None, seed = None):
	# returning time
	curr = random.Random(seed)
	if low == None:
		low = mode / 3.
	if high == None:
		high = 5 * mode / 3.
	while True:
		yield curr.triangular(low, high, mode)

def Tri0(mode, low = None, high = None, seed = None):
	# returning time
	curr = random.Random(seed)
	if low == None:
		low = .5 * mode
	if high == None:
		high = 1.5 * mode
	while True:
		yield curr.triangular(low, high, mode)

def Tria(mode, low = None, high = None, seed = None):
	# returning time
	curr = random.Random(seed)
	if low == None:
		low = 0
	if high == None:
		high = 2. * mode
	while True:
		yield curr.triangular(low, high, mode)

def Cons(mean, seed = None):
	# returning mean
	while True:
		yield mean

def Disc(mean, seed = None):
	# returning mean
	curr = random.Random(seed)
	while True:
		yield curr.choice([ 0.9 * mean, 0.95 * mean, mean, 1.05 * mean, 1.1 * mean])


class event():

	def __init__(self, time, vehicle, typ):
		if time == None:
			import pdb; pdb.set_trace()
		self.time = time
		self.vehicle = vehicle
		self.typ = typ

	def __eq__(self, other):
		return (float(self.time) == float(other.time)) & (self.typ == other.typ) & (self.vehicle.idx == other.vehicle.idx) 

	def __gt__(self, other):
		if float(self.time) == float(other.time):
			if (self.typ == other.typ):
				if (self.vehicle.idx is None) | (other.vehicle.idx is None):
					assert (self.vehicle.idx is None) and (other.vehicle.idx is None)
					assert (self.typ == 'enter_system')
					return False
				assert self.vehicle.j is not None and other.vehicle.j is not None
				if (self.vehicle.j == other.vehicle.j):
					if not self.typ in ['start_pulling_in', 'start_second_enter', 'leave_system']:
						import pdb;pdb.set_trace()
					return (self.vehicle.idx > other.vehicle.idx)
				return (self.vehicle.j < other.vehicle.j)
			return (event_priority[self.typ] > event_priority[other.typ])
		return (float(self.time) > float(other.time))

	def __lt__(self, other):
		if float(self.time) == float(other.time):
			if self.typ == other.typ:
				if (self.vehicle.idx is None) | (other.vehicle.idx is None):
					assert (self.vehicle.idx is None) and (other.vehicle.idx is None)
					assert (self.typ == 'enter_system')
					return False
				assert self.vehicle.j is not None and other.vehicle.j is not None
				if (self.vehicle.j == other.vehicle.j):
					if not self.typ in ['start_pulling_in', 'start_second_enter', 'prepare_pulling_out','leave_system']:
						import pdb;pdb.set_trace()
					return (self.vehicle.idx < other.vehicle.idx)
				return (self.vehicle.j > other.vehicle.j)
			return (event_priority[self.typ] < event_priority[other.typ])
		return (float(self.time) < float(other.time))

	def __ne__(self, other):
		return not self == other

	def __ge__(self, other):
		return not self < other

	def __le__(self, other):
		return not self > other 
		


def sortHelper(headNode, tlen1, tlen2):

	curr = headNode
	print ('starting...')
	while curr.nex:
		print (curr.data)
		curr = curr.nex
	print (curr.data)
	print ('target length is:', tlen1, tlen2 )

	if headNode != None and headNode.nex != None:
		prev = headNode
		curr = headNode.nex
		assert tlen1 >= 1
		assert tlen2 >= 1
		assert tlen1 >= tlen2
		assert tlen1 - tlen2 <= 1
		for _ in range(tlen1 - 1):
			prev = curr
			curr = curr.nex
		prev.nex = None
			
		curr1 = sortHelper(headNode, tlen1 - tlen1//2, tlen1//2)
		curr2 = sortHelper(curr, tlen2 - tlen2//2, tlen2//2)

		if curr1.data <= curr2.data:
			headNode = curr1
			curr1 = curr1.nex
		else:
			headNode = curr2
			curr2 = curr2.nex
		curr = headNode 

		while curr1 is not None and curr2 is not None:
			if curr1.data <= curr2.data:
				curr.nex = curr1
				curr1 = curr1.nex
			else:
				curr.nex = curr2
				curr2 = curr2.nex
			curr = curr.nex

		if curr1 is not None:
			assert curr2 is None
			curr.nex = curr1
		if curr2 is not None:
			assert curr1 is None
			curr.nex = curr2
	return headNode



class Node:

	def __init__(self, data = None):
		self.data = data
		self.prev = None
		self.nex = None



class DLinkedList:

	def __init__(self):
		self.head = None
		self.end = None
		self.length = 0
		
	def __len__(self):
		return self.length

	def print(self):
		curr = self.head
		while (curr):
			print (curr.data)
			curr = curr.nex

	def isEmpty(self):
		assert (self.head is None) == (self.length == 0)
		return (self.head is None)

	def addHead(self, new):
		NewNode = Node(new)
		NewNode.nex = self.head
		self.head = NewNode
		self.end = NewNode
		self.length += 1

	def addEnd(self, new):
		NewNode = Node(new)
		if self.head is None:
			self.head = NewNode
			return
		curr = self.head
		while (curr.nex):
			curr = curr.nex
		curr.nex = NewNode
		NewNode.prev = curr
		self.end = NewNode
		self.length += 1

	def insert(self, loc, new):
		if loc == 0:
			self.addHead(new)
		elif loc >= self.length:
			self.addEnd(new)
		else:
			NewNode = Node(new)
			prev = self.head
			for _ in range(loc - 1):
				prev = prev.nex
			curr = prev.nex
			prev.nex = NewNode
			prev.nex.prev = prev
			prev.nex.nex = curr
			curr.prev = prev.nex
			self.length += 1

	# def sort(self):
	# 	self.head = sortHelper(self.head, self.length - self.length//2, self.length//2)

	def getVal(self, idx):
		curr = self.head
		while (nstep == idx > 0):
			curr = curr.nex
		return curr.data

	def __contains__(self, data):
		curr = self.head
		while curr != None:
			if curr.data == data:
				return True
			else:
				return False
