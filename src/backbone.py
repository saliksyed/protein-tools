import pyrr
from pyrr.matrix44 import Matrix44
import numpy as np

class Joint:
	def __init__(self):
		self.children = []
		self.state = {}
		self.position = np.array([0,0,0])
		self.relative_position = np.array([0,0,0])

	def get_state(self):
		return self.state

	def type(self):
		return self.__class__.__name__

	def set_name(self, name):
		self.name = name
		self.state[name] = self.position

	def set_relative_position(self, pos):
		# update children base positions:
		self.relative_position = pos

	def set_position(self, pos):
		if self.name:
			self.state[self.name] = pos
		self.position = pos

	def current_position(self):
		return self.position

	def add_child(self, joint):
		pos = self.get_next_child_position(joint)
		joint.set_relative_position(pos)
		self.children.append(joint)
	
	def get_force(self, joint):
		raise 'Abstract'

	def get_next_child_position(self, joint):
		raise 'Abstract'

	def update(self, t, dt, state):
		# extend this class for a new base
		raise 'Abstract'

	def edges(self):
		edges = {}
		edges[self.name] = map(lambda x : x.name, self.children)
		for child in self.children:
			edges.update(child.edges())
		return edges

	def get_axis(self):
		raise 'Abstract'

	# private methods
	def _tick(self, t, dt, state, root_transform=None):
		if root_transform == None:
			root_transform = Matrix44.identity()
		# apply the root transform to the current frame:
		v4pos = np.array([self.relative_position[0], self.relative_position[1], self.relative_position[2], 1.0])
		new_pos = pyrr.matrix44.apply_to_vector(root_transform, v4pos)
		self.set_position(np.array([new_pos[0], new_pos[1], new_pos[2]]))
		translate = pyrr.matrix44.create_from_translation(self.relative_position)
		next_transform = pyrr.matrix44.multiply(translate, root_transform)
		for child in self.children:
			self.state.update(child.tick(t, dt, state, next_transform))
		return self.state

class FixedJoint(Joint):
	def __init__(self, position):
		Joint.__init__(position)

	def update(self, t, dt, state):
		pass
	
	def get_axis(self):
		return None

	def tick(self, t, dt, state, transform=None):
		return self._tick(t, dt, state)

class RevoluteJoint(Joint):
	def __init__(self):
		Joint.__init__(self)

	def tick(self, t, dt, state, transform=None):
		# compute the transform at time t
		r = self.update(t, dt, state)
		mat = pyrr.matrix44.create_from_axis_rotation(self.get_axis(), r)
		if transform != None:
			mat = pyrr.matrix44.multiply(mat, transform)
		return self._tick(t, dt, state, mat)

class TranslationalJoint(Joint):
	def __init__(self):
		Joint.__init__(self)

	def tick(self, t, dt, state, transform=None):
		# compute the transform at time t
		d = self.update(t, dt, state)
		mat = pyrr.matrix44.create_from_translation(d * self.get_axis())
		if transform != None:
			mat = pyrr.matrix44.multiply(mat, transform)
		return self._tick(t, dt, state, mat)

class Simulator:
	def __init__(self, symbols):
		# symbols is a mapping from symbol => Class
		self.symbols = symbols
		self.symbol_counts = {}
		self.t = 0

	def _construct_from_symbol(self, symbol):
		if not symbol in self.symbol_counts:
			self.symbol_counts[symbol] = 0
		self.symbol_counts[symbol] += 1
		name = "%s_%d" % (symbol, self.symbol_counts[symbol])
		joint = self.symbols[symbol]()
		joint.set_name(name)
		return joint

	def parse(self, code, recursive=False):
		if isinstance(code, list):
			root = self._construct_from_symbol(code[0])
			children = code[1:]
			for child in children:
				parsed_child = self.parse(child, True)
				root.add_child(parsed_child)
		else:
			root = self._construct_from_symbol(code)

		if not recursive:
			self.root = root
		return root

	def tick(self, dt):
		self.t += dt
		state = self.root.get_state()
		self.root.tick(self.t, dt, state)
		ret = {}
		ret['nodes'] = self.root.get_state()
		ret['edges'] = self.root.edges()
		return ret



