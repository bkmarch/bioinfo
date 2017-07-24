"""
Biocircuits and transcription networks.
"""

class gene(object):	#node representation of a gene.
	"""
	gene is a node with input and output functions.
	(inp) is list of afferent edges. transcription factors, promoters, etc.
	(out) is list of efferent edges. protein, mRNA, etc.
	"""
	def __init__(self, inp, out, seq):
		self.inp = inp
		self.out = out
		self.seq = seq	#sequence of gene
		
	def add_inp(self, inp_edge):	#inp_edge is an object. tc_f or tc_oo
		self.inp.append(inp_edge)
		
	def add_out(self, out_edge):	#same as inp_edge
		self.out.append(out_edge)
		
class tc_f(object):	#transcription smooth function
	
	"""
	X -> Y
	input is X* which is concentration of active transcription factor.
	Y is molecules of protein
	"""
	
	def __init__(self, inp):
		self.inp = inp
	
	"""
	Hill function relates concentration of TF to protein levels.
	
	(K) is activation coefficient. relates affinity of X to promoter site
		-defines conc of X needed to activate expression
		-units of concentration
	(beta) is promoter level maximum.
		-as X* >> K , f(X*) = beta
	(n) = Hill coefficient. steepness of function
		-values range from 1-4
	(beta0) is basal expression level if needed
		
	"""
	def hill_activator(self, k, beta, n, beta0 = None):		
		f_act = (beta * (self.inp ** n)) / ((k**n) + (self.inp **n)) + beta0

	def hill_repressor(self, k, beta, n, beta0 = None):
		f_rep = beta / (1 + ((self.inp/ k )**n)) + beta0

class tc_oo(object):	#transcription off/on
	
	"""
	Gene originally set to 0 expression. If we reach minimum level of k,
	turn on to maximum level beta.
	(inp) is X*
	"""
	
	def __init__(self, inp):
		self.inp = inp
		self.state = 0
	
	def hill_b(self, k, beta):
		if self.inp >= k:
			self.state = beta
					
def tc_and(list):	#Transcription AND. All genes on simultaneously.
	"""
	use the off/on for simplicity purposes
	list contains objects that have states.
	if state != 0, then gene is on (1)
	if state == 0, then gene is off (0)
	"""
	on_off = []
	gene_state = None
	for i in range(len(list)):
		if list[i].state != 0
			on_off.append(1)
		else:
			on_off.append(0)
			gene_state = False
	if gene_state = False:
		return gene_state
	else:
		gene_state = True
		return gene_state
	
	
