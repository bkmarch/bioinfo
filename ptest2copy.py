from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


#work on the differential equation portion first.
#https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.integrate.odeint.html#scipy.integrate.odeint
#DE methods
#http://tutorial.math.lamar.edu/Classes/DE/IntroFirstOrder.aspx

"""
Network Motifs- Given a set of nodes and edges we find network motifs by 
comparing it to a randomized graph with same number of nodes and edges. 
Significant patterns denoted as 'network motifs'. 

Negative Autoregulation - self edges. type of network motif
Feed forward motif. X-> Y, X->Z, Y->Z. all arrows signify promoter activity.

'Erdos and Renyi' randomized network. 

"""

min_level = {'x':[],
			'y':[['x',2, None]],
			'z':[['x',2, None],['y',1, None]]
			}	#thresholds that must be hit to activate certain proteins
				# None is a place holder for the time at which a parent reaches threshold.
beta = {'x':5,
		'y':3,
		'z':3
		}
alpha = {'x':2,
		'y':4,
		'z':5
		}				
				
class p_time():	#calculate protein levels at various times.

	"""
	time model of protein levels
	beta is production rate of protein
	alpha is degradation rate
	y is protein amount
	"""

	def __init__(self , beta, alpha, y0, t, minlvl):	
		self.beta = beta
		self.alpha = alpha
		self.y0 = y0
		self.t = t
			
	def p_level(self):	#protein level solver
		y00 = np.array(self.y0)	#initial concentration		
		x_out = odeint(self.ode, y00, self.t, args=(self.beta,self.alpha))	#ODE solution
		return x_out
	
	def ode(self, y, t, beta, alpha):	#ODE equation
		dydt = np.array(self.beta - self.alpha * y)
		return dydt
		
	def p_t_level(self, lvl):	#returns the Time at which protein reaches specified lvl.
		#used first order linear ODE. integrating factor method. floats to keep decimals
		time = np.log(float(float(-1 * float(lvl) * self.alpha / self.beta) + 1)) / (-1 * self.alpha)
		return time
	
	
def threshold(protein, dic_min):	#check at what times predecessors hit minimum activation lv	
	for i in dic_min:
		if bool(dic_min[i]) == False:	#empty set
			p_curve_calculate(0,100)		
		else:
			pass#calculate p level curve starting at last threshold time.

								
def p_curve_calculate(start,datapts):						
	
	t = np.linspace(start, start+3 ,datapts)	#time interval. saying stop taking data after 3seconds.
	k = p_time(5, 2, [0], t, min_level)	#(beta, alpha, initial y0 array, time interval, threshold dictionary)	
	i = k.p_level()	#calculation of protein level
	plt.plot(t, i[:, 0], 'y', label='protein level')
	plt.legend(loc='best')
	plt.xlabel('t')
	plt.grid()
	plt.show()


	
threshold('x' , min_level)
threshold('y' , min_level)



#the start function needs to pair p_time class with s_activate class. only when s_activate
# gives the go ahead, can we start using p_time.

"""
modelling a feed forward loop. x->y. x->z. y->z. x and y needed.
def ff3(x,y,z):	#feed forward 3 node model
	
	
	time = np.linspace(0, 10, 101)
	nodes = [x,y,z]	#instances of p_time
	
	
	x = 0	#initialize x to concentration 0.
"""

#x= p_time(True)
#y= p_time(False)
#z= p_time(False)

