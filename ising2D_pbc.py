#
# Code to compute average magnetization of 2D Ising magnet with no applied field
#
# USAGE:
# python ising2D_pbc.py [config file name]
#
# CONFIG FILE FORMAT:
#
# logFile = [log file name]
# temperature = [temperature in units of J/k (float)]
# latticeSize = [lattice size (int)]
# nIter = [number of Monte Carlo iterations (int)]
# deltaWrite = [frequency to write running averages (int)]
#
#
# The Onsager analytic solution in 2D gives Tc=2.3J/k
# The Mean-field model in 2D predicts Tc=4J/k
#


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import sys


cfg_file = sys.argv[1] # read the config file name

#################################################################################################
###############################            SUBROUTINES          #################################
#################################################################################################

# read the configuration file and populate the global variables
def parseConfigFile(cfg_file):
	global N, nIter, deltaItt, invT, log_file, T
	f = open(cfg_file)
	log_file = ""
	N = ""
	nIter = ""
	deltaItt = ""
	T = ""
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			# check value
			if option.lower()=='logfile':
				log_file = value
			elif option.lower()=='latticesize':
				N = int(value)
			elif option.lower()=='niter':
				nIter = int(value)
			elif option.lower()=='deltawrite':
				deltaItt = int(value)
			elif option.lower()=='temperature':
				T = float(value)
				invT = 1.0/T
			else :
				print "Option:", option, " is not recognized"
	f.close()
	# check that all variables are defined, assign defaults if need be and then print values
	if log_file:
		print "log file name:", log_file
	else:
		print "Must define log file name in config file with following line:"
		print "logFile = [log file name]"
		sys.exit(000)
	if N:
		print "Lattice size:", N
	else:
		# assign default value of 10
		N = 10
		print "Lattice size assigned default value of:", N
	if nIter:
		print "Number of Monte Carlo iterations:", nIter
	else:
		#assign default value of 100000
		nIter = 100000
		print "Number of Monte Carlo iterations assigned default value of:", nIter
	if deltaItt:
		print "Writing running averages every ", deltaItt, " iterations"
	else:
		# assign defualt
		deltaItt = 100
		print "Writing running averages every ", deltaItt, " iterations"
	if T:
		print "Temperature (in units of J/k):", T
	else:
		# assign defualt
		T = 10.0
		invT = 1.0/T
		print "Temperature assigned default value (in units of J/k):", T

# compute the total energy of a lattice configuration for the Ising model with zero applied field
def TotalEnergy(lattice):
	totalEnergy=0
	for i in range(len(lattice)):
		for j in range(len(lattice)):
			S = lattice[i,j]
			# sum spins of nearest neighbors (mod N comes from applying periodic boundary conditions)
			WF = lattice[(i+1)%N, j] + lattice[i,(j+1)%N] + lattice[(i-1)%N,j] + lattice[i,(j-1)%N]
			totalEnergy += -WF*S # Each neighbor gives energy 1.0
	return totalEnergy/2. # Each par counted twice

# generate a random configuation of square 2D lattice with binary occupation
def RandomL(N):
	lattice = np.zeros((N,N), dtype=int)
	for i in range(N):
		for j in range(N):
			lattice[i,j] = np.sign(2*np.random.rand()-1)
	return lattice


# generate a random configuation of square 2D lattice with binary occupation
def DividedStart(N):
	lattice = np.zeros((N,N), dtype=int)
	for i in range(N):
		for j in range(N):
			if j<N/2:
				lattice[i,j] = 1
			else:
				lattice[i,j] = -1
	return lattice

#################################################################################################
###############################            MAIN PROGRAM         #################################
#################################################################################################

#read parameters from cfg file
parseConfigFile(cfg_file)

# Initialize grid
#lattice = RandomL(N)
lattice = DividedStart(N)
print "Initial lattice:"
print lattice
# initial lattice energy
energy = TotalEnergy(lattice)
N2 = N*N
# Measurements
Mn = sum(sum(lattice))
print "Initial magnetism:", Mn
avgE = 0
avgMn = 0
count = 0
# open a log file
log = open(log_file,'w')
print "starting Monte Carlo loop"
# Run Monte Carlo iterations
for itt in range (0, nIter):
	# pick a random integer between 1 and NxN
	t = int(N2*np.random.rand())
	# determine x and y values of random integer
	i = t%N
	j = int(t/N)
	# determine spin of randomly selected particle
	S = lattice[i,j]
        # spins of  nearest neighbors 
	WF = lattice[(i+1)%N, j] + lattice[i,(j+1)%N] + lattice[(i-1+N)%N,j] + lattice[i,(j-1+N)%N]
	deltaE = 2*S*WF
	rand = np.random.rand()
	# check to see if we should accept the move based on the Metropolis criteria
	if deltaE<0 or math.exp(-float(deltaE)*invT)>rand: # accept move based on metropolis criteria
		# change spin
		lattice[i,j] = -S;
		# modify total energy
 	        energy += deltaE
		# update magnetism
		Mn -= 2*S;
	if itt%deltaItt==0:
		avgE += energy
		avgMn += Mn
		count += 1
		print "Step ", itt, ":\n", lattice, "\n"
		log.write("%10d %10.5f %10.5f\n" % (itt, float(avgE)/float(count), float(avgMn)/float(count)))
		log.flush()



# Finish averages
avgE = float(avgE)/float(count)
avgMn = float(avgMn)/float(count)
log.write("%10d %10.5f %10.5f\n" % (itt, avgE, avgMn))
log.flush()
log.close()
# print final lattice to standard out
print "Final lattice:"
print lattice
# plot
data = np.loadtxt(log_file)
#plt.figure(1)
#plt.subplot(211)
#plt.title('PBC Monte Carlo Ising Model for grid %d at temperature %d' %(N, T))
#plt.ylabel('Energy')
#plt.plot(data[:,0],data[:,1])
#plt.subplot(212)
#plt.ylabel('Magnetization')
#plt.xlabel('Monte Carlo Step')
#plt.plot(data[:,0],data[:,2])
#plt.show()
