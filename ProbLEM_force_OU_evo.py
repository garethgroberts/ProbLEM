import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition
import os

# gareth roberts, 12 may 2023, caltech & imperial college london

print("=============================================================================================")
N = 1000 		#number of blocks (elements in spatial array) along river
dx = 1			#width of each block
dt = 1			#time step length

# critical value(s) above which erosion occurs
cmean = 2					#mean
cstd = 0					#standard deviation 
cmean = np.full(N,cmean)	#define along entire profile
cstd = np.full(N,cstd)		#define along entire profile
cvar = cstd**2				#variance

x = np.arange(0,N,dx)		#set up spatial grid
z = np.arange(N,0,-1)	    #set up initial elevation without dz

# assume erosional force is normally distributed
# forcemean=x/100						#mean
# forcestd=1							#standard deviation
# forcestd = np.full(N,forcestd)		#define along entire profile
# forcevar = forcestd**2				#variance
# force = np.random.normal(forcemean,forcestd, size=(N,))	#force along entire profile

#assume force is determined by Ornstein-Uhlenbeck process
theta=1e-3	#drift/reversion coefficient
sigma=2e-2	#volatility
mu=5		#long term mean
xsize = np.size(x) 
dw = np.random.normal(0, scale = np.sqrt(dx), size = xsize)
force = np.zeros(xsize)
initialforce= 1 #initial value for force (at river head)
force[0] = initialforce 		

for i in range(1, xsize):
	force[i] = force[i-1] + theta*(mu - force[i-1])*dx + sigma*dw[i] #Euler-M scheme: OU process

var = ((sigma**2) /(2*theta)) * (1-np.exp(-2*theta*x))	#Gardiner04 Eqn 4.4.29, with zero variance for x(0)
exp = mu + (force[0] - mu)*np.exp(-theta * x)
prob = 1 - norm.cdf((cmean-exp)/np.sqrt(var))	
probflip=np.flip(prob)
rundispa = np.cumsum(probflip) * dx					#cumulative expected displacement
runvar = (dx**2) * np.cumsum(probflip*(1-probflip))	#cumulative variance

# plt.plot(x, exp, color='grey', linewidth=2)
# plt.plot(x, exp+np.sqrt(var), color='grey', linestyle='dashed')
# plt.plot(x, exp-np.sqrt(var), color='grey', linestyle='dashed')
# plt.plot(x,force,color='black', linewidth=2)
# plt.plot(x,cmean, color='red', linewidth=2)
# plt.plot(x,prob, color='green', zorder=10)
# plt.plot(x,rundispa, color='blue', zorder=10)
# plt.xlabel("Distance, [L]")
# plt.ylabel("Erosional force, [M][L]/[T]$^2$")
# plt.show()
# exit()

#print("Expected displacement is", format(dispa,'.2f'), "[L], with variance", format(var,'.2f'), "[L], in", dt * N, "[T]")


# the basic control sequence for this model is as follows. note that indexing is from 
# 0 at the mouth to N at the head. 
# dz = z[:-1] - z[1:]			#calculate local, point-to-point, differences in relief
# dz = np.pad(dz, (0, 1))		#pad most downstream dz with 0
# dz[dz < c] = 0				#set all points with dz<c to 0 (i.e. don't erode)

#plot starting condition 
# fig, ax1 = plt.subplots()
# ax1.bar(x,z,color='gray')
# ax1.set_xlabel('Distance, [L]')
# ax1.set_ylabel('Elevation, [L]')
# ax1.set_ylim(ymin=-50,ymax=1100)
# ax1.text(0, 1025, "Time step 0")
# ax1.plot(x, 100*(exp), color='red', linewidth=2)
# ax1.plot(x, 100*(exp+np.sqrt(var)), color='red', linestyle='dashed')
# ax1.plot(x, 100*(exp-np.sqrt(var)), color='red', linestyle='dashed')
# ax1.plot(x, 100*force, color='pink', label='F, [M][L]/[T]$^2$', linewidth=1)
# ax1.plot(x, 100*cmean, color='orange', label='c, [L]', linewidth=2)
# ax1.plot(x, 100*prob, color='black', label='P(F>c)', linewidth=3)
# ax1.errorbar(N - rundispa[0], 0, xerr=runvar[0], color='black')		#plot variance
# ax1.scatter(N - rundispa[0], 0, color='red', edgecolors='black', zorder = 10)	#plot expected values]
# ax1.legend()
# 
# ax2 = plt.axes([0,0,1,1])	#plot inset histogram of relief
# ip = InsetPosition(ax1, [0.4,0.73,0.25,0.25])
# ax2.set_axes_locator(ip)
# ax2.set_xlabel('$\Delta z$, [L]')
# ax2.set_ylabel('PDF')
# ax2.set_xlim(xmin=0,xmax=100)
# ax2.set_ylim(ymin=0,ymax=0.1)
# binwidth=1
# dz = z[:-1] - z[1:] #actual dz
# ax2.hist(dz, bins=np.arange(min(dz), max(dz) + binwidth, binwidth), color='gray', density=True)
# plt.savefig('./riv_evo_OU/fig0.png', dpi=300)
# plt.show()
# plt.close()
# exit()


fig=plt.figure(figsize=(6,8), facecolor='white')

ax = fig.add_subplot(211)
plt.text(25, 1000, "Time step 0", size=14, horizontalalignment='left')
plt.bar(x,z,color='gray',width=1)
plt.xlabel('Distance, [L]')
plt.ylabel('Elevation, [L]')
plt.ylim(ymin=-50,ymax=1100)
plt.plot(x, 400*prob, color='black', label='P(F>c)', linewidth=3)
plt.text(1015, 400, "- p = 1", size=8, horizontalalignment='left')
plt.text(1015, 0, "- p = 0", size=8, horizontalalignment='left')
plt.xlim([-15, 1015])
plt.errorbar(N - rundispa[0], 0, xerr=runvar[0], color='black')		#plot variance
plt.scatter(N - rundispa[0], 0, color='red', edgecolors='black', zorder = 10)	#plot expected values]

axs_ins = inset_axes(ax, 
                    width="25%", 	# width relative to width of parent_bbox
                    height=1, 		# height in inches
                    loc=1)
plt.xlabel('$\Delta z$, [L]')
plt.ylabel('PDF')
plt.xlim(xmin=0,xmax=99)
plt.ylim(ymin=0,ymax=0.1)
dz = z[:-1] - z[1:] #actual dz
binwidth=1
plt.hist(dz, bins=np.arange(min(dz), max(dz) + binwidth, binwidth), color='gray', density=True)

ax = fig.add_subplot(212)
plt.ylim([-0.2,5])
plt.plot(x, cmean, color='orange', label='c, [L]', linewidth=2)
plt.plot(x, exp, color='grey', linewidth=3, label='$<F_x>$')
plt.plot(x, exp+np.sqrt(var), color='grey', linestyle='dashed')
plt.plot(x, exp-np.sqrt(var), color='grey', linestyle='dashed')
plt.plot(x, force, color='black', label='F, [M][L]/[T]$^2$', linewidth=1)
plt.xlabel('Distance, [L]')
plt.ylabel('Force, [L][M]/[T]$^2$')
#plt.legend()

plt.savefig('./riv_evo_OU/fig0.png', dpi=300)
#plt.show()
plt.close()
#exit()


runs = 1000			#number of time steps that the model will take

j=50
for i in range(runs):
	dw = np.random.normal(0, scale = np.sqrt(dx), size = xsize)
	force = np.zeros(xsize)
	force[0] = initialforce #initial value for force (at river head)		

	for l in range(1, xsize):
		force[l] = force[l-1] + theta*(mu - force[l-1])*dx + sigma*dw[l] #Euler-M scheme: OU process

	for k in range(N-1):
		if force[k] > cmean[k]:
			z[k] = z[k+1]

	if i == j :	
		dz = z[:-1] - z[1:] #actual dz
# 		fig, ax1 = plt.subplots()
# 		
# 		#plot evolving river profile and variance
# 		ax1.bar(x,z, color='gray')
# 		ax1.set_xlabel('Distance, [L]')
# 		ax1.set_ylabel('Elevation, [L]')
# 		ax1.set_ylim(ymin=-50,ymax=1100)
# 		ax1.text(0, 1025, "Time step {}".format(i))
# 		ax1.plot(x, 100*exp, color='red', linewidth=2)
# 		ax1.plot(x, 100*(exp+np.sqrt(var)), color='red', linestyle='dashed')
# 		ax1.plot(x, 100*(exp-np.sqrt(var)), color='red', linestyle='dashed')
# 		ax1.plot(x, 100*force, color='pink', label='F, [M][L]/[T]$^2$', linewidth=1)
# 		ax1.plot(x, 100*cmean, color='orange', label='c, [L]', linewidth=2)
# 		ax1.plot(x, 100*prob, color='black', label='P(F>c)', linewidth=3)
# 		ax1.errorbar(N - rundispa[i], 0, xerr=runvar[i], color='black')		#plot variance
# 		ax1.scatter(N - rundispa[i], 0, color='red', edgecolors='black', zorder = 10)	#plot expected values]
# 		ax1.legend()
# 
# 		#plot histogram of relief
# 		ax2 = plt.axes([0,0,1,1])
# 		ip = InsetPosition(ax1, [0.4,0.73,0.25,0.25])
# 		ax2.set_axes_locator(ip)
# 		ax2.set_xlabel('$\Delta z$, [L]')
# 		ax2.set_ylabel('PDF')
# 		ax2.set_xlim(xmin=0,xmax=100)
# 		ax2.set_ylim(ymin=0,ymax=0.1)
# 		binwidth=1
# 		ax2.hist(dz,color='gray', bins=np.arange(min(dz), max(dz) + binwidth, binwidth), density=True)
		
		
		fig=plt.figure(figsize=(6,8), facecolor='white')

		ax = fig.add_subplot(211)
		plt.text(25, 1000, "Time step {}".format(i), size=14, horizontalalignment='left')
		plt.bar(x,z,color='gray',width=1)
		plt.xlabel('Distance, [L]')
		plt.ylabel('Elevation, [L]')
		plt.plot(x, 400*prob, color='black', label='P(F>c)', linewidth=3)
		plt.text(1015, 400, "- p = 1", size=8, horizontalalignment='left')
		plt.text(1015, 0, "- p = 0", size=8, horizontalalignment='left')
		plt.xlim([-15, 1015])
		plt.ylim(ymin=-50,ymax=1100)
		plt.errorbar(N - rundispa[i], 0, xerr=runvar[i], color='black')		#plot variance
		plt.scatter(N - rundispa[i], 0, color='red', edgecolors='black', zorder = 10)	#plot expected values]

		
		axs_ins = inset_axes(ax, 
                    	width="25%", 	# width relative to width of parent_bbox
                    	height=1, 		# height in inches
                    	loc=1)
		plt.xlabel('$\Delta z$, [L]')
		plt.ylabel('PDF')
		plt.xlim(xmin=0,xmax=99)
		plt.ylim(ymin=0,ymax=0.1)
		dz = z[:-1] - z[1:] #actual dz
		binwidth=1
		plt.hist(dz, bins=np.arange(min(dz), max(dz) + binwidth, binwidth), color='gray', density=True)

		ax = fig.add_subplot(212)
		plt.ylim([-0.2,5])
		plt.plot(x, cmean, color='orange', label='c, [L]', linewidth=2)
		plt.plot(x, exp, color='grey', linewidth=3, label='$<F_x>$')
		plt.plot(x, exp+np.sqrt(var), color='grey', linestyle='dashed')
		plt.plot(x, exp-np.sqrt(var), color='grey', linestyle='dashed')
		plt.plot(x, force, color='black', label='F, [M][L]/[T]$^2$', linewidth=1)
		plt.xlabel('Distance, [L]')
		plt.ylabel('Force, [L][M]/[T]$^2$')
#		plt.legend()


#		plt.show()
		plt.savefig('./riv_evo_OU/fig%s.png' % i, dpi=300)
		plt.close()
		print('fig%s.png' % i)


		j+=50
		
exit()

			
			
	