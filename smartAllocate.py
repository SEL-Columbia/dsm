import subfunctions_smartAllocate as sf
import pandas as p
import numpy as np
import scipy.io as spio
import datetime as dt
import matplotlib.pyplot as plt


## load in mat weather data
mat = spio.loadmat('resourceSolTim1.mat')
data = mat['MaliNTS']

## create dates from weather data using 'list comprehension'
dates = [dt.datetime(w[0], w[1], w[2], w[3]-1) for w in data]
r , c = dates.shape
## extract resource data from weather data
resource = np.array([w[4] for w in data])	

## load up demand data
mat2 = spio.loadmat('subMetDem.mat') # This is a 8760x20 demand profile that does not yet exist
demand = mat2['subMetDem'] # This is a 8760x20 demand profile that does not yet exist
r2, c2 = demand.shape
## Subfunction resourceCalc inputs
phi_c = 0.
lats = 13.45
sigma = lats
rho = 0.2

## Loop here to create vector I_C
I_C = np.zeros(len(demand))
for i, date in enumerate(dates):
	I_C[i] = resourceCalc(date,sigma,phi_c,resource[i],lats,rho)

## Microgrid hardware parameters 
minA = 150. # W-hr. Minimum amount of energy allocated to each submeter
minAT = np.sum(minA*ones(c2)) # W-hr. Minimum allocated to the entire micro-grid
pvCap = 2000. # W. PV module peak output
batCap = 20000. # W-hr. Peak battery storage capacity
batMin = batCap * 0.50 #W-hr. Minmum Charge state of the battery.
supply = I_C/1000. * pvCap # W. Converting insolation
# on collector into generated electricity.

## Parameters for Script Functionality
HL = 336  #number of hours I want the histor of customers to be stored
batChar = zeros(r2)
supplyM = zeros((r2,c2))

## Begin time dependent energy balance and allocation simulation
for ix in np.arange(1,r2):

	# What to do before there is sufficient history.
	# Everyone gets what they want until the battery 
	# runs out

	if (HL-2) < ix:
		batChar[ix], supplyM[ix-1,:] = sf.energyBalance(batChar[ix-1],supply[ix-1],demand[ix-1,:],batCap,batMin) 
	
	# The history based algorithm will start on day with (hour number) >= (HL + 1)
	# Once the sun begins to rise
	# Indexing was to make sure that I capture midnight in the elif statement
	
	elif ((HL-2)< ix < (HL + 12)) and (supply == 0):
		# Calculate most recent 14 days of energy usage history
		# with days starting at midnight
		
		batChar[ix], supplyM[ix-1,:] = sf.energyBalance(batChar[ix-1],supply[ix-1],demand[ix-1,:],batCap,batMin) 
	
		if dates[ix-1,4] == 1:
			H = supplyM[(ix-1-336):(336-1),:] 
		
	else:
		
		# The energy supplied to each meter is equal
		# to the demand at each meter. If this is not 
		# true, this section of code will correct it.
		supplyM[ix-1,:] = demand[ix-1,:]

		# Calculate most recent 14 days of energy usage history
		# with days starting at midnight
		if dates[ix-1,4] == 1:
			H = supplyM[(ix-1-336):(336-1),:] 
		
		# If it is sunrise
		if (supply[ix-1] !=0) and (supply[ix-1-1] == 0):
			solDayStart = ix-1 # start my solar day
			jx = 1
			batChar[ix], supplyM[ix-1,:] = sf.energyBalance(batChar[ix-1],supply[ix-1],demand[ix-1,:],batCap,batMin)
			
		elif (supply[ix-1] !=0) and (supply[ix-1-1] != 0):
			jx = jx + 1
			batChar[ix], supplyM[ix-1,:] = sf.energyBalance(batChar[ix-1],supply[ix-1],demand[ix-1,:],batCap,batMin)
			
		# The sun just went down
		# Square away all accounts so that excess capacity is allotted accordingly
		elif (supply[ix-1] == 0) and (supply[ix-1-1] != 0):
		
			# Calculate energy Generated for this day
			HGD = np.sum(supply[solDayStart:(ix-1)]) # History of Generation for Day
			
			# if the daily generated capacity is greater than the
			# combined minimum allocation for all meters do this.
			if HGD >= minA*c2: 
				ERem = batChar[ix-1] - batMin #Energy Remaining
				EUsed = supplyM[solDayStart:(ix-1),:] #Energy Used
				
				# Net Energy above Minimum allotment 
				# Remaining energy - (guaranteed energy - used guarenteed energy)
				EBonus = ERem - (minAT - np.sum(EUsed))  
				
				#Individual Meter Allotments
				# guaranteed energy + percentage of bonus weight by 
				# amount of energy used by meter during H
				ELim = minA.*ones(c2) + EBonus.*H.sum(axis=0)/sum(H) 
				
			# if the daily generated electricity is less than the
			# combined minimum allocation for all meters 
			# make the Energy Limit on all meters equal
			
			if HGD <= minAT:
				ERem = batChar[ix-1] - batMin #Energy Remaining
				EUsed = supplyM[(ix-jx-1):(ix-1),:] #Energy Used
				ELim =  HGD/c2*ones(c2)
				
			# shut off any meters that have surpassed energy limit
			EUsedSum = EUsed.sum(axis=0)
			for kx in np.arange(0,c2-1):
				if EUsedSum[kx] > ELim[kx]:
					supplyM[ix-1,kx] = -999
					
			# Energy balance
			batChar[ix+1-1] = batChar[ix-1] + supply[ix-1] + np.sum(demand[ix-1,:]*(supplyM[ix-1,:] != -999))
			if batChar[ix+1-1] > batCap:
				batChar[ix+1-1] = batCap
			elif batChar[ix+1-1] < batMin:
				batChar[ix+1-1] = batMin
				# Distributing remaining electricity accordingly
				for kx in np.arange(0,c2-1):
					if supplyM[ix-1,kx] != -999:
						supplyM[ix-1,kx] = supplyM[ix-1,kx]/sum(supplyM*[ix-1,:]*(supplyM[ix-1,:]!=-999))*(batChar[ix-1]-batMin)
		
		# The sun is down and it has been down for at more than one hour
		elif (supply[ix-1] == 0) and (supply[ix-1-1] == 0): 
		
			#	If the meter was shut off during the previous hour
			#	the meter is still shut off
			for kx in np.arange(0,c2-1):
				if supplyM[ix-1-1,kx] == -999:
					supplyM[ix-1,kx] = -999
			
			batChar[ix+1-1] = batChar[ix-1] + supply[ix-1] + sum(demand[ix-1,:]*(supplyM[ix-1,:] != -999))
			if batChar[ix+1-1] > batCap:
				batChar[ix+1-1] = batCap
			elif batChar[ix+1-1] < batMin:
				batChar[ix+1-1] = batMin
				# Distributing remaining electricity accordingly
				for kx in np.arange(0,c2-1):
					if supplyM[ix-1,kx] != -999:
						supplyM[ix-1,kx] = supplyM[ix-1,kx]/sum(supplyM[ix-1,:]*(supplyM[ix-1,:]!=-999))*(batChar[ix-1]-batMin)
	
			#	See if any meters have just fallen short
			#	If they have fallen short, then turn them off
			for kx in np.arange(0,c2-1):
				if	np.sum(supplyM[solDayStart:(ix-1),kx]) >= ELim[kx]:
					supplyM[ix-1,kx] = -999