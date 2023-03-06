# Gerardo Cano Perea  
"""
Created on Sun Oct 11 22:40:23 2020
AirfoilPrep.py from NREL / Airfoil Data Extrapolation.
Before Run the Script Check the Installation. 
"""
import numpy as np
import pandas as pd
from airfoilprep import Polar
from airfoilprep import Airfoil

Ref_file = 'NACA_63824' # Reference CSV File 
Data = pd.read_csv(Ref_file +'.csv') #Polar File 
Re = 50000 # Reynolds Number - 75%. [-]
alpha = np.array(Data.AoA) # Angle of Attack. [-]
cl = np.array(Data.CL) # Lift Coefficient. [-]
cd = np.array(Data.CD) # Drag Coefficient. [-] 
cm = np.zeros([len(Data.AoA)]) # Momentum Coefficient. [-] 
# If the cm was not calculated then is replaced with a zeros vector. 

# Creating Polar and Airfoil Objects. 
p1 = Polar(Re, alpha, cl, cd, cm) # First Set of Polar Data. 
af = Airfoil([p1]) # Contains all Sets of Polar Data.  

AR = 10 # Aspect Ratio Value. [-]
cdmax = 1.11 + 0.018 * AR # Maximum Drag Coefficient [-] 
cdmin = 0.001  # minimum drag coefficient.  Viterna's method can occasionally produce
               # negative drag coefficients.  A minimum is used to prevent unphysical data.
               # The passed in value is used to override the defaul. 
af_extrapolate = af.extrapolate(cdmax, AR=AR, cdmin=cdmin) # Sub-routine for Extrapolation 

# Plotting Airfoil Data  
af.plot(single_figure=False)
af_extrapolate.plot(single_figure=False)

af_extrapolate.writeToAerodynFile(Ref_file + '_Extended.dat')
