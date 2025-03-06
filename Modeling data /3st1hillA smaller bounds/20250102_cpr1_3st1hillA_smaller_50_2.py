import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

#t = np.linspace(0,34800, num=6960)
t=np.linspace(0,61200,num=12241)
rawdata = np.transpose(np.delete(np.genfromtxt('cpr1longdata_50.csv',delimiter=','),0,0))
newdata = list(rawdata.flatten())
#new = newdata[0:13920:240]
new = newdata[0:12240:240]
lightdata = np.transpose(np.delete(np.genfromtxt('longlight_50.csv', delimiter=','),0,0))

def sqres(ptuple):
    return np.sum((np.asarray(new)-func(t,*ptuple))**2)

#def func(t,Kd,n,d2,k2,k3):
def func(t,d2,k2,k3,k7):
    inivalues = [1,0,0,0,0,0,0]
    arrayvalues = np.asarray([])
   
    for i in range(1):
        def I(t):
            tindex = t/5
            if tindex > 12241:
                tindex = 12240
            return lightdata[int(tindex)]
  
        #def odes(z,t,Kd,n,d2,k2,k3):
        def odes(z,t,d2,k2,k3,k7):
            Pu,Pb,Pa,mRNA,mCherry1,mCherry2,mCherry3 = z
            d1 = 0.019905
            k1 = 0.08299
            Kd = 90.41
            n = 0.964487
            #d2= 486.67 
            #k2= 6.597 
            #k3= 0.0539 


            d3 = 0.000077
            k4 = 1.25
            d4 = 0.000031
            k5 = 0.00283
            k6 = 0.00283

            Pu = z[0]
            Pb = z[1]
            Pa = z[2]
            mRNA = z[3]
            mCherry1 = z[4]
            mCherry2 = z[5]
            mCherry3 = z[6]
            
            dPudt = d1*Pb - k1*I(t)**n/(Kd**n+I(t)**n)*Pu 
            dPbdt = k1*I(t)**n/(Kd**n+I(t)**n)*Pu - k2*Pb - d1*Pb + d2*Pa
            dPadt = k2*Pb - d2*Pa 
            dmRNAdt = k3*Pa - d3*mRNA
            dmCherry1dt = k4*mRNA-(d4 + k5)*mCherry1
            dmCherry2dt = k5*mCherry1-(d4+k6)*mCherry2
            dmCherry3dt = k6*mCherry2 - d4*mCherry3



            return [dPudt,dPbdt,dPadt,dmRNAdt,dmCherry1dt,dmCherry2dt,dmCherry3dt]  


        #solver = odeint(odes,inivalues,t,args = (Kd,n,d2,k2,k3),hmax=1)
        solver = odeint(odes,inivalues,t,args = (d2,k2,k3,k7),hmax=0.1)

        mCherryout = solver[:,6]
        #mCherry = mCherryout[0:13920:240]
        mCherry = mCherryout[0:12240:240]
        arrayvalues = np.hstack((arrayvalues,mCherry))  
    return list(arrayvalues)