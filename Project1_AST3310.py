from numpy import *
from matplotlib import pyplot as plt
plt.rcParams['font.family'] = 'serif'
from scipy import interpolate


dynamic_steplength = 'on' #here turn the dynamic step size on and off <-----------------------------------------------------------------


def read_kappa(T, rho): #function reading opacity.txt, interpolating or extrapolating values if necessary
   rho_cgs = 1e-3*rho
   R = 1e18*rho_cgs*T**(-3.)
   infile = open('opacity.txt', 'r')
   rows = infile.readlines()
   infile.close()
   log_R = rows[0].split()
   log_R = log_R[1:]
   log_R = [float(logR) for logR in log_R]
   
   log_T = []
   log_kappa = []
   for r in rows[2:]: #skip log_R row and a blank row
      r_list = r.split()
      r_list = [float(rs) for rs in r_list]
      temperature = r_list[0]
      kappas = r_list[1:]
      log_T.append(temperature)
      log_kappa.append(kappas)
   
   log_T = array(log_T) #[K]
   log_R = array(log_R) #cgs
   log_kappa = array(log_kappa) #cgs
   
   T_log = log10(T)
   R_log = log10(R)
   
   #when the input value can be found in the table:
   for i in range(len(log_T)):  
      if T_log == log_T[i]:
         kappa_1index = i
         #print "log(T) found in the table"
         break
      else:
         kappa_1index = 'not found'
   for i in range(len(log_R)):  
      if R_log == log_R[i]:
         kappa_2index = i
         #print "log(R) found in the table"
         break
      else:
         kappa_2index = 'not found'
   
   #Extrapolation warning
   if T_log < 3.75 or T_log > 8.7:
      print 'Extrapolation: log(T) outside the bounds of the table'
   if R_log < -8.0 or R_log > 1.:
      print 'Extrapolation: log(R) outside the bounds of the table'
   
   
   if kappa_1index != 'not found' and kappa_2index != 'not found':
      log_kappa_value = log_kappa[kappa_1index][kappa_2index]
      kappa_value = 10**log_kappa_value #cgs
   else: 
      log_kappa_value_function = interpolate.interp2d(log_R, log_T, log_kappa, kind='linear', copy=False, bounds_error=False, fill_value=None)
      log_kappa_value = log_kappa_value_function(R_log, T_log)
      kappa_value = 10**log_kappa_value[0] #cgs
   return kappa_value/10. #kappa (SI)

     
  

def energy_PP(T, rho): #function returning energy generation through PP chain reactions
   #fractional abundances by weight
   X = 0.7 #H
   Y_3 = 1e-10 #He_3
   Y = 0.29 #He_3 + He_4
   Z = 0.01 #metals
   Z_Li = 1e-7 #Li
   Z_Be = 1e-7 #Be
   
   MeV = 1.60218*10**(-13.) #[J]
   m_u = 1.6605*10**(-27.) #[kg]
   Q_pp = (0.15 + 1.02)*MeV #[J]
   Q_Dp = 5.49*MeV #[J]
   Q_33 = 12.86*MeV #[J]
   Q_34 = 1.59*MeV #[J]
   Q_7e = 0.05*MeV #[J]
   Q_71_ = 17.35*MeV #[J]
   Q_71 = (0.14 + 1.02 + 6.88 + 3.0)*MeV #[J]
 

   T_9 = T*10**(-9.) # temperature conversion to units 10^9K
   T_9_ = T_9/(1+4.95*10**(-2.)*T_9)
   T_9__ = T_9/(1+0.759*T_9)
   N_A = 6.022*10**23 #[1/mol]
   
   N_A_Lamb_pp = 4.01*10**(-15.)*T_9**(-2./3.)*exp(-3.38*T_9**(-1./3.))*(1. + 0.123*T_9**(1./3.) + 1.09*T_9**(2./3.) + 0.938*T_9)
   N_A_Lamb_33 = 6.04e10*T_9**(-2./3.)*exp(-12.276*T_9**(-1./3.))*(1+0.034*T_9**(1./3.) - 0.522*T_9**(2./3.) - 0.124*T_9 + 0.353*T_9**(4./3.) + 0.213*T_9**(-5./3.))
   N_A_Lamb_34 = 5.61e6*T_9_**(5./6.)*T_9**(-3./2.)*exp(-12.826*T_9_**(-1./3.))
   N_A_Lamb_7e = 1.34*10**(-10.)*T_9**(-1./2.)*(1. - 0.537*T_9**(1./3.) + 3.86*T_9**(2./3.) + 0.0027*T_9**(-1.)*exp(2.515*10**(-3.)*T_9**(-1.)))
   N_A_Lamb_71_ = 1.096e9*T_9**(-2./3.)*exp(-8.472*T_9**(-1./3.)) - 4.83e8*T_9__**(5./6.)*T_9**(-3./2.)*exp(-8.472*T_9__**(-1./3.)) + 1.06e10*T_9**(-3./2.)*exp(-30.442*T_9**(-1.))
   N_A_Lamb_71 = 3.11e5*T_9**(-2./3.)*exp(-10.262*T_9**(-1./3.)) + 2.53e3*T_9**(-3./2.)*exp(-7.306*T_9**(-1.))

   Lamb_pp = 1e-6*N_A_Lamb_pp/N_A #[m^3/s]
   Lamb_33 = 1e-6*N_A_Lamb_33/N_A #[m^3/s]
   Lamb_34= 1e-6*N_A_Lamb_34/N_A  #[m^3/s]
   if T < 1e6:
      if N_A_Lamb_7e > 1.57e-7/n_e:
         N_A_Lamb_7e = 1.57e-7/n_e
   Lamb_7e = 1e-6*N_A_Lamb_7e/N_A #[m^3/s] #Be + e
   Lamb_71_ = 1e-6*N_A_Lamb_71_/N_A #[m^3/s] #Li + H
   Lamb_71 = 1e-6*N_A_Lamb_71/N_A #[m^3/s] #Be + H

   n_p = rho*X/(1.*m_u)
   n_He = rho*Y/(4.*m_u) #both He_4 and He_3
   n_He_3 = rho*Y_3/(3.*m_u)
   n_He_4 = n_He - n_He_3
   n_Be = rho*Z_Be/(7.*m_u)
   n_Li = rho*Z_Li/(7.*m_u)
   n_e = n_p + 2.*n_He_3 + 2.*n_He_4 + 2.*n_Be + 1.*n_Li

   r_pp = n_p*n_p*Lamb_pp/(rho*2.)
   r_33 = n_He_3*n_He_3*Lamb_33/(rho*2.)
   r_34 = n_He_3*n_He_4*Lamb_34/rho
   if r_pp < (r_33*2. + r_34):
      rate1 = r_pp/(2.*r_33 + r_34)
      r_33 *= rate1
      r_34 *= rate1
      
   r_7e = n_Be*n_e*Lamb_7e/rho
   r_71_ = n_Li*n_p*Lamb_71_/rho
   r_71 = n_Be*n_p*Lamb_71/rho
   if r_34 < (r_7e + r_71):
      rate2 = r_34/(r_7e + r_71)
      r_7e *= rate2
      r_71 *= rate2
   if r_7e < r_71_:
      rate3 = r_7e/r_71_
      r_71_ *= rate3
 
 
   #return r_pp*(Q_pp + Q_Dp)*rho, r_33*Q_33*rho, r_34*Q_34*rho, r_7e*Q_7e*rho, r_71_*Q_71_*rho, r_71*Q_71*rho
   eps = r_pp*(Q_pp + Q_Dp) + r_33*Q_33 + r_34*Q_34 + r_7e*Q_7e + r_71_*Q_71_ + r_71*Q_71
   return eps


def calculate_P(rho,T): #calculate pressure
   sigma = 5.67e-8 #[Wm^-2K^-4]
   c = 2.998e8 #[m/s]
   k_B = 1.382e-23 #[m^2kgs^-2K^-1]
   m_u = 1.6605*10**(-27.) #[kg]
   a = (4.*sigma)/c
   Z = 0.01
   X = 0.7
   Y = 0.29
   mu = 1./(2.*X + 3.*Y/4. + Z/2.)
   P_rad = (a*T**4.)/3.   
   P_G = rho*k_B*T/(mu*m_u)
   P = P_rad + P_G
   return P

def calculate_rho(P,T): #calculate density
   sigma = 5.67e-8 #[Wm^-2K^-4]
   c = 2.998e8 #[m/s]
   k_B = 1.382e-23 #[m^2kgs^-2K^-1]
   m_u = 1.6605*10**(-27.) #[kg]
   a = (4.*sigma)/c
   Z = 0.01
   X = 0.7
   Y = 0.29
   mu = 1./(2.*X + 3.*Y/4. + Z/2.)
   P_rad = (a*T**4.)/3.   
   P_G = P - P_rad
   rho = P_G*mu*m_u/(k_B*T)
   return rho
   

R_Sun = 6.96e8 #[m]
rho_Sun = 1.408e3 #[kgm^-3], average density of the Sun
M_Sun = 1.989e30 #[kg] #mass of the entire Sun
L_Sun = 3.846e26 #[W]
G = 6.672e-11 #[Nm^2kg^-2]
sigma = 5.67e-8 #[Wm^-2K^-4]

#Starting point - best-fit parameters:
rho_0 = 5.9e3 #[kgm^-3]
T_0 = 8.6e6 #[K]
R_0 = 0.3711*R_Sun
M_0 = 0.8*M_Sun

if dynamic_steplength == 'off':
   M_final = 0.
   dm = -1e-4*M_Sun #vary fixed steplength here <-----------------
   n = 1 + (M_0 - M_final)/abs(dm) #number of datapoints
   n = int(n)
   mass = linspace(M_0, M_final, n)
   
else:
   n = 1e4 #set the max number of datapoints, if mass is negative or step length gets too small -  it cuts the loop and plotting before that
   n = int(n)
   mass = zeros(n)
   mass[0] = M_0

radius = zeros(n)
pressure = zeros(n)
luminosity = zeros(n)
temperature = zeros(n)
density = zeros(n)
epsilon = zeros(n)

radius[0] = R_0
pressure[0] = calculate_P(rho_0,T_0) 
luminosity[0] = L_Sun
temperature[0] = T_0
density[0] = rho_0
epsilon[0] = energy_PP(T_0, rho_0)

breaking_point = n #plot the whole array if the mass doesn't get negative

print_counter = 0.
for i_ in range(1, n):
   if dynamic_steplength == 'off':
      radius[i_] = radius[i_-1] + dm*(1./(4.*pi*(radius[i_-1]**2)*density[i_-1]))
      pressure[i_] = pressure[i_-1] + dm*((-G*mass[i_-1])/(4.*pi*(radius[i_-1]**4)))
      luminosity[i_] = luminosity[i_-1] + epsilon[i_-1]*dm
      kappa = read_kappa(temperature[i_-1], density[i_-1])
      temperature[i_] = temperature[i_-1] + (-3.*kappa*luminosity[i_-1]/(256.*pi**2*sigma*radius[i_-1]**4.*temperature[i_-1]**3))*dm
      density[i_] = calculate_rho(pressure[i_], temperature[i_])
      epsilon[i_] = energy_PP(temperature[i_],density[i_])
   else:
      p = 0.1 #allowed fraction of change
      f1 = 1./(4.*pi*(radius[i_-1]**2)*density[i_-1])
      dm1 = p*radius[i_-1]/f1
      f2 = -G*mass[i_-1]/(4.*pi*(radius[i_-1]**4))
      dm2 = p*pressure[i_-1]/f2
      f3 = epsilon[i_-1]
      dm3 = p*luminosity[i_-1]/f3
      kappa = read_kappa(temperature[i_-1], density[i_-1])
      f4 = -3.*kappa*luminosity[i_-1]/(256.*pi**2*sigma*radius[i_-1]**4.*temperature[i_-1]**3)
      dm4 = p*temperature[i_-1]/f4
      dm = -1.*(min(abs(dm1), abs(dm2), abs(dm3), abs(dm4)))
      mass[i_] = mass[i_-1] + dm
      radius[i_] = radius[i_-1] + dm*f1
      pressure[i_] = pressure[i_-1] + dm*f2
      luminosity[i_] = luminosity[i_-1] + f3*dm
      temperature[i_] = temperature[i_-1] + f4*dm
      density[i_] = calculate_rho(pressure[i_], temperature[i_])
      epsilon[i_] = energy_PP(temperature[i_],density[i_])
   if mass[i_] < 0.:
      breaking_point = i_ #cut the loop if mass turns negative
      break
   if abs(dm) < 1e-6*M_0: #break before the step gets so low that the majority of values is around the final one
      breaking_point = i_ 
      break
   print_counter += 1.
   if print_counter > 10.: #print every 10th value
      #print dm
      #print 'M/M_Sun:', mass[i_]/M_Sun, 'R/R_Sun:', radius[i_]/R_Sun, 'P:', pressure[i_], 'L/L_Sun:', luminosity[i_]/L_Sun, 'kappa:', kappa, 'T:', temperature[i_], \
      #'rho/rho_Sun:', density[i_]/rho_Sun, 'epsilon:', epsilon[i_]
      print_counter = 0.

 
print 'Number of datapoints used:',  breaking_point

#logarithmic scale
l_eps = log10(epsilon[:breaking_point])
l_pressure = log10(pressure[:breaking_point])
l_density = log10(density[:breaking_point])


#check the radius of the core
for l in range(len(luminosity)):
   if luminosity[l]/L_Sun < 0.995:
      core_starts = l
      print 'Radius of the core:', (radius[l]/radius[0])*100, '% of initial radius.'
      break
      

#print final values of parameters
print mass[breaking_point-1], radius[breaking_point-1], luminosity[breaking_point-1]
print 'final radius:', (radius[breaking_point-1]/radius[0])*100, '%', 'final mass:', (mass[breaking_point-1]/mass[0])*100, '%', 'final luminosity:',(luminosity[breaking_point-1]/luminosity[0])*100, '%'


#Final plots

plt.plot(radius[:breaking_point]/R_Sun, l_eps, color='firebrick')
plt.title('Energy production vs. radius', fontsize=20)
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$log(\epsilon)$ $[Jkg^{-1}s^{-1}]$', fontsize=18)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, l_eps[0], l_eps[breaking_point-1]])
plt.savefig('final_epsilon.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, l_pressure, color='darkcyan')
plt.title('Pressure vs. radius', fontsize=20)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, l_pressure[0], l_pressure[breaking_point-1]])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('log(P) [Pa]', fontsize=18)
plt.savefig('final_pressure.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, mass[:breaking_point]/M_Sun, color='darkgreen')
#plt.scatter(radius[:breaking_point]/R_Sun, mass[:breaking_point]/M_Sun, color='darkgreen')
plt.title('Mass vs. radius', fontsize=20)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, 0,0.85])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$M/M_{\odot}$', fontsize=18)
plt.savefig('final_mass.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, luminosity[:breaking_point]/L_Sun, color='darkorange')
plt.title('Luminosity vs. radius', fontsize=20)
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$L/L_{\odot}$', fontsize=18)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, 0,1])
plt.axvline(x=radius[core_starts]/R_Sun, linestyle='dashed', color='gold')
plt.savefig('final_luminosity.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, l_density, color='slateblue')
plt.title('Density vs. radius', fontsize=20)
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$log(\\rho)$ $[kgm^{-3}]$', fontsize=18)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, l_density[0],l_density[breaking_point-1]])
plt.savefig('final_density.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, temperature[:breaking_point]/1e6, color='darkred')
#plt.scatter(radius[:breaking_point]/R_Sun, temperature[:breaking_point]/1e6, color='cyan')
plt.title('Temperature vs. radius', fontsize=20)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, T_0/1e6,temperature[breaking_point-1]/1e6])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('T [MK]', fontsize=16)
plt.savefig('final_temperature.png')
plt.show()

print 'initial pressure:', calculate_P(rho_0, T_0)
print 'final temperature:', temperature[breaking_point-1]
print 'final pressure:', pressure[breaking_point-1]
print 'final density:', density[breaking_point-1]
print 'initial epsilon:', epsilon[0]
print 'final epsilon:', epsilon[breaking_point-1]


#dm-static test
'''
plt.plot(mass[:breaking_point]/M_Sun, radius[:breaking_point]/R_Sun, color='gold')
#plt.scatter(mass[:breaking_point]/M_Sun, radius[:breaking_point]/R_Sun, color='darkgreen')
plt.title('Radius vs. mass, $\Delta m = -10^{-4}M_{\odot}$', fontsize=20)
plt.axis([0, 0.8, radius[breaking_point-1]/R_Sun, radius[0]/R_Sun])
plt.xlabel('$M/M_{\odot}$', fontsize=18)
plt.ylabel('$R/R_{\odot}$', fontsize=18)
plt.savefig('step_size_static.png')
plt.show()
'''
'''
Printout of the program:
Number of datapoints used: 91
Radius of the core: 81.9376803913 % of initial radius.
2.74213201217e+28 2731223.80358 6.08872121416e+23
final radius: 1.0574433122 % final mass: 1.72331071655 % final luminosity: 0.158313084092 %
initial pressure: 6.86558833929e+14
final temperature: 20049659.2927
final pressure: 3.02119130176e+17
final density: 1115729.14938
initial epsilon: 6.5524279599e-06
final epsilon: 0.0427071269495
'''
