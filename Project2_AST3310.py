from numpy import *
from matplotlib import pyplot as plt
plt.rcParams['font.family'] = 'serif'
from scipy import interpolate

#Functions needed in the computations---------------------------------------------------------------------

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
   MeV = 1.60218*10**(-13.) #[J]
   
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
   
   N_A_Lamb_pp = 4.01*10**(-15.)*T_9**(-2./3.)*exp(-3.38*T_9**(-1./3.))*(1. + 0.123*T_9**(1./3.) + 1.09*T_9**(2./3.) + 0.938*T_9)
   N_A_Lamb_33 = 6.04e10*T_9**(-2./3.)*exp(-12.276*T_9**(-1./3.))*(1+0.034*T_9**(1./3.) - 0.522*T_9**(2./3.) - 0.124*T_9 + 0.353*T_9**(4./3.) + 0.213*T_9**(-5./3.))
   N_A_Lamb_34 = 5.61e6*T_9_**(5./6.)*T_9**(-3./2.)*exp(-12.826*T_9_**(-1./3.))
   N_A_Lamb_7e = 1.34*10**(-10.)*T_9**(-1./2.)*(1. - 0.537*T_9**(1./3.) + 3.86*T_9**(2./3.) + 0.0027*T_9**(-1.)*exp(2.515*10**(-3.)*T_9**(-1.)))
   N_A_Lamb_71_ = 1.096e9*T_9**(-2./3.)*exp(-8.472*T_9**(-1./3.)) - 4.83e8*T_9__**(5./6.)*T_9**(-3./2.)*exp(-8.472*T_9__**(-1./3.)) + 1.06e10*T_9**(-3./2.)*exp(-30.442*T_9**(-1.))
   N_A_Lamb_71 = 3.11e5*T_9**(-2./3.)*exp(-10.262*T_9**(-1./3.)) + 2.53e3*T_9**(-3./2.)*exp(-7.306*T_9**(-1.))
   
   n_p = rho*X/(1.*m_u)
   n_He = rho*Y/(4.*m_u) #both He_4 and He_3
   n_He_3 = rho*Y_3/(3.*m_u)
   n_He_4 = n_He - n_He_3
   n_Be = rho*Z_Be/(7.*m_u)
   n_Li = rho*Z_Li/(7.*m_u)
   n_e = n_p + 2.*n_He_3 + 2.*n_He_4 + 2.*n_Be + 1.*n_Li

   Lamb_pp = 1e-6*N_A_Lamb_pp/N_A #[m^3/s]
   Lamb_33 = 1e-6*N_A_Lamb_33/N_A #[m^3/s]
   Lamb_34= 1e-6*N_A_Lamb_34/N_A  #[m^3/s]
   if T < 1e6:
      if N_A_Lamb_7e > 1.57e-7/n_e:
         N_A_Lamb_7e = 1.57e-7/n_e
   Lamb_7e = 1e-6*N_A_Lamb_7e/N_A #[m^3/s] #Be + e
   Lamb_71_ = 1e-6*N_A_Lamb_71_/N_A #[m^3/s] #Li + H
   Lamb_71 = 1e-6*N_A_Lamb_71/N_A #[m^3/s] #Be + H

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
 
   eps = r_pp*(Q_pp + Q_Dp) + r_33*Q_33 + r_34*Q_34 + r_7e*Q_7e + r_71_*Q_71_ + r_71*Q_71 #total energy production
   eps1 = (2.*(Q_pp + Q_Dp) + Q_33)*r_33 #PP1 branch
   eps2 = (eps1 + Q_pp + Q_Dp + Q_34)*r_34 + (eps1 + Q_pp + Q_Dp + Q_34 + Q_7e)*r_7e + (eps1 + Q_pp + Q_Dp + Q_34 + Q_7e + Q_71_)*r_71_ #PP2 branch
   eps3 = (eps1 + Q_pp + Q_Dp + Q_34)*r_34 + (eps1 + Q_pp + Q_Dp + Q_34 + Q_71)*r_71 #PP3 branch
   return [eps, eps1/(eps1 +eps2 +eps3), eps2/(eps1 +eps2 +eps3), eps3/(eps1 +eps2 +eps3)]


def calculate_P(rho,T): #calculate pressure
   P_rad = (a*T**4.)/3.   
   P_G = rho*k_B*T/(mu*m_u)
   P = P_rad + P_G
   return P

def calculate_rho(P,T): #calculate density
   P_rad = (a*T**4.)/3.   
   P_G = P - P_rad
   rho = P_G*mu*m_u/(k_B*T)
   return rho

def Hp(_R_,_T_,_M_): #pressure scale height
 return (k_B*_T_*_R_**2.)/(mu*m_u*G*_M_)

def U(_T_, kappa, _rho_, _Hp_, _R_, _M_): #function needed later in calculating temperature gradient
   return ((64.*sigma*_T_**3.)/(3.*kappa*(_rho_/(m_u*mu*N_A))*_rho_*c_p))*((_Hp_*_R_**2.)/(G*_M_*delta))**(1./2.)

def calculate_xi(Hp_, U_,Nabla_ad_,Nabla_stable_): #calculating xi, needed for temperature gradient
   xi_ = roots([(alpha*Hp_)**2/U_, 1., 4.*U_/(alpha*Hp_)**2, Nabla_ad_ - Nabla_stable_]) #find all roots of the polynomial
   for i in range(len(xi_)):
      if xi_[i].imag < 1e-5:
         xi_s = xi_[i].real #choose the real solution
   return xi_s


#Constants and parameters used-------------------------------------------------------------------------------------------------------------------

#fractional abundances by weight
X = 0.7 #H
Y_3 = 1e-10 #He_3
Y = 0.29 #He_3 + He_4
Z = 0.01 #metals
Z_Li = 1e-13 #Li
Z_Be = 1e-13 #Be
mu = 1./(2.*X + 3.*Y/4. + Z/2.)

sigma = 5.67e-8 #[Wm^-2K^-4]
c = 2.998e8 #[m/s]
k_B = 1.382e-23 #[m^2kgs^-2K^-1]
m_u = 1.6605*10**(-27.) #[kg]  
a = (4.*sigma)/c 
G = 6.672e-11 #[Nm^2kg^-2]
sigma = 5.67e-8 #[Wm^-2K^-4]
N_A = 6.022*10**23 #[1/mol]
R_const = 8.314 #[Jmol^-1K^-1]

alpha = 1.
delta = 1.
c_p = 5.*R_const/2.

R_Sun = 6.96e8 #[m]
rho_Sun = 1.408e3 #[kgm^-3], average density of the Sun
M_Sun = 1.989e30 #[kg] #mass of the entire Sun
L_Sun = 3.846e26 #[W]

#First sanity check:---------------------------------------------------------------------------------------------------------------------------
T_sanity = 0.9*10**6 #K
rho_sanity = 55.9 #kgm^-3
R_sanity = 0.84*R_Sun
M_sanity = 0.99*M_Sun
kappa_sanity = 3.98 #m^2kg^-1

Hp_sanity = Hp(R_sanity, T_sanity, M_sanity)
Nabla_stable_sanity = (3.*kappa_sanity*rho_sanity*Hp_sanity*L_Sun)/(64.*pi*R_sanity**2.*sigma*T_sanity**4.)
Nabla_ad_sanity = (2.*alpha)/(5.*delta)
U_sanity = U(T_sanity, kappa_sanity, rho_sanity, Hp_sanity, R_sanity, M_sanity)
xi_s = calculate_xi(Hp_sanity, U_sanity, Nabla_ad_sanity, Nabla_stable_sanity)
Nabla_star_sanity = xi_s**2. + (4*U_sanity/(alpha*Hp_sanity)**2)*xi_s + Nabla_ad_sanity
v_sanity = ((G*M_sanity*Hp_sanity)**(1./2.)/(2.*R_sanity))*xi_s
F_rad_sanity = 16.*sigma*T_sanity**4.*Nabla_star_sanity/(3.*kappa_sanity*rho_sanity*Hp_sanity)
F_C_sanity = 16.*sigma*T_sanity**4.*Nabla_stable_sanity/(3.*kappa_sanity*rho_sanity*Hp_sanity) - F_rad_sanity
frac1_sanity = F_C_sanity/(F_C_sanity + F_rad_sanity)
frac2_sanity = F_rad_sanity/(F_C_sanity + F_rad_sanity)
Nabla_p_sanity = Nabla_star_sanity - xi_s**2.


#Starting point:------------------------------------------------------------------------------------------------------------------------------
#initial values used in the second sanity check, with plots:
rho_0 = 1.42e-7*rho_Sun #[kgm^-3]
T_0 = 5770. #[K]
R_0 = R_Sun
#p = 0.01 #allowed fraction of change

#Best-fit:-----------comment this part to perform sanity check----
rho_0 = 217*rho_0
T_0 = 1.4*T_0
R_0 = 0.85*R_0
M_0 = M_Sun
L_0 = L_Sun
p = 0.075 
#-----------------------------------------------------------------

n = 1e4 #set the max number of datapoints, if mass is negative or step length gets too small -  it cuts the loop and plotting before that
n = int(n)
breaking_point = n #plot the whole array if the mass doesn't get negative

#Create arrays to be updated in the loop
mass = zeros(n)
radius = zeros(n)
pressure = zeros(n)
luminosity = zeros(n)
temperature = zeros(n)
density = zeros(n)
epsilon = zeros(n) #total energy production
epsilon1 = zeros(n) #PP1 branch
epsilon2 = zeros(n) #PP2 branch
epsilon3 = zeros(n) #PP3 branch
Nabla_ad = zeros(n) #adiabatic temperature gradient
Nabla_star = zeros(n)
Nabla_stable = zeros(n)
F_C_ = zeros(n) #fraction of the energy transported by convection
F_rad_ = zeros(n) #fraction of the energy transported by radiation

#Set initial values
mass[0] = M_0
radius[0] = R_0
pressure[0] = calculate_P(rho_0,T_0) 
luminosity[0] = L_0
temperature[0] = T_0
density[0] = rho_0
epsilon[0] = energy_PP(T_0, rho_0)[0] 
epsilon1[0] = energy_PP(T_0, rho_0)[1] 
epsilon2[0] = energy_PP(T_0, rho_0)[2] 
epsilon3[0] = energy_PP(T_0, rho_0)[3] 
Nabla_ad[0] = (2.*alpha)/(5.*delta)
Nabla_stable[0] = (3.*read_kappa(T_0, rho_0)*rho_0*Hp(R_0,T_0,M_0)*L_0)/(64.*pi*R_0**2.*sigma*T_0**4.)
Hp_0 = Hp(R_0,T_0,M_0)
U_0 = U(T_0, read_kappa(T_0,rho_0),rho_0,Hp_0,R_0,M_0)
xi_0 = calculate_xi(Hp_0, U_0,Nabla_ad[0], Nabla_stable[0])
Nabla_star[0] = xi_0**2. + (4*U_0/(alpha*Hp_0)**2)*xi_0 + Nabla_ad[0]
F_rad0 = (16.*sigma*T_0**4.*Nabla_star[0]/(3.*read_kappa(T_0, rho_0)*rho_0*Hp_0))
F_C0 = (16.*sigma*T_0**4.*Nabla_stable[0]/(3.*read_kappa(T_0, rho_0)*rho_0*Hp_0)) - F_rad0
F_C_[0] = F_C0/(F_C0 + F_rad0)
F_rad_[0] = F_rad0/(F_C0 + F_rad0)


#Update arrays - integrate toward the center of the star with mass as an independent variable
for i_ in range(1, n):
      f1 = 1./(4.*pi*(radius[i_-1]**2)*density[i_-1])
      dm1 = p*radius[i_-1]/f1
      f2 = -G*mass[i_-1]/(4.*pi*(radius[i_-1]**4))
      dm2 = p*pressure[i_-1]/f2
      f3 = epsilon[i_-1]
      dm3 = p*luminosity[i_-1]/f3
      kappa = read_kappa(temperature[i_-1], density[i_-1])
      Hp_before = Hp(radius[i_-1], temperature[i_-1], mass[i_-1])
      if Nabla_stable[i_-1] > Nabla_ad[i_-1]: #convective instability, convection
         f4 = -(temperature[i_-1]*Nabla_star[i_-1]*f1)/Hp_before
      else: #if stable, no convection
         f4 = -(temperature[i_-1]*Nabla_stable[i_-1]*f1)/Hp_before
      dm4 = p*temperature[i_-1]/f4
      dm = -1.*(min(abs(dm1), abs(dm2), abs(dm3), abs(dm4))) #choose the smallest mass-step
      mass[i_] = mass[i_-1] + dm
      radius[i_] = radius[i_-1] + dm*f1
      pressure[i_] = pressure[i_-1] + dm*f2
      luminosity[i_] = luminosity[i_-1] + f3*dm
      temperature[i_] = temperature[i_-1] + f4*dm
      density[i_] = calculate_rho(pressure[i_], temperature[i_])
      epsilon[i_] = energy_PP(temperature[i_],density[i_])[0]
      epsilon1[i_] = energy_PP(temperature[i_],density[i_])[1]
      epsilon2[i_] = energy_PP(temperature[i_],density[i_])[2]
      epsilon3[i_] = energy_PP(temperature[i_],density[i_])[3]
      if mass[i_] < 0.: #cut the loop if mass turns negative
         breaking_point = i_ 
         break
      Hp_new = Hp(radius[i_], temperature[i_], mass[i_])
      kappa_new = read_kappa(temperature[i_], density[i_])
      Nabla_stable[i_] = (3.*kappa_new*density[i_]*Hp_new*luminosity[i_])/(64.*pi*radius[i_]**2.*sigma*temperature[i_]**4.)
      Nabla_ad[i_] = Nabla_ad[i_-1]
      U_new = U(temperature[i_], kappa_new, density[i_], Hp_new, radius[i_], mass[i_])
      xi_new = calculate_xi(Hp_new, U_new,Nabla_ad[i_], Nabla_stable[i_])
      if Nabla_stable[i_] > Nabla_ad[i_]: #convection
         Nabla_star[i_] = xi_new**2. + (4*U_new/(alpha*Hp_new)**2)*xi_new + Nabla_ad[i_]
      else: #no convection
         Nabla_star[i_] = Nabla_stable[i_]
      F_rad = 16.*sigma*temperature[i_]**4.*Nabla_star[i_]/(3.*kappa_new*density[i_]*Hp_new)
      F_C = 16.*sigma*temperature[i_]**4.*Nabla_stable[i_]/(3.*kappa_new*density[i_]*Hp_new) - F_rad
      F_C_[i_] = F_C/(F_C + F_rad)
      F_rad_[i_] = F_rad/(F_C + F_rad)

 
#check the radius of the core
for l in range(len(luminosity)):
   if luminosity[l]/L_Sun < 0.995:
      core_starts = l
      break

#check the width of layers
outer_convection_layer = []
outer_radiation_layer = []
inner_radiation_layer = []
inner_convection_layer = []
for k_ in range(len(F_C_[:breaking_point])):
   if luminosity[k_]/L_Sun > 0.995: #outside the core
      if F_C_[k_] > 0.0: #convection
         outer_convection_layer.append(radius[k_]/R_0)
      else: #radiation
         outer_radiation_layer.append(radius[k_]/R_0)
   else: #inside the core
      if F_C_[k_] > 0.0: #convection
         inner_convection_layer.append(radius[k_]/R_0)
      else: #radiation
         inner_radiation_layer.append(radius[k_]/R_0)
outer_convection_layer = array(outer_convection_layer)
inner_convection_layer = array(inner_convection_layer)
outer_radiation_layer = array(outer_radiation_layer)
inner_radiation_layer = array(inner_radiation_layer)


#logarithmic scale for pressure and density
l_pressure = log10(pressure[:breaking_point])
l_density = log10(density[:breaking_point])


#Printouts---------------------------------------------------------------------------------------------------------------------------
print 'First sanity check (in order of example 5.1):'
print Nabla_stable_sanity, Nabla_ad_sanity, Hp_sanity*1e-6,'Mm', U_sanity, xi_s, Nabla_star_sanity, v_sanity, 'm/s', frac1_sanity, frac2_sanity
print Nabla_ad_sanity, '<', Nabla_p_sanity, '<', Nabla_star_sanity, '<', Nabla_stable_sanity

print 'Number of datapoints used:',  breaking_point
print 'Final radius:', (radius[breaking_point-1]/radius[0])*100, '% of initial radius.'
print 'Final mass:', (mass[breaking_point-1]/mass[0])*100, '% of initial mass.'
print 'Final luminosity:',(luminosity[breaking_point-1]/luminosity[0])*100, '% of initial luminosity.'
print 'Radius of the core:', (radius[l]/radius[0])*100, '% of initial radius.'

print 'Width of the outer convection layer:', (outer_convection_layer[0] - outer_convection_layer[-1])*100, '% of initial radius.'
print 'Width of the outer radiation layer:', (outer_radiation_layer[0] - outer_radiation_layer[-1])*100, '% of initial radius.'
print 'Width of the inner radiation layer:', (inner_radiation_layer[0] - inner_radiation_layer[-1])*100, '% of initial radius.'
print 'Width of the inner convection layer:', (inner_convection_layer[0] - inner_convection_layer[-1])*100, '% of initial radius.'

#Plots-------------------------------------------------------------------------------------------------------------------------------
#Plot cross-section
R_values = radius[:breaking_point]/R_Sun
L_values = luminosity[:breaking_point]/L_Sun 
R0 = radius[0]/R_Sun
show_every = 5
core_limit = 0.995

plt.figure()
plt.fig = plt.gcf() # get current figure
ax = plt.gca()  # get current axis
rmax = 1.2*R0
ax.set_xlim(-rmax,rmax)
ax.set_ylim(-rmax,rmax)
ax.set_aspect('equal')	# make the plot circular
j = show_every
for k in range(0, breaking_point-1):
	j += 1
	if j >= show_every:	# don't show every step 
		if(L_values[k] > core_limit):	# outside core
			if(F_C_[k] > 0.0):		# convection
				circR = plt.Circle((0,0),R_values[k],color='firebrick',fill=False)
				ax.add_artist(circR)
			else:				# radiation
				circY = plt.Circle((0,0),R_values[k],color='yellow',fill=False)
				ax.add_artist(circY)
		else:				# inside core
			if(F_C_[k] > 0.0):		# convection
				circB = plt.Circle((0,0),R_values[k],color='navy',fill = False)
				ax.add_artist(circB)
			else:				# radiation
				circC = plt.Circle((0,0),R_values[k],color='c',fill = False)
				ax.add_artist(circC)
		j = 0
circR = plt.Circle((2*rmax,2*rmax),0.1*rmax,color='firebrick',fill=True)	# These are for the legend (drawn outside the main plot)
circY = plt.Circle((2*rmax,2*rmax),0.1*rmax,color='yellow',fill=True)
circC = plt.Circle((2*rmax,2*rmax),0.1*rmax,color='c',fill=True)
circB = plt.Circle((2*rmax,2*rmax),0.1*rmax,color='navy',fill=True)
ax.legend([circR, circY, circC, circB], ['Convection outside core', 'Radiation outside core', 'Radiation inside core', 'Convection inside core']) # only add one (the last) circle of each colour to legend
plt.legend(loc=2)
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$R/R_{\odot}$', fontsize=18)
plt.title('Cross-section of star', fontsize=22)
plt.show()


plt.figure()
plt.plot(radius[:breaking_point]/R_Sun, Nabla_stable[:breaking_point],color='navy' )
plt.plot(radius[:breaking_point]/R_Sun, Nabla_star[:breaking_point], color='firebrick')
plt.plot(radius[:breaking_point]/R_Sun, Nabla_ad[:breaking_point], color='darkcyan')
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$\\nabla$', fontsize=18)
plt.legend(['$\\nabla_{stable}$', '$\\nabla *}$','$\\nabla_{ad}$'], fontsize=16, loc='upper center')
plt.title('Temperature gradients', fontsize=18)
ax = plt.gca()
ax.set_yscale('log')
plt.axis([radius[breaking_point - 1]/R_Sun - 0.01, radius[0]/R_Sun +0.01, 0, max(Nabla_stable)])
plt.savefig('nablas.png')
plt.show()


plt.figure()
plt.scatter(outer_convection_layer*R_0/R_Sun, Nabla_star[:len(outer_convection_layer)],color='firebrick', s=5)
plt.plot(outer_convection_layer*R_0/R_Sun, Nabla_ad[:len(outer_convection_layer)],color='darkcyan')
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$\\nabla$', fontsize=18)
plt.legend(['$\\nabla_{stable}$', '$\\nabla *}$','$\\nabla_{ad}$'], fontsize=16, loc='upper center')
plt.title('Temperature gradients - outer convection layer', fontsize=18)
plt.axis([outer_convection_layer[-1]*R_0/R_Sun, outer_convection_layer[0]*R_0/R_Sun +0.005, 0.395, 0.425])
plt.savefig('nablas_zoom.png')
plt.show()

plt.figure()
plt.plot(radius[:breaking_point]/R_Sun, epsilon1[:breaking_point],color='navy' )
plt.plot(radius[:breaking_point]/R_Sun, epsilon2[:breaking_point], color='firebrick')
plt.plot(radius[:breaking_point]/R_Sun, epsilon3[:breaking_point], color='darkcyan')
plt.plot(radius[:breaking_point]/R_Sun, epsilon[:breaking_point]/max(epsilon), color='gold')
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$\\varepsilon$', fontsize=18)
plt.axis([radius[breaking_point - 1]/R_Sun, radius[0]/R_Sun,0,1])
plt.legend(['$\\varepsilon_{PPI}$', '$\\varepsilon_{PPII}}$','$\\varepsilon_{PPIII}$', '$ \\varepsilon / \\varepsilon_{max} $'], fontsize=16)
plt.title('Energy production', fontsize=18)
plt.savefig('final_epsilon.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, luminosity[:breaking_point]/L_Sun, color='darkorange')
plt.title('Luminosity vs. radius', fontsize=20)
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$L/L_{\odot}$', fontsize=18)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, 0,1.05])
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
plt.title('Temperature vs. radius', fontsize=20)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, T_0/1e6,temperature[breaking_point-1]/1e6])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('T [MK]', fontsize=16)
plt.savefig('final_temperature.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, l_pressure, color='darkcyan')
plt.title('Pressure vs. radius', fontsize=20)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, l_pressure[0], l_pressure[breaking_point-1]])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('log(P) [Pa]', fontsize=18)
plt.savefig('final_pressure.png')
plt.show()


plt.plot(radius[:breaking_point]/R_Sun, mass[:breaking_point]/M_Sun, color='darkgreen')
plt.title('Mass vs. radius', fontsize=20)
plt.axis([radius[breaking_point-1]/R_Sun, radius[0]/R_Sun, 0,1.05])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.ylabel('$M/M_{\odot}$', fontsize=18)
plt.savefig('final_mass.png')
plt.show()

plt.plot(radius[:breaking_point]/R_Sun, F_C_[:breaking_point], color='cyan')
plt.plot(radius[:breaking_point]/R_Sun, F_rad_[:breaking_point], color='blue')
plt.title('Fraction of energy transported by convection and radiation', fontsize=17)
plt.axis([radius[breaking_point-1]/R_Sun-0.004, radius[0]/R_Sun+0.004, -0.05,1.05])
plt.legend(['$F_C/F$','$F_R/F$' ])
plt.xlabel('$R/R_{\odot}$', fontsize=18)
plt.savefig('final_flux.png')
plt.show()



#Printouts - output from the program:
'''
First sanity check (in order of example 5.1):
3.17477413853 0.4 31.6189478301 Mm 602985.321789 0.00118726165962 0.400001409593 65.4432165172 m/s 0.874006341195 0.125993658805
0.4 < 0.400000000003 < 0.400001409593 < 3.17477413853
Number of datapoints used: 337
Final radius: 4.93097592322 % of initial radius.
Final mass: 0.0105530142232 % of initial mass.
Final luminosity: 4.78654626411 % of initial luminosity.
Radius of the core: 29.747970883 % of initial radius.
Width of the outer convection layer: 19.3889555346 % of initial radius.
Width of the outer radiation layer: 49.8865635525 % of initial radius.
Width of the inner radiation layer: 20.3588634326 % of initial radius.
Width of the inner convection layer: 4.15961839927 % of initial radius.
'''
      


