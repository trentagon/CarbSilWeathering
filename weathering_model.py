###############################################################################
# 7/9/18
# 
# This file contains all of the functions necessary to simulate pCO2 in the 
# atmosphere of a planet with liquid water at the surface. The model is written
# using the equations of "Constraining the climate and ocean pH of the early
# Earth with a geological carbon cycle model" by Krissansen-Totten et al. 
# (2018). Unless otherwise noted, all equation references in this code refer 
# to the equations of that paper, which is abbreviated as JKT in this code.
#
# This code was written by Owen Lehmer, questions can be addressed to:
#     owen.r.lehmer@nasa.gov
# 
#                            RUNNING THE MODEL
# Steps: 
# Simply call runWeatheringModel() to run the model. See the end of this file 
# for an example call to runWeatheringModel(). Uncomment that example and run
# this file to see the model in action. It's also recommended that you set DEBUG
# to True for additional model output.
###############################################################################

from math import exp, log, log10
from scipy.integrate import ode
from scipy.optimize import fsolve
import numpy as np

DEBUG = False #show prints during running when True

class ModelInputs:
    """
    This class represents the structure for changing model parameters. The
    default parameters will be for the modern Earth. The values were taken from
    the mean values of the parameter ranges given in Table S1 and Table S2 of
    JKT. To change the parameters create a new ModelInputs object and pass it to 
    runWeatheringModel().
    """
    def __init__(self):
        self.Mo = 1.35E21 #ocean mass [kg]
        self.f_bio = 1.0 #biological weathering fraction 
        self.f_land = 1.0 #land fraction compared to modern Earth
        self.oceanpH = 8.2 #pH of modern ocean
        self.Hmod_mol = 10.0**(-self.oceanpH) #equation S16, initial H conc.
        self.pCO2 = 0.000280 #pCO2 on modern Earth [bar]
        self.CO2_alpha = 0.3 #alpha term in equation 1
        self.CO2_xi = 0.3 #xi term in equation S2
        self.Te = 25.0 #e-folding temp in equations 1, S2 [K]
        self.Fmod_out = 6.0E12 #modern outgassing rate [mol C yr-1]
        self.Fmod_carb = 10.0E12 #modern carbonate weathering rate [mol C yr-1]
        self.carb_n = 1.75 #carbonate precipitation coefficient 
        self.diss_x = 1.0 #modern seafloor dissolution relative to prec.
        self.grad = 1.075 #temperature gradient from surface to ocean depth
        self.gamma = 0.2 #pH dependence of seafloor weathering
        self.E_bas = 90000.0 #temp dependence of seafloor weathering [J mol-1]
        self.beta = 0.1 #spreading rate dependence
        self.out_m = 1.5 #outgassing exponent
        self.sed_depth = 1.0 #sediment thickness relative to modern Earth
        self.Pmod_pore = 0.45E12 #modern pore space precipitation [mol C yr-1]
        self.Fmod_diss = 0.45E12 #modern seafloor dissolution rate [mol C yr-1]
        self.Ts_mod = 288.0 #modern (preindustrial) surface temp [K]
        self.ca = 0.01 #modern Ca abundance [mol kg-1]
        self.s = 1.8E20/self.Mo #correction factor for mass balance
                                #the 1.8E20 is the mass of the atmosphere
        self.Q = 1.0 #internal heat flow compared to modern
        self.K = 77.8 #conductivity of the sediments [m K-1]
        self.lum = 1.0 #the luminosity compared to the modern Earth

    def print(self):
        for key, val in self.__dict__.items():
            print("%s = %2.3e"%(key, val))


def runWeatheringModel(inputs=ModelInputs(), guess=None, chem_tol=1.0E-5, 
        chem_max_iter=10):
    """
    This is the top level function to run the model. The model inputs can be
    changed by creating an instance of the ModelInputs() class then changing
    the values as desired, i.e.:
        my_inputs = ModelInputs()
        my_inputs.oceanpH = 7.0

    To run the model for the modern Earth simply call runWeatheringModel() with
    no parameters. 

    The model can be run with a fast numerical solver (scipy.optimize.fsolve)
    or via manual integration (very slow). If the model is not converging try
    running the manual mode.

    Inputs:
        inputs        - the model inputs parameters, default to the modern Earth
        guess         - a guess for the alkalinity and carbon concentrations in
                        the system. The guess must be an array of the form:
                        [Co,Ao,Cp,Ap]
        chem_tol      - the tolerance (% difference between runs) that will be 
                        used with the carbon chemistry
        chem_max_iter - the maximum number of iterations allowed during 
                        chemistry calculations

    Returns:
        [pCO2, pH, Ts] - array of values where the values are:
                            pCO2   - the atmospheric CO2 partial pressure [bar]
                            pH     - the ocean pH
                            Ts     - the surface temperature [K]
        times  - the time values that the model was run at
        y_vals - an array with the results of the model in the form:
                 [Co, Ao] where Co is the ocean dissolved inorganic carbon 
                 concentration Ao is the alkalinity
        status - the status of the model run: 
                 -1 = solver had an error during integration
                  0 = model reached 10 Gyr without converging
                  1 = model converged successfully 
    """

    #calculate all the parameters needed
    F_out = globalOutgassing(inputs.Fmod_out, inputs.Q, inputs.out_m)

    [K_diss, K_ocean, K_pore, Fmod_sil, Co0, Ao0] = initialParameters(inputs)

    if guess is not None:
        #was given a guess
        [Co0, Ao0] = guess

    #parameters we want to record, initialize them here
    pCO2 = inputs.pCO2
    pH = inputs.oceanpH
    Ts = 0.0 #just initialize the variable here

    def calculateModelDerivatives(_, params):
        """
        The equation that will be passed to the ODE solver. The solver will take 
        initial guesses for Co, and Ao then iterate on each variable 
        until convergence is reached or the maximum number of iterations is 
        reached.
        """
        nonlocal pCO2, pH, Ts

        [Co, Ao] = params #Co is the ocean DIC, Ao is the ocean alkalinity

        #correct the ocean DIC for mass balance
        DIC = Co - inputs.s*pCO2
        

        Ts_old = 0.0
        num = 0 #track the number of runs
        diff = 1.0 #track the relative difference in surface temperature
        omega = 0.0
        T_pore = 0.0
        while diff > chem_tol and num < chem_max_iter:
            """
            Solve for the equilibrium carbon chemistry. The solution is 
            temperature dependent, so we'll iterate the process until the 
            temperature and chemistry become steady (measured by the % change
            in Ts between runs).
            """
            
            Ts = surfaceTempFromClimateModel(pCO2, inputs.lum) 
            T_do = deepOceanTemperature(Ts, inputs.grad)
            T_pore = poreSpaceTemperature(T_do, inputs.Q, inputs.sed_depth, 
                                          inputs.K)

            omega, pCO2, pH = equilibriumChemistry(Ts, Ao, DIC, inputs.s, 
                                                     Ao0, inputs.ca)


            diff = np.abs(Ts-Ts_old)/Ts 
            num += 1
            Ts_old = Ts 

        deltaTs = Ts - inputs.Ts_mod

        #calculate the terms needed
        P_ocean = 0
        P_pore = 0
        #only precipitate when omega > 1
        if omega > 1:
            P_ocean = K_ocean*inputs.f_land*(omega - 1.0)**inputs.carb_n #S19
            P_pore = K_pore*(omega - 1.0)**inputs.carb_n #S19

        F_carb = continentalCarbonateWeathering(inputs.f_bio, inputs.f_land, 
            inputs.Fmod_carb, pCO2, inputs.pCO2, inputs.CO2_xi, 
            deltaTs, inputs.Te)

        F_sil = continentalSilicateWeahtering(inputs.f_bio, inputs.f_land, 
                Fmod_sil, pCO2, inputs.pCO2, inputs.CO2_alpha, deltaTs, 
                inputs.Te)

        r_sr = spreadingRate(inputs.Q, inputs.beta)
        H_mol = 10**(-pH)
        F_diss = seafloorBasaltDissolution(K_diss, r_sr, inputs.E_bas, 
                T_pore, H_mol, inputs.Hmod_mol, inputs.gamma)

        derivatives = carbonAndAlkalinityEquations(inputs.Mo, F_out, 
                F_carb, P_ocean, F_sil, F_diss, P_pore)

        return derivatives #derivatives = [dtCo, dtAo]

    #prepare the model inputs
    initial_val = [Co0, Ao0]
       
    r = ode(calculateModelDerivatives).set_integrator('dopri5', method='bdf')
    r.set_initial_value(initial_val, 0)

    old_co2 = 0
    dt = 10000
    cont = True
    status = 0
    while cont and r.successful() and r.t < 1.0E10:
        r.integrate(r.t + dt)

        if r.t > 1.0E5:
            #make sure nothing odd happens with the initial conditions
            slope = (pCO2 - old_co2)/dt
            change = (pCO2 + slope*10**9)/pCO2
            if change < 1.01: 
                #less than a 1% relative change, 
                #return a positive status value to end
                cont = False 
                status = 1
                if DEBUG:
                    print("converged at time: %2.3e"%(r.t))
            old_co2 = pCO2
    if not r.successful():
        if DEBUG:
            print("Failed to converge with inputs:")
            inputs.print()
        status = -1


    if DEBUG:
        print("Solver finished with pCO2=%2.2e, Ts=%0.2f, pH=%0.2f"%
                (pCO2, Ts, pH))
        inputs.print()

    return [pCO2, pH, Ts], status


def initialParameters(inputs):
    """
    Several parameters need to be initialized. Here we use the (assumed) values
    for modern seafloor carbonate precipitation, outgassing, and the ratio of
    seafloor dissolution to carbonate precipitation to calculate the rate 
    constants needing by the model.
    """
    T_s = 288 #we initialize from the modern Earth
    T_do = deepOceanTemperature(T_s, inputs.grad)
    T_pore = poreSpaceTemperature(T_do, inputs.Q, inputs.sed_depth, inputs.K)

    [K1, K2, H_CO2] = equilibriumRateConstants(T_s)

    partition = inputs.Fmod_diss/inputs.Fmod_out
    Fmod_diss = partition*inputs.Fmod_out*inputs.diss_x
    Fmod_sil = (1.0 - partition)*inputs.Fmod_out + \
            (1 - inputs.diss_x)*partition*inputs.Fmod_out
    Pmod_pore = partition*inputs.Fmod_out

    #initial conditions for atmosphere-ocean system (eqns S12 to S14)
    CO2aq_o = H_CO2*inputs.pCO2
    HCO3_o = K1*CO2aq_o/(10**-inputs.oceanpH)
    CO3_o = K2*HCO3_o/(10**-inputs.oceanpH)
    DIC_o = CO3_o + HCO3_o + CO2aq_o #total dissolved inorganic carbon
    ALK_o = 2.0*CO3_o + HCO3_o #carbonate alkalinity

    #assume steady state for modern, so ocean precip is equal to inputs minus pore
    Pmod_ocean = inputs.Fmod_out + inputs.Fmod_carb - Pmod_pore

    omega_o = inputs.ca*CO3_o/carbonateSolubility(T_do)

    omega_p = inputs.ca*CO3_o/carbonateSolubility(T_pore)

    K_pore = Pmod_pore/(omega_p - 1.0)**inputs.carb_n #for pore precip.
    K_ocean = Pmod_ocean/(omega_o - 1.0)**inputs.carb_n #for ocean precip.

    K_diss = Fmod_diss/(2.88*10**-14*10**(-inputs.gamma*inputs.oceanpH)*
                   exp(-inputs.E_bas/(8.314*T_pore)))

    Co = DIC_o + inputs.pCO2*inputs.s
    Ao = ALK_o

    return [K_diss, K_ocean, K_pore, Fmod_sil, Co, Ao]


def carbonAndAlkalinityEquations(Mo, F_out, F_carb, P_ocean, F_sil, F_diss, 
        P_pore):
    """
    This function calculates the derivatives of concentrations of carbon in the 
    ocean-atmosphere system. It will also solve for the alkalinity of the ocean. 
    These functions corresponds to equation S1 without the pore space. This is 
    the top-level function that represents the model described by 
    Krissansen-Totten et al. (2018).

    Inputs:
        Mo      - mass of the ocean [kg]
        F_out   - global outgassing flux [mol C yr-1]
        F_carb  - continental carbonate weathering rate [mol C yr-1]
        P_ocean - precipitation flux of carbonates in the ocean [mol C yr-1]
        F_sil   - continental silicate weathering flux [ mol C yr-1]
        F_diss  - seafloor weathering from basalt dissolution [mol eq yr-1]
        P_pore  - carbonate precipitation flux in the pore space [mol C yr-1]

    Returns:
        dCo_dt - change in atmosphere-ocean carbon concentration [mol C yr-1]
        dAo_dt - change in atmosphere-ocean alkalinity [mol eq yr-1]
    """

    dCo_dt = (F_out + F_carb - P_ocean - P_pore)/Mo
    dAo_dt = 2*(F_sil + F_carb - P_ocean + F_diss - P_pore)/Mo

    return [dCo_dt, dAo_dt]


def continentalCarbonateWeathering(f_bio, f_land, Fmod_carb, pCO2, pCO2mod,
        eps, deltaTs, Te):
    """
    The rate of carbon liberated by continental weathering will be returned from 
    this function. This is function S2 in JKT.

    Inputs:
        f_bio     - biological enhancement of weathering, set to 1 for the 
                    modern Earth [dimensionless]
        f_land    - land fraction compared to the modern Earth [dimensionless]
        Fmod_carb - Earth's modern carbonate weathering rate [mol yr-1]
        pCO2      - partial pressure of CO2 [Pa]
        pCO2mod   - Earth's modern (preindustrial) CO2 partial pressure [Pa]
        eps       - empirical constant [dimensionless]
        deltaTs   - difference in global mean surface temperature [K]
        Te        - defines temperature dependence of weathering [K]

    Returns:
        F_carb - the carbonate weathering rate [mol yr-1]
    """

    F_carb = f_bio*f_land*Fmod_carb*(pCO2/pCO2mod)**eps
    if Te > 0:
        F_carb = f_bio*f_land*Fmod_carb*(pCO2/pCO2mod)**eps*exp(deltaTs/Te)

    return F_carb

def continentalSilicateWeahtering(f_bio, f_land, Fmod_sil, pCO2, pCO2mod, 
        alpha, deltaTs, Te):
    """
    The rate of silicate weathering  from continents. This function corresponds
    to equation 1 from JKT.

    Inputs:
        f_bio    - biological enhancement of weathering, set to 1 for the 
                   modern Earth [dimensionless]
        f_land   - land fraction compared to the modern Earth [dimensionless]
        Fmod_sil - Earth's modern silicate weathering rate [mol yr-1]
        pCO2     - partial pressure of CO2 [Pa]
        pCO2mod  - Earth's modern (preindustrial) CO2 partial pressure [Pa]
        alpha    - empirical constant [dimensionless]
        deltaTs  - difference in global mean surface temperature [K]
        Te       - defines temperature dependence of weathering [K]

    Returns:
        F_sil - the silicate weathering rate [mol yr-1]
    """

    F_sil = f_bio*f_land*Fmod_sil*(pCO2/pCO2mod)**alpha
    if Te > 0: 
        F_sil = f_bio*f_land*Fmod_sil*(pCO2/pCO2mod)**alpha*exp(deltaTs/Te)

    return F_sil

def seafloorBasaltDissolution(k_diss, r_sr, E_bas, T_pore, H_mol, Hmod_mol,
        gamma):
    """
    This function will calculate the rate of basalt dissolution on the 
    seafloor. This function represents equation S3 of JKT.

    Inputs:
        k_diss   - proportionality constant chosen to match modern flux 
                   [dimensionless]
        r_sr     - spreading rate compared to modern [dimensionless] 
        E_bas    - effective activation energy of dissolution [J mol-1]
        T_pore   - temperature of the pore space [K]
        H_mol    - hydrogen ion molality in the pore space [mol kg-1]
        Hmod_mol - the modern H ion molality in pre space [mol kg-1]
        gamma    - empirical scaling parameter [dimensionless]

    Returns:
        F_diss - rate of seafloor basalt dissolution [mol eq yr-1]
    """

    Rg = 8.314 #universal gas constant [J mol-1 K-1]
    F_diss = k_diss*r_sr*exp(-E_bas/(Rg*T_pore))*2.88*10**-14*\
            10**(-gamma*(-log10(H_mol))) #see code from JKT

    return F_diss

def poreSpaceTemperature(T_D, Q, S_thick, K):
    """
    This function will calculate the temperature of the pore space. This is 
    based on equation S4.

    Inputs:
        T_D     - deep ocean temperature [K]
        Q       - pore space heat flow relative to modern Earth [dimensionless]
        S_thick - the thickness of the ocean sediment relative to modern Earth 
                  [dimensionless]
        K       - conductivity of the pore space sediments [m K-1]

    Returns:
    T_pore - the temperature of the pore space [K]
    """

    sed = S_thick*700 #modern Earth sediment thickness is ~700 m
    T_pore = T_D + Q*sed/K

    return T_pore

def globalOutgassing(Fmod_out, Q, m):
    """
    This function will calculate the outgassing flux (equation S9).

    Inputs:
        Fmod_out - the modern Earth's outgassing rate [mol C yr-1]
        Q        - pore space heat flow relative to modern Earth [dimensionless]
        m        - scaling parameter [dimensionless]

    Returns:
        F_out - the global outgassing flux [mol C yr-1]
    """

    F_out = Fmod_out*Q**m
    
    return F_out

def spreadingRate(Q, beta):
    """
    Calculates the spreading rate on the planet.

    Inputs:
        Q    - pore space heat flow relative to modern Earth [dimensionless]
        beta - scaling parameter [dimensionless]

    Returns:
        r_sr - the spreading rate relative to the modern Earth [dimensionless]
    """

    r_sr = Q**beta

    return r_sr

def equilibriumChemistry(T, alk, carb, s, alk_init, Ca_init):
    """
    Calculate the carbonate equilibrium and alkalinity. This can be used for
    either the atmosphere-ocean or the pore space. This function represents 
    equations S11-S18.

    Inputs:
        T        - the temperature of the system [K]
        alk      - the alkalinity of the system [mol eq]
        carb     - the carbon abundance in the system [mol]
        s        - correction factor for mass balance
        alk_init - initial alkalinity of the system [mol eq]
        Ca_init  - initial calcium ion concentration in the system [mol]

    Returns:
        omega - the saturation state of the system
        pCO2  - the partial pressure of CO2 [bar]
        pH    - pH of the system
    """

    #get the rate constants and Henry's constant

    [K1, K2, H_CO2] = equilibriumRateConstants(T)
    
    #use equation S15 to first calculate the H+ ion concentration
    roots = np.roots([alk/(K1*K2)*(1.0+s/H_CO2),
                      (alk-carb)/K2,
                      alk-2.0*carb])

    H_ion = np.max(roots) #just take the positive root
    pH = -log10(H_ion) #equation S16 (aka pH definition)

    CO3 = alk/(2.0+H_ion/K2) #S14 with S11
    HCO3 = alk - 2.0*CO3 #S11
    CO2_aq = H_ion*HCO3/K1 #S13
    pCO2 = CO2_aq/H_CO2 #S12
    Ca_ion = 0.5*(alk - alk_init) + Ca_init #S17
    K_sp = carbonateSolubility(T)
    omega = Ca_ion*CO3/K_sp # S18


    return [omega, pCO2, pH]


def carbonateSolubility(T):
    """
    Calculates carbonate solubility rate constant as a function of temperature.
    See Appendix A of JKT 2018 for further details (you'll need to look at
    their 2017 paper for these actual equations - but Appendix A tells you 
    that).

    Inputs:
        T - the temperature of the system [K]

    Returns:
        result - the solubility rate constant
    """
    bo = -0.77712
    b1 = 0.0028426
    b2 = 178.34
    co = -0.07711
    do = 0.0041249
    S = 35.0
    logK0=-171.9065-0.077993*T+2839.319/T+71.595*log10(T) 
    logK=logK0+(bo+b1*T+b2/T)*S**0.5+co*S+do*S**1.5

    result = 10.0**logK

    return result


def equilibriumRateConstants(T):
    """
    Calculates the carbon chemistry equilibrium constants as a function of 
    temperature following the method in Appendix A of JKT 2018 (you actually
    have to look at their 2017 paper for these equations).

    Inputs:
        T - the temperature of the system [K]

    Returns:
        K1    - the first apparent dissociation rate constant of carbonic acid
        K2    - the second apparent dissociation rate constant of carbonic acid
        H_CO2 - Henry's law constant for CO2
    """

    pK1=17.788 - .073104 *T - .0051087*35 + 1.1463*10**-4*T**2
    pK2=20.919 - .064209 *T - .011887*35 + 8.7313*10**-5*T**2
    H_CO2=exp(9345.17/T - 167.8108 + 23.3585 * log(T) + 
            (.023517 - 2.3656*10**-4*T+4.7036*10**-7*T**2)*35)
    
    K1 = 10.0**-pK1
    K2 = 10.0**-pK2

    return [K1,K2,H_CO2]


def deepOceanTemperature(Ts, gradient, min_temp=271.15):
    """
    Determine the deep ocean temperature based on the surface temperature. The
    intercept term is chosen so that gradient*Ts+intercept gives the correct
    surface temperature. In the case of the modern Earth, that would be the
    modern average surface temperature. This function corresponds to equation
    S20.

    Inputs:
        Ts        - surface temperature [K]
        gradient  - total temperature gradient in the ocean [dimensionless]
        min_temp  - the minimum allowable temperature at the bottom of the 
                    ocean. For an Earth-like planet below 271.15 K (the default
                    value) the ocean would freeze.

    Returns:
        Td - the temperature at the bottom of the ocean [K]
    """

    # intercept chosen to reproduce initial (modern) temperature
    intercept = 274.037 - gradient*Ts 
    Td = np.max([np.min([gradient*Ts+intercept, Ts]), min_temp])

    return Td


def surfaceTempFromClimateModel(pCO2, flux):
    """
    This function will return the surface temperature of an Earth-like planet for
    the given partial pressure of CO2 and incident flux (normalized to modern
    Earth's). The function is defined for CO2 levels between >1.0E-7 and <10 
    bar. The flux definitions are from 1.05 to 0.31 (the HZ for a Sun-like
    star). The fit is to a 4th order polynomial over CO2 and flux.

    Inputs:
        pCO2 - the CO2 partial pressure of the atmosphere [bar]
        flux - the incident flux normalised to the modern Earths (i.e. divided
               by ~1360 [W m-2])

    Returns:
        the surface temperature of the planet [K]
    """
    #the fit was done in log space for CO2
    x = np.log(pCO2)
    y = flux

    coeffs = np.array([4.8092693271e+00, -2.2201836059e+02, -6.8437057004e+01, 
        -6.7369814833e+00, -2.0576569974e-01, 1.4144615786e+03, 
        4.4638645525e+02, 4.4412679359e+01, 1.3641352778e+00, -2.9643244170e+03, 
        -9.7844390774e+02, -9.8858815404e+01, -3.0586461777e+00, 
        2.6547903068e+03, 9.0749599550e+02, 9.2870700889e+01, 2.8915352308e+00, 
        -8.6843290311e+02, -3.0464088878e+02, -3.1476199768e+01, 
        -9.8478712084e-01, 1.0454688611e+03, -1.4964888001e+03, 
        1.0637917601e+03, -2.8114373919e+02])

    p4_in = np.array([1, x, x**2, x**3, x**4, x*y, x**2*y, x**3*y, x**4*y, 
        x*y**2, x**2*y**2, x**3*y**2, x**4*y**2, x*y**3, x**2*y**3, x**3*y**3, 
        x**4*y**3, x*y**4, x**2*y**4, x**3*y**4, x**4*y**4, y, y**2, y**3, 
        y**4])


    return np.sum(p4_in*coeffs)



#UNCOMMENT the below lines to run the model. Use the ModelInputs() class to 
#alter the inputs for the model.
#inpts = ModelInputs()
#inpts.lum = 0.9
#[pCO2_out, pH_out, Ts_out], status = runWeatheringModel(inputs=inpts,
#        guess=None, chem_tol=1.0E-5, chem_max_iter=10)
#
#print("Flux: %0.2f"%(inpts.lum))
#print("pH=%2.3e"%(pH_out))
#print("pCO2=%2.3e"%(pCO2_out))
#print("Ts=%2.3e"%(Ts_out))
