"""
##############################################################################
##############################################################################
TMS model is built by Bin Xu for the line 500-836, 1390-1631
##############################################################################
##############################################################################
"""

"""
##############################################################################
##############################################################################
Propulsion Model is the Pythonic copy of NREL's FASTSim
(Future Automotive Systems Technology Simulator)
Input Arguments
1) cyc: dictionary defining drive cycle to be simulated
        cyc['cycSecs']: drive cycle time in seconds (begins at zero)
        cyc['cycMps']: desired vehicle speed in meters per second
        cyc['cycGrade']: road grade
        cyc['cycRoadType']: Functional Class of GPS data
2) veh: dictionary defining vehicle parameters for simulation
Output Arguments
A dictionary containing scalar values summarizing simulation
(mpg, Wh/mi, avg engine power, etc) and/or time-series results for component
power and efficiency. Currently these values are user clarified in the code by assigning values to the 'output' dictionary.
    List of Abbreviations
    cur = current time step
    prev = previous time step

    cyc = drive cycle
    secs = seconds
    mps = meters per second
    mph = miles per hour
    kw = kilowatts, unit of power
    kwh = kilowatt-hour, unit of energy
    kg = kilograms, unit of mass
    max = maximum
    min = minimum
    avg = average
    fs = fuel storage (eg. gasoline/diesel tank, pressurized hydrogen tank)
    fc = fuel converter (eg. internal combustion engine, fuel cell)
    mc = electric motor/generator and controller
    ess = energy storage system (eg. high voltage traction battery)

    chg = charging of a component
    dis = discharging of a component
    lim = limit of a component
    regen = associated with regenerative braking
    des = desired value
    ach = achieved value
    in = component input
    out = component output
##############################################################################
##############################################################################
"""

### Import necessary python modules
import numpy as np
import pandas as pd
import warnings
import csv
import CoolProp.CoolProp as CP
import math
import time
t0=time.time()
from scipy import interpolate
warnings.simplefilter('ignore')

def get_standard_cycle(cycle_name):
    csv_path = 'cycles/'+cycle_name+'.csv'
    data = dict()
    dkeys=[]
    with open(csv_path) as csvfile:
        csv_reader = csv.reader(csvfile)
        for row in csv_reader:
            if len(data)==0: # initialize all elements in dictionary based on header
                for ii in range(len(row)):
                    data[row[ii]] = []
                    dkeys.append( row[ii] )
            else: # append values
                for ii in range(len(row)):
                    try:
                        data[dkeys[ii]].append( float(row[ii]) )
                    except:
                        data[dkeys[ii]].append( np.nan )
    for ii in range(len(dkeys)):
        data[dkeys[ii]] = np.array(data[dkeys[ii]])
    return data

def get_veh(vnum):
    with open('docs/FASTSim_py_veh_db.csv','r') as csvfile:
#     with open('//Users//Mingjue//EV_EMS//docs//FASTSim_py_veh_db.csv','r') as csvfile:

        reader = csv.reader(csvfile)
        vd = dict()
        data = dict()
        z=0

        for i in reader:
            data[z]=i
            z=z+1

        variables = data[0]
        del data[0] # deletes the first list, which corresponds to the header w/ variable names
        vd=data

        ### selects specified vnum from vd
        veh = dict()
        variables = ['selection','name', 'vehPtType', 'dragCoef', 'frontalAreaM2', 'gliderKg', 'vehCgM', 'driveAxleWeightFrac', 'wheelBaseM', 'cargoKg', 'vehOverrideKg', 'maxFuelStorKw', 'fuelStorSecsToPeakPwr', 'fuelStorKwh', 'fuelStorKwhPerKg', 'maxFuelConvKw', 'fcEffType', 'fcAbsEffImpr', 'fuelConvSecsToPeakPwr', 'fuelConvBaseKg', 'fuelConvKwPerKg', 'maxMotorKw', 'motorPeakEff', 'motorSecsToPeakPwr', 'mcPeKgPerKw', 'mcPeBaseKg', 'maxEssKw', 'maxEssKwh', 'essKgPerKwh', 'essBaseKg', 'essRoundTripEff', 'essLifeCoefA', 'essLifeCoefB', 'wheelInertiaKgM2', 'numWheels', 'wheelRrCoef', 'wheelRadiusM', 'wheelCoefOfFric', 'minSoc', 'maxSoc', 'essDischgToFcMaxEffPerc', 'essChgToFcMaxEffPerc', 'maxAccelBufferMph', 'maxAccelBufferPercOfUseableSoc', 'percHighAccBuf', 'mphFcOn', 'kwDemandFcOn', 'altEff', 'chgEff', 'auxKw', 'forceAuxOnFC', 'transKg', 'transEff', 'compMassMultiplier', 'essToFuelOkError', 'maxRegen', 'valUddsMpgge', 'valHwyMpgge', 'valCombMpgge', 'valUddsKwhPerMile', 'valHwyKwhPerMile', 'valCombKwhPerMile', 'valCdRangeMi', 'valConst65MphKwhPerMile', 'valConst60MphKwhPerMile', 'valConst55MphKwhPerMile', 'valConst45MphKwhPerMile', 'valUnadjUddsKwhPerMile', 'valUnadjHwyKwhPerMile', 'val0To60Mph', 'valEssLifeMiles', 'valRangeMiles', 'valVehBaseCost', 'valMsrp', 'minFcTimeOn', 'idleFcKw','fuelKwhPerKg']
        if vnum in vd:
            for i in range(len(variables)):
                vd[vnum][i]=str(vd[vnum][i])
                if vd[vnum][i].find('%') != -1:
                    vd[vnum][i]=vd[vnum][i].replace('%','')
                    vd[vnum][i]=float(vd[vnum][i])
                    vd[vnum][i]=vd[vnum][i]/100.0
                elif vd[vnum][i].find('TRUE') != -1 or vd[vnum][i].find('True') != -1 or vd[vnum][i].find('true') != -1:
                    vd[vnum][i]=1
                elif vd[vnum][i].find('FALSE') != -1 or vd[vnum][i].find('False') != -1 or vd[vnum][i].find('false') != -1:
                    vd[vnum][i]=1
                else:
                    try:
                        vd[vnum][i]=float(vd[vnum][i])
                    except:
                        pass
                veh[variables[i]]=vd[vnum][i]

    ######################################################################
    ### Append additional parameters to veh structure from calculation ###
    ######################################################################

    ### Build roadway power lookup table
    veh['MaxRoadwayChgKw_Roadway'] = range(6)
    veh['MaxRoadwayChgKw'] = [0]*len(veh['MaxRoadwayChgKw_Roadway'])
    veh['chargingOn'] = 0

     # Checking if a vehicle has any hybrid components
    if veh['maxEssKwh']==0 or veh['maxEssKw']==0 or veh['maxMotorKw']==0:
        veh['noElecSys'] = 'TRUE'

    else:
        veh['noElecSys'] = 'FALSE'

    # Checking if aux loads go through an alternator
    if veh['noElecSys']=='TRUE' or veh['maxMotorKw']<=veh['auxKw'] or veh['forceAuxOnFC']=='TRUE':
        veh['noElecAux'] = 'TRUE'

    else:
        veh['noElecAux'] = 'FALSE'

    veh['vehTypeSelection'] = np.copy( veh['vehPtType'] ) # Copying vehPtType to additional key

    ### Defining Fuel Converter efficiency curve as lookup table for %power_in vs power_out
    ### see "FC Model" tab in FASTSim for Excel

    if veh['maxFuelConvKw']>0:


        # Discrete power out percentages for assigning FC efficiencies
        fcPwrOutPerc = np.array([0, 0.005, 0.015, 0.04, 0.06, 0.10, 0.14, 0.20, 0.40, 0.60, 0.80, 1.00])

        # Efficiencies at different power out percentages by FC type
        eff_si = np.array([0.10, 0.12, 0.16, 0.22, 0.28, 0.33, 0.35, 0.36, 0.35, 0.34, 0.32, 0.30])
        eff_atk = np.array([0.10, 0.12, 0.28, 0.35, 0.375, 0.39, 0.40, 0.40, 0.38, 0.37, 0.36, 0.35])
        eff_diesel = np.array([0.10, 0.14, 0.20, 0.26, 0.32, 0.39, 0.41, 0.42, 0.41, 0.38, 0.36, 0.34])
        eff_fuel_cell = np.array([0.10, 0.30, 0.36, 0.45, 0.50, 0.56, 0.58, 0.60, 0.58, 0.57, 0.55, 0.54])
        eff_hd_diesel = np.array([0.10, 0.14, 0.20, 0.26, 0.32, 0.39, 0.41, 0.42, 0.41, 0.38, 0.36, 0.34])
        veh['fcPwrOutPerc'] =fcPwrOutPerc 
        veh['eff_si']=eff_si

        if veh['fcEffType']==1:
            eff = np.copy( eff_si ) + veh['fcAbsEffImpr']

        elif veh['fcEffType']==2:
            eff = np.copy( eff_atk ) + veh['fcAbsEffImpr']

        elif veh['fcEffType']==3:
            eff = np.copy( eff_diesel ) + veh['fcAbsEffImpr']

        elif veh['fcEffType']==4:
            eff = np.copy( eff_fuel_cell ) + veh['fcAbsEffImpr']

        elif veh['fcEffType']==5:
            eff = np.copy( eff_hd_diesel ) + veh['fcAbsEffImpr']

        inputKwOutArray = fcPwrOutPerc * veh['maxFuelConvKw']
        fcPercOutArray = np.r_[np.arange(0,3.0,0.1),np.arange(3.0,7.0,0.5),np.arange(7.0,60.0,1.0),np.arange(60.0,105.0,5.0)] / 100
        fcKwOutArray = veh['maxFuelConvKw'] * fcPercOutArray
        fcEffArray = np.array([0.0]*len(fcPercOutArray))

        for j in range(0,len(fcPercOutArray)-1):

            low_index = np.argmax(inputKwOutArray>=fcKwOutArray[j])
            fcinterp_x_1 = inputKwOutArray[low_index-1]
            fcinterp_x_2 = inputKwOutArray[low_index]
            fcinterp_y_1 = eff[low_index-1]
            fcinterp_y_2 = eff[low_index]
            fcEffArray[j] = (fcKwOutArray[j] - fcinterp_x_1)/(fcinterp_x_2 - fcinterp_x_1)*(fcinterp_y_2 - fcinterp_y_1) + fcinterp_y_1

        fcEffArray[-1] = eff[-1]
        veh['fcEffArray'] = np.copy(fcEffArray)
        veh['fcKwOutArray'] = np.copy(fcKwOutArray)
        veh['maxFcEffKw'] = np.copy(veh['fcKwOutArray'][np.argmax(fcEffArray)])
        veh['fcMaxOutkW'] = np.copy(max(inputKwOutArray))
        veh['minFcTimeOn'] = 30

    else:
        veh['fcKwOutArray'] = np.array([0]*101)
        veh['maxFcEffKw'] = 0
        veh['fcMaxOutkW'] = 0
        veh['minFcTimeOn'] = 30

    ### Defining MC efficiency curve as lookup table for %power_in vs power_out
    ### see "Motor" tab in FASTSim for Excel
    if veh['maxMotorKw']>0:

        maxMotorKw = veh['maxMotorKw']

        mcPwrOutPerc = np.array([0.00, 0.02, 0.04, 0.06, 0.08,	0.10,	0.20,	0.40,	0.60,	0.80,	1.00])
        large_baseline_eff = np.array([0.88, 0.90,	0.93,	0.95,	0.97,	0.98,	0.98,	0.98,	0.98,	0.98,	0.97])
#         large_baseline_eff = np.array([0.83, 0.85,	0.87,	0.89,	0.90,	0.91,	0.93,	0.94,	0.94,	0.93,	0.92])
        small_baseline_eff = np.array([0.12,	0.16,	 0.21, 0.29, 0.35, 0.42, 0.75, 0.92, 0.93,	0.93,	0.92])

        modern_max = 0.95
#         modern_max = 0.95
        modern_diff = modern_max - max(large_baseline_eff)

        large_baseline_eff_adj = large_baseline_eff + modern_diff

        mcKwAdjPerc = max(0.0,min((maxMotorKw - 7.5)/(75.0-7.5),1.0))
        mcEffArray = np.array([0.0]*len(mcPwrOutPerc))

        for k in range(0,len(mcPwrOutPerc)):
            mcEffArray[k] = mcKwAdjPerc*large_baseline_eff_adj[k] + (1-mcKwAdjPerc)*(small_baseline_eff[k])

        mcInputKwOutArray = mcPwrOutPerc * maxMotorKw

        mcPercOutArray = np.linspace(0,1,101)
        mcKwOutArray = np.linspace(0,1,101) * maxMotorKw

        mcFullEffArray = np.array([0.0]*len(mcPercOutArray))

        for m in range(1,len(mcPercOutArray)-1):
            low_index = np.argmax(mcInputKwOutArray>=mcKwOutArray[m])

            fcinterp_x_1 = mcInputKwOutArray[low_index-1]
            fcinterp_x_2 = mcInputKwOutArray[low_index]
            fcinterp_y_1 = mcEffArray[low_index-1]
            fcinterp_y_2 = mcEffArray[low_index]

            mcFullEffArray[m] = (mcKwOutArray[m] - fcinterp_x_1)/(fcinterp_x_2 - fcinterp_x_1)*(fcinterp_y_2 - fcinterp_y_1) + fcinterp_y_1

        mcFullEffArray[0] = 0
        mcFullEffArray[-1] = mcEffArray[-1]

        mcKwInArray = mcKwOutArray / mcFullEffArray
        mcKwInArray[0] = 0

        veh['mcKwInArray'] = np.copy(mcKwInArray)
        veh['mcKwOutArray'] = np.copy(mcKwOutArray)
        veh['mcMaxElecInKw'] = np.copy(max(mcKwInArray))
        veh['mcFullEffArray'] = np.copy(mcFullEffArray)
        veh['mcEffArray'] = np.copy(mcEffArray)

    else:
        veh['mcKwInArray'] = np.array([0.0] * 101)
        veh['mcKwOutArray'] = np.array([0.0]* 101)
        veh['mcMaxElecInKw'] = 0

    veh['mcMaxElecInKw'] = max(veh['mcKwInArray'])

    ### Specify shape of mc regen efficiency curve
    ### see "Regen" tab in FASTSim for Excel
    veh['regenA'] = 500.0
    veh['regenB'] = 0.99

    ### Calculate total vehicle mass
    if veh['vehOverrideKg']==0 or veh['vehOverrideKg']=="":
        if veh['maxEssKwh']==0 or veh['maxEssKw']==0:
            ess_mass_kg = 0.0
        else:
            ess_mass_kg = ((veh['maxEssKwh']*veh['essKgPerKwh'])+veh['essBaseKg'])*veh['compMassMultiplier']
        if veh['maxMotorKw']==0:
            mc_mass_kg = 0.0
        else:
            mc_mass_kg = (veh['mcPeBaseKg']+(veh['mcPeKgPerKw']*veh['maxMotorKw']))*veh['compMassMultiplier']
        if veh['maxFuelConvKw']==0:
            fc_mass_kg = 0.0
        else:
            fc_mass_kg = (((1/veh['fuelConvKwPerKg'])*veh['maxFuelConvKw']+veh['fuelConvBaseKg']))*veh['compMassMultiplier']
        if veh['maxFuelStorKw']==0:
            fs_mass_kg = 0.0
        else:
            fs_mass_kg = ((1/veh['fuelStorKwhPerKg'])*veh['fuelStorKwh'])*veh['compMassMultiplier']
        veh['vehKg'] = veh['cargoKg'] + veh['gliderKg'] + veh['transKg']*veh['compMassMultiplier'] + ess_mass_kg + mc_mass_kg + fc_mass_kg + fs_mass_kg
#         print('L281, vehKg:',veh['vehKg'])
    else:
        veh['vehKg'] = np.copy( veh['vehOverrideKg'] )

    return veh

def sim_drive( cyc , veh ):

    if veh['vehPtType']==1:

        # If no EV / Hybrid components, no SOC considerations.

        initSoc = 0.0
        output = sim_drive_sub( cyc , veh , initSoc )

    elif veh['vehPtType']==2:

        #####################################
        ### Charge Balancing Vehicle SOC ###
        #####################################

        # Charge balancing SOC for PHEV vehicle types. Iterating initsoc and comparing to final SOC.
        # Iterating until tolerance met or 30 attempts made.
        ##########################################################################
        ##########################################################################
        ##########################################################################
        initSoc = 0.6
#         initSoc = (veh['maxSoc']+veh['minSoc'])/2.0        
#         print('004/',initSoc)
        ess2fuelKwh = 1.0
        sim_count = 0
        ##########################################################################
        ##########################################################################
        ##########################################################################   
        cyc_num = 'one'
#         cyc_num = 'mul'
        if cyc_num=='one':
            # one cycle
            sim_count += 1
            output = sim_drive_sub( cyc , veh , initSoc )
            ess2fuelKwh = abs( output['ess2fuelKwh'] )
        else:
        ##########################################################################
        ##########################################################################
        ##########################################################################
            # multiple cycle
            while ess2fuelKwh>veh['essToFuelOkError'] and sim_count<30:
                sim_count += 1
                output = sim_drive_sub( cyc , veh , initSoc )
                ess2fuelKwh = abs( output['ess2fuelKwh'] )
                initSoc = min(1.0,max(0.0,output['final_soc']))
        np.copy( veh['maxSoc'] )
        output = sim_drive_sub( cyc , veh , initSoc )        

    elif veh['vehPtType']==3 or veh['vehPtType']==4:

        # If EV, initializing initial SOC to maximum SOC.

        initSoc = np.copy( veh['maxSoc'] )
        output = sim_drive_sub( cyc , veh , initSoc )

    return output

def sim_drive_sub( cyc , veh , initSoc):

    # sim_drive_sub receives second-by-second cycle information,
    # vehicle properties, and an initial state of charge and performs
    # a backward facing powertrain simulation. The function returns an 
    # output dictionary starting at approximately line 1030. Powertrain
    # variables of interest (summary or time-series) can be added to the 
    # output dictionary for reference.
    
    ############################
    ###   Define Constants   ###
    ############################
    if veh['maxFuelConvKw']>0:
        fcPwrOutPerc  = veh['fcPwrOutPerc'] 
        eff_si = veh['eff_si']
        f = interpolate.interp1d(fcPwrOutPerc,eff_si)

    airDensityKgPerM3 = 1.2 # Sea level air density at approximately 20C
    gravityMPerSec2 = 9.81
    mphPerMps = 2.2369
    kWhPerGGE = 33.7
    metersPerMile = 1609.00
    maxTracMps2 = ((((veh['wheelCoefOfFric']*veh['driveAxleWeightFrac']*veh['vehKg']*gravityMPerSec2)/(1+((veh['vehCgM']*veh['wheelCoefOfFric'])/veh['wheelBaseM']))))/(veh['vehKg']*gravityMPerSec2))*gravityMPerSec2
    maxRegenKwh = 0.5*veh['vehKg']*(27**2)/(3600*1000)

    #############################
    ### Initialize Variables  ###
    #############################

    ### Drive Cycle
    cycSecs = np.copy( cyc['cycSecs'] )
    cycMps = np.copy( cyc['cycMps'] )
    cycGrade = np.copy( cyc['cycGrade'] )
    cycRoadType = np.copy( cyc['cycRoadType'] )
    cycMph = [x * mphPerMps for x in cyc['cycMps']]
    secs = np.insert(np.diff(cycSecs),0,0)

    ### Component Limits
    curMaxFsKwOut = [0]*len(cycSecs)
    fcTransLimKw = [0]*len(cycSecs)
    fcFsLimKw = [0]*len(cycSecs)
    fcMaxKwIn = [0]*len(cycSecs)
    curMaxFcKwOut = [0]*len(cycSecs)
    essCapLimDischgKw = [0]*len(cycSecs)
    curMaxEssKwOut = [0]*len(cycSecs)
    curMaxAvailElecKw = [0]*len(cycSecs)
    essCapLimChgKw = [0]*len(cycSecs)
    curMaxEssChgKw = [0]*len(cycSecs)
    curMaxRoadwayChgKw = np.interp( cycRoadType, veh['MaxRoadwayChgKw_Roadway'], veh['MaxRoadwayChgKw'] )
    curMaxElecKw = [0]*len(cycSecs)
    mcElecInLimKw = [0]*len(cycSecs)
    mcTransiLimKw = [0]*len(cycSecs)
    curMaxMcKwOut = [0]*len(cycSecs)
    essLimMcRegenPercKw = [0]*len(cycSecs)
    essLimMcRegenKw = [0]*len(cycSecs)
    curMaxMechMcKwIn = [0]*len(cycSecs)
    curMaxTransKwOut = [0]*len(cycSecs)

    ### Drive Train
    cycDragKw = [0]*len(cycSecs)
    cycAccelKw = [0]*len(cycSecs)
    cycAscentKw = [0]*len(cycSecs)
    cycTracKwReq = [0]*len(cycSecs)
    curMaxTracKw = [0]*len(cycSecs)
    spareTracKw = [0]*len(cycSecs)
    cycRrKw = [0]*len(cycSecs)
    cycWheelRadPerSec = [0]*len(cycSecs)
    cycTireInertiaKw = [0]*len(cycSecs)
    cycWheelKwReq = [0]*len(cycSecs)
    regenContrLimKwPerc = [0]*len(cycSecs)
    cycRegenBrakeKw = [0]*len(cycSecs)
    cycFricBrakeKw = [0]*len(cycSecs)
    cycTransKwOutReq = [0]*len(cycSecs)
    cycMet = [0]*len(cycSecs)
    transKwOutAch = [0]*len(cycSecs)
    transKwInAch = [0]*len(cycSecs)
    curSocTarget = [0]*len(cycSecs)
    minMcKw2HelpFc = [0]*len(cycSecs)
    mcMechKwOutAch = [0]*len(cycSecs)
    mcElecKwInAch = [0]*len(cycSecs)
    auxInKw = [0]*len(cycSecs)

    #roadwayMaxEssChg = [0]*len(cycSecs)
    roadwayChgKwOutAch = [0]*len(cycSecs)
    minEssKw2HelpFc = [0]*len(cycSecs)
    essKwOutAch = [0]*len(cycSecs)
    fcKwOutAch = [0]*len(cycSecs)
    fcKwOutAch_pct = [0]*len(cycSecs)
    fcKwInAch = [0]*len(cycSecs)
    fsKwOutAch = [0]*len(cycSecs)
    fsKwhOutAch = [0]*len(cycSecs)
    essCurKwh = [0]*len(cycSecs)
    soc = [0]*len(cycSecs)

    # Vehicle Attributes, Control Variables
    regenBufferSoc = [0]*len(cycSecs)
    essRegenBufferDischgKw = [0]*len(cycSecs)
    maxEssRegenBufferChgKw = [0]*len(cycSecs)
    essAccelBufferChgKw = [0]*len(cycSecs)
    accelBufferSoc = [0]*len(cycSecs)
    maxEssAccelBufferDischgKw = [0]*len(cycSecs)
    essAccelRegenDischgKw = [0]*len(cycSecs)
    mcElectInKwForMaxFcEff = [0]*len(cycSecs)
    electKwReq4AE = [0]*len(cycSecs)
    canPowerAllElectrically = [0]*len(cycSecs)
    desiredEssKwOutForAE = [0]*len(cycSecs)
    essAEKwOut = [0]*len(cycSecs)
    erAEKwOut = [0]*len(cycSecs)
    essDesiredKw4FcEff = [0]*len(cycSecs)
    essKwIfFcIsReq = [0]*len(cycSecs)
    curMaxMcElecKwIn = [0]*len(cycSecs)
    fcKwGapFrEff = [0]*len(cycSecs)
    erKwIfFcIsReq = [0]*len(cycSecs)
    mcElecKwInIfFcIsReq = [0]*len(cycSecs)
    mcKwIfFcIsReq = [0]*len(cycSecs)
    fcForcedOn = np.full(len(cycSecs),False)
    fcForcedState = [0]*len(cycSecs)
    mcMechKw4ForcedFc = [0]*len(cycSecs)
    fcTimeOn = [0]*len(cycSecs)
    prevfcTimeOn = [0]*len(cycSecs)

    ### Additional Variables
    mpsAch = [0]*len(cycSecs)
    mphAch = [0]*len(cycSecs)
    distMeters = [0]*len(cycSecs)
    distMiles = [0]*len(cycSecs)
    highAccFcOnTag = [0]*len(cycSecs)
    reachedBuff = [0]*len(cycSecs)
    maxTracMps = [0]*len(cycSecs)
    addKwh = [0]*len(cycSecs)
    dodCycs = [0]*len(cycSecs)
    essPercDeadArray = [0]*len(cycSecs)
    dragKw = [0]*len(cycSecs)
    essLossKw = [0]*len(cycSecs)
    accelKw = [0]*len(cycSecs)
    ascentKw = [0]*len(cycSecs)
    rrKw = [0]*len(cycSecs)
    motor_index_debug = [0]*len(cycSecs)
    debug_flag = [0]*len(cycSecs)

    ############################
    ###  Assign First Value  ###
    ############################

    ### Drive Train
    cycMet[0] = 1
    curSocTarget[0] = veh['maxSoc']
    essCurKwh[0] = initSoc*veh['maxEssKwh']
    soc[0] = initSoc
    
    ############################
    ###         TMS          ###
    ############################
    T_init          = veh['T_amb_ref']
    T_cab_REF       = veh['T_cab_ref']
    T_ess_ref       = veh['T_ess_ref']
    cycL=len(cycSecs)
    Cy=veh['Cy']
    car=veh['car']
    AC_switch  = veh['AC_switch']
    fluid      = 'R134a'

    ### constant
    dt         = 0.25
    steps      = int(len(cycSecs)/dt+1)

    ## essT control
    kp_essT=5;ki_essT=1/1000
    p_term_essT      = np.zeros(steps)
    i_term_essT      = np.zeros(steps)
    e_essT           = np.zeros(steps)

    ## chill rfg superheat control 
    kp_chillTrfg=0.1;ki_chillTrfg=0.1
    p_term_chillTrfg = np.zeros(steps)
    i_term_chillTrfg = np.zeros(steps)
    e_chillTrfg      = np.zeros(steps)

    ## evap rfg superheat control 
    kp_evapTrfg=0.01;ki_evapTrfg=0.01
    p_term_evapTrfg  = np.zeros(steps)
    i_term_evapTrfg  = np.zeros(steps)
    e_evapTrfg       = np.zeros(steps)

    ## cond rfg subcool control 
    kp_condTrfg=0.1;ki_condTrfg=0.001
    p_term_condTrfg  = np.zeros(steps)
    i_term_condTrfg  = np.zeros(steps)
    e_condTrfg       = np.zeros(steps)

    ## Tcab control
    kp_Tcab=0.01;ki_Tcab=0.001
    p_term_Tcab      = np.zeros(steps)
    i_term_Tcab      = np.zeros(steps)
    e_Tcab           = np.zeros(steps)

    # constant
    a          = np.array([1,1])
    if car=='model3':
        essm       = 350  # [kg]
    if car=='modelS':
        essm       = 400  # [kg]
    essc       = 1000 # [J/K/kg]
    essmcool   = 2    # [kg]
    essccool   = 3323 # [J/K/kg] ref: 2016-ATE-Comparison of different cooling methods for lithium ion battery cells
    chillmcool = 0.5    # [kg]
    chillccool = 3323 # [J/K/kg]
    chillmw    = 1    # [kg]
    chillcw    = 871  # [J/K/kg] aluminum
    chillmrfg  = 0.5  # [kg]
    chillcrfg  = 3000 # [J/K/kg]
    radmcool   = 1    # [kg]
    radccool   = 3000 # [J/K/kg]
    radmw      = 2    # [kg]
    radcw      = 871  # [J/K/kg]
    radcair    = 1006  # [J/K/kg]
    evapmrfg   = 1  # [kg]
    evapcair   = 1006  # [J/K/kg]
    condmrfg   = 2  # [kg]
    resmrfg    = 1  # [kg]
    evapmw     = 2    # [kg]
    evapcw     = 871  # [J/K/kg]
    condmw     = 2    # [kg]
    condcw     = 871  # [J/K/kg]
    cabcair    = 1006  # [J/K/kg]
    cabmair    = 3  # [kg]
    condcair   = 1006  # [J/K/kg]
    condmair   = 0.1  # [kg]
    rfgVcomp   = 50  #     rfgVcomp=32
    condeff    = 0.9
    evapeff    = 0.9
    radeff     = 0.9
    chilleff   = 0.9
    essHXeff   = 0.9
    R               = 0.15   # 2020-reseaerch on control strategy for a battery TMS for EV based on secondary loop cooling
    if car=='model3':
        N_cell          = 444*11
    if car=='modelS':
        N_cell          = 444*14
    V_nominal       = 3.8
    
    #cab const
    R_sun      = 800  # [W/m^2]
    A_sun      = 2    # [m^2]
    M_driver   = 85   # [W/m^2]
    A_driver   = 1.74 # [m^2]
    K_amb      = 30   # [W/m^2]

    ## compressor eff
    # eff_is
    a1 =      1.1171;
    a2 =    -0.07598;
    a3 =    -0.01733;
    a4 =     0.01593;
    a5 =   -0.001453;
    a6 =   -0.001061;
    # eff_vol
    b1 =       1.644;
    b2 =        -0.2;
    b3 =     -0.2246;
    b4 =     0.06626;
    b5 =     0.01856;
    b6 =   -0.005932; 
    eff_mech_comp=0.9
    eff_ele_comp =0.88

    ## cool loop    
    coolNpump     = 1000  # [rpm]
    coolpdownpump = 3e5   # [par]
    coolpuppump   = 1.1e5 # [par]
    cooldens      = 1000  # [kg/m^3]
    cooliseffpump = 0.25       
    chillhcoolw   = 140   # [w/m^2K]  2006-battey thermal management system design modeling
    chillAcoolw   = 0.2   # [m^2] 
    chillhrfgw    = 1000   # [w/m^2K]
    chillArfgw    = 0.2   # [m^2] 
    radhcoolw     = 140   # [w/m^2K]
    radAcoolw     = 0.5     # [m^2] 
    radhairw      = 60   # [w/m^2K] 2006-battey thermal management system design modeling
    radAairw      = 0.5     # [m^2]     
    evaphrfgw     = 3000  # [w/m^2K] 2020 - Structure design and effect analysis on refrigerant cooling enhancement of battery thermal management system for electric vehicles
    evapArfgw     = 2     # [m^2]     
    evaphairw     = 60   # [w/m^2K]
    evapAairw     = 2     # [m^2]   
    cabhairw      = 35   # [w/m^2K] # 2021 ECM A novel electric vehicle thermal management system based on cooling and heating of batteries by refrigerant
    cabAairw      = 3     # [m^2]   
    condhrfgw     = 1500  # [w/m^2K] # 2018 - Battery thermal management system for electric vehicle using heat pipes
    condArfgw     = 3     # [m^2]
    condhairw     = 60   # [w/m^2K]
    condAairw     = 3     # [m^2]

    ## array
    essQheat      = np.zeros(steps)
    essT          = np.zeros(steps)
    essTcool      = np.zeros(steps)
    essTcoolin    = np.zeros(steps)
    essTcoolout   = np.zeros(steps)
    essmdotcool   = np.zeros(steps)
    coolmdot      = np.zeros(steps)
    coolPpump     = np.zeros(steps)
    rfgPpump      = np.zeros(steps)
    rfgPcomp      = np.zeros(steps)
    chillmdotrfg  = np.zeros(steps)
    chillTcool    = np.zeros(steps)
    chillTcoolin  = np.zeros(steps)
    chillTcoolout = np.zeros(steps)
    chillTw       = np.zeros(steps)
    chillTrfg     = np.zeros(steps)
    resTrfgout    = np.zeros(steps)
    chillTrfgout  = np.zeros(steps)
    radmdotcool   = np.zeros(steps)
    radmdotair    = np.zeros(steps)
    radTcool      = np.zeros(steps)
    radTcoolin    = np.zeros(steps)
    radTcoolout   = np.zeros(steps)
    radTw         = np.zeros(steps)
    radTair       = np.zeros(steps)
    radTairin     = np.zeros(steps)
    radTairout    = np.zeros(steps)
    plowrfg       = np.zeros(steps)
    phighrfg      = np.zeros(steps)
    evapmdotrfg   = np.zeros(steps)
    evapmdotair   = np.zeros(steps)
    evaphrfg      = np.zeros(steps)
    evaphrfgin    = np.zeros(steps)
    evaphrfgout   = np.zeros(steps)
    chillhrfg     = np.zeros(steps)
    chillhrfgout  = np.zeros(steps)
    mixhrfg       = np.zeros(steps)
    condhrfg      = np.zeros(steps)
    condhrfgin    = np.zeros(steps)
    condhrfgout   = np.zeros(steps)
    evapTrfg      = np.zeros(steps)
    evapTrfgin    = np.zeros(steps)
    evapTrfgout   = np.zeros(steps)
    evapTw        = np.zeros(steps)
    evapTair      = np.zeros(steps)
    evapTairout   = np.zeros(steps)
    evapPfanair   = np.zeros(steps)
    condTrfg      = np.zeros(steps)
    condTrfgin    = np.zeros(steps)
    condTrfgout   = np.zeros(steps)
    condTw        = np.zeros(steps)
    condTair      = np.zeros(steps)
    condTairout   = np.zeros(steps)
    condmdotair   = np.zeros(steps)
    resTrfg       = np.zeros(steps)
    reshrfg       = np.zeros(steps)
    cabTair       = np.zeros(steps)
    cabTairout    = np.zeros(steps)
    cabmdotair    = np.zeros(steps)
    cabTw         = np.zeros(steps)
    rfgPcomp      = np.zeros(steps)
    rfgNcomp      = np.zeros(steps)
    sumPTMS       = np.zeros(steps)

    # initial condition load for TMS
    data_dir='data/'
    if T_ess_ref==25:
        if car=='model3':
            if Cy=='udds':
                if AC_switch=='on':
                    init        = pd.read_csv(data_dir+'EV_udds_model3_init.csv')
                if AC_switch=='off':
                    init        = pd.read_csv(data_dir+'EV_udds_off_model3_init.csv')
            if Cy=='hwfet':
                if AC_switch=='on':
                    init        = pd.read_csv(data_dir+'EV_hwfet_model3_init.csv')
                if AC_switch=='off':
                    init        = pd.read_csv(data_dir+'EV_hwfet_off_model3_init.csv')    
            if Cy=='wltp':
                init        = pd.read_csv(data_dir+'EV_wltp_model3_init.csv')
            if Cy=='highspeed_calm' or Cy=='highspeed_aggressive':
                init        = pd.read_csv(data_dir+'EV_hwfet_model3_init.csv')
        if car=='modelS':
            if Cy=='udds':
                if AC_switch=='on':
                    init        = pd.read_csv(data_dir+'EV_udds_modelS_init.csv')
                if AC_switch=='off':
                    init        = pd.read_csv(data_dir+'EV_udds_off_modelS_init.csv')
            if Cy=='hwfet':
                if AC_switch=='on':
                    init        = pd.read_csv(data_dir+'EV_hwfet_modelS_init.csv')
                if AC_switch=='off':
                    init        = pd.read_csv(data_dir+'EV_hwfet_off_modelS_init.csv')    
            if Cy=='highspeed_calm' or Cy=='highspeed_aggressive':
                init        = pd.read_csv(data_dir+'EV_hwfet_modelS_init.csv')
    if T_ess_ref==22:
        if car=='model3':
            if Cy=='udds':
                init        = pd.read_csv(data_dir+'EV_udds_model3_Tess22_init.csv')
            if Cy=='hwfet':
                init        = pd.read_csv(data_dir+'EV_hwfet_model3_Tess22_init.csv')
        if car=='modelS':
            if Cy=='udds':
                init        = pd.read_csv(data_dir+'EV_udds_modelS_Tess22_init.csv')
            if Cy=='hwfet':
                init        = pd.read_csv(data_dir+'EV_hwfet_modelS_Tess22_init.csv')
    if T_ess_ref==28:
        if car=='model3':
            if Cy=='udds':
                init        = pd.read_csv(data_dir+'EV_udds_model3_Tess28_init.csv')
            if Cy=='hwfet':
                init        = pd.read_csv(data_dir+'EV_hwfet_model3_Tess28_init.csv')
        if car=='modelS':
            if Cy=='udds':
                init        = pd.read_csv(data_dir+'EV_udds_modelS_Tess28_init.csv')
            if Cy=='hwfet':
                init        = pd.read_csv(data_dir+'EV_hwfet_modelS_Tess28_init.csv')

    # initial condition assignment for TMS 
    i_term_essT[0]     = init['i_term_essT'].iloc[0]
    i_term_chillTrfg[0]= init['i_term_chillTrfg'].iloc[0]
    i_term_evapTrfg[0] = init['i_term_evapTrfg'].iloc[0]
    i_term_condTrfg[0] = init['i_term_condTrfg'].iloc[0]
    i_term_Tcab[0]     = init['i_term_Tcab'].iloc[0]
    essQheat[0]     = init['essQheat'].iloc[0]
    essT[0]         = init['essT'].iloc[0]
    essTcool[0]     = init['essTcool'].iloc[0]
    essTcoolin[0]   = init['essTcoolin'].iloc[0]
    essTcoolout[0]  = init['essTcoolout'].iloc[0]
    essmdotcool[0]  = init['essmdotcool'].iloc[0]
    coolmdot[0]     = init['coolmdot'].iloc[0]
    coolPpump[0]    = init['coolPpump'].iloc[0]
    rfgPpump[0]     = init['rfgPpump'].iloc[0]
    rfgPcomp[0]     = init['rfgPcomp'].iloc[0]
    chillmdotrfg[0] = init['chillmdotrfg'].iloc[0]
    chillTcool[0]   = init['chillTcool'].iloc[0]
    chillTcoolin[0] = init['chillTcoolin'].iloc[0]
    chillTcoolout[0]= init['chillTcoolout'].iloc[0]
    chillTw[0]      = init['chillTw'].iloc[0]
    chillTrfg[0]    = init['chillTrfg'].iloc[0]
    resTrfgout[0]   = init['resTrfgout'].iloc[0]
    chillTrfgout[0] = init['chillTrfgout'].iloc[0]
    radmdotcool[0]  = init['radmdotcool'].iloc[0]
    radmdotair[0]   = init['radmdotair'].iloc[0]
    radTcool[0]     = init['radTcool'].iloc[0]
    radTcoolin[0]   = init['radTcoolin'].iloc[0]
    radTcoolout[0]  = init['radTcoolout'].iloc[0]
    radTw[0]        = init['radTw'].iloc[0]
    radTair[0]      = init['radTair'].iloc[0]
    radTairin[0]    = init['radTairin'].iloc[0]
    radTairout[0]   = init['radTairout'].iloc[0]
    plowrfg[0]      = init['plowrfg'].iloc[0]
    phighrfg[0]     = init['phighrfg'].iloc[0]
    evapmdotrfg[0]  = init['evapmdotrfg'].iloc[0]
    evapmdotair[0]  = init['evapmdotair'].iloc[0]
    evaphrfg[0]     = init['evaphrfg'].iloc[0]
    evaphrfgin[0]   = init['evaphrfgin'].iloc[0]
    evaphrfgout[0]  = init['evaphrfgout'].iloc[0]
    chillhrfg[0]    = init['chillhrfg'].iloc[0]
    chillhrfgout[0] = init['chillhrfgout'].iloc[0]
    mixhrfg[0]      = init['mixhrfg'].iloc[0]
    condhrfg[0]     = init['condhrfg'].iloc[0]
    condhrfgin[0]   = init['condhrfgin'].iloc[0]
    condhrfgout[0]  = init['condhrfgout'].iloc[0]
    evapTrfg[0]     = init['evapTrfg'].iloc[0]
    evapTrfgin[0]   = init['evapTrfgin'].iloc[0]
    evapTrfgout[0]  = init['evapTrfgout'].iloc[0]
    evapTw[0]       = init['evapTw'].iloc[0]
    evapTair[0]     = init['evapTair'].iloc[0]
    evapTairout[0]  = init['evapTairout'].iloc[0]
    evapPfanair[0]  = init['evapPfanair'].iloc[0]
    condTrfg[0]     = init['condTrfg'].iloc[0]
    condTrfgin[0]   = init['condTrfgin'].iloc[0]
    condTrfgout[0]  = init['condTrfgout'].iloc[0]
    condTw[0]       = init['condTw'].iloc[0]
    condTair[0]     = init['condTair'].iloc[0]
    condTairout[0]  = init['condTairout'].iloc[0]
    condmdotair[0]  = init['condmdotair'].iloc[0]
    resTrfg[0]      = init['resTrfg'].iloc[0]
    reshrfg[0]      = init['reshrfg'].iloc[0]
    cabTair[0]      = init['cabTair'].iloc[0]
    cabTairout[0]   = init['cabTairout'].iloc[0]
    cabmdotair[0]   = init['cabmdotair'].iloc[0]
    cabTw[0]        = init['cabTw'].iloc[0]
    rfgPcomp[0]     = init['rfgPcomp'].iloc[0]
    rfgNcomp[0]     = init['rfgNcomp'].iloc[0]

    ## inputs
    Tamb            = 273+T_init
    radTairin       = np.ones(steps)*(273+T_init)
    plowrfg         = np.ones(steps)*(2*1e5)
    phighrfg        = np.ones(steps)*(14*1e5)
    condTairin      = np.ones(steps)*(273+T_init)
    cabTw           = np.ones(steps)*(273+T_init)        
    ############################
    ###   End of TMS initial ###
    ############################    
    
    ############################
    ###   Loop Through Time  ###
    ############################
    for i in range(1,len(cycSecs)):
        ### Misc calcs
        if veh['noElecAux']=='TRUE':
            auxInKw[i] = veh['auxKw']/veh['altEff']
        else:
            auxInKw[i] = veh['auxKw']+sumPTMS[(i-1)*2]
#             print(i,sumPTMS[(i-1)*2])
            
#             auxInKw[i] = veh['auxKw']

        if soc[i-1]<(veh['minSoc']+veh['percHighAccBuf']):
            reachedBuff[i] = 0
        else:
            reachedBuff[i] = 1
        if soc[i-1]<veh['minSoc'] or (highAccFcOnTag[i-1]==1 and reachedBuff[i]==0):
            highAccFcOnTag[i] = 1
        else:
            highAccFcOnTag[i] = 0
		##### speed? BX
        maxTracMps[i] = mpsAch[i-1]+(maxTracMps2*secs[i])

        ### Component Limits
        curMaxFsKwOut[i] = min( veh['maxFuelStorKw'] , fsKwOutAch[i-1] + ((veh['maxFuelStorKw']/veh['fuelStorSecsToPeakPwr'])*(secs[i])))
        fcTransLimKw[i] = fcKwOutAch[i-1] + ((veh['maxFuelConvKw']/veh['fuelConvSecsToPeakPwr'])*(secs[i]))

        fcMaxKwIn[i] = min(curMaxFsKwOut[i], veh['maxFuelStorKw'])
        fcFsLimKw[i] = veh['fcMaxOutkW']
        curMaxFcKwOut[i] = min(veh['maxFuelConvKw'],fcFsLimKw[i],fcTransLimKw[i])

        if veh['maxEssKwh']==0 or soc[i-1]<veh['minSoc']:
            essCapLimDischgKw[i] = 0.0

        else:
            essCapLimDischgKw[i] = (veh['maxEssKwh']*np.sqrt(veh['essRoundTripEff']))*3600.0*(soc[i-1]-veh['minSoc'])/(secs[i])
        curMaxEssKwOut[i] = min(veh['maxEssKw'],essCapLimDischgKw[i])

        if  veh['maxEssKwh'] == 0 or veh['maxEssKw'] == 0:
            essCapLimChgKw[i] = 0

        else:
            essCapLimChgKw[i] = max(((veh['maxSoc']-soc[i-1])*veh['maxEssKwh']*(1/np.sqrt(veh['essRoundTripEff'])))/((secs[i])*(1/3600.0)),0)

        curMaxEssChgKw[i] = min(essCapLimChgKw[i],veh['maxEssKw'])

        if veh['fcEffType']==4:
            curMaxElecKw[i] = curMaxFcKwOut[i]+curMaxRoadwayChgKw[i]+curMaxEssKwOut[i]-auxInKw[i]

        else:
            curMaxElecKw[i] = curMaxRoadwayChgKw[i]+curMaxEssKwOut[i]-auxInKw[i]
                       
        curMaxAvailElecKw[i] = min(curMaxElecKw[i], veh['mcMaxElecInKw'])

        if curMaxElecKw[i]>0:
            if curMaxAvailElecKw[i] == max(veh['mcKwInArray']):
                mcElecInLimKw[i] = min(veh['mcKwOutArray'][len(veh['mcKwOutArray'])-1],veh['maxMotorKw'])
            else:
                mcElecInLimKw[i] = min(veh['mcKwOutArray'][np.argmax(veh['mcKwInArray']>min(max(veh['mcKwInArray'])-0.01,curMaxAvailElecKw[i]))-1],veh['maxMotorKw'])
        else:
            mcElecInLimKw[i] = 0.0

        mcTransiLimKw[i] = abs(mcMechKwOutAch[i-1])+((veh['maxMotorKw']/veh['motorSecsToPeakPwr'])*(secs[i]))
        curMaxMcKwOut[i] = max(min(mcElecInLimKw[i],mcTransiLimKw[i],veh['maxMotorKw']),-veh['maxMotorKw'])

        if curMaxMcKwOut[i] == 0:
            curMaxMcElecKwIn[i] = 0
        else:
            if curMaxMcKwOut[i] == veh['maxMotorKw']:
                curMaxMcElecKwIn[i] = curMaxMcKwOut[i]/veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]
            else:
                curMaxMcElecKwIn[i] = curMaxMcKwOut[i]/veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,curMaxMcKwOut[i]))-1)]

        if veh['maxMotorKw']==0:
            essLimMcRegenPercKw[i] = 0.0

        else:
            essLimMcRegenPercKw[i] = min((curMaxEssChgKw[i]+auxInKw[i])/veh['maxMotorKw'],1)
        if curMaxEssChgKw[i]==0:
            essLimMcRegenKw[i] = 0.0

        else:
            if veh['maxMotorKw'] == curMaxEssChgKw[i]-curMaxRoadwayChgKw[i]:
                essLimMcRegenKw[i] = min(veh['maxMotorKw'],curMaxEssChgKw[i]/veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1])
            else:
                essLimMcRegenKw[i] = min(veh['maxMotorKw'],curMaxEssChgKw[i]/veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,curMaxEssChgKw[i]-curMaxRoadwayChgKw[i]))-1)])

        curMaxMechMcKwIn[i] = min(essLimMcRegenKw[i],veh['maxMotorKw'])
        curMaxTracKw[i] = (((veh['wheelCoefOfFric']*veh['driveAxleWeightFrac']*veh['vehKg']*gravityMPerSec2)/(1+((veh['vehCgM']*veh['wheelCoefOfFric'])/veh['wheelBaseM'])))/1000.0)*(maxTracMps[i])

        if veh['fcEffType']==4:

            if veh['noElecSys']=='TRUE' or veh['noElecAux']=='TRUE' or highAccFcOnTag[i]==1:
                curMaxTransKwOut[i] = min((curMaxMcKwOut[i]-auxInKw[i])*veh['transEff'],curMaxTracKw[i]/veh['transEff'])
                debug_flag[i] = 1

            else:
                curMaxTransKwOut[i] = min((curMaxMcKwOut[i]-min(curMaxElecKw[i],0))*veh['transEff'],curMaxTracKw[i]/veh['transEff'])
                debug_flag[i] = 2

        else:

            if veh['noElecSys']=='TRUE' or veh['noElecAux']=='TRUE' or highAccFcOnTag[i]==1:
                curMaxTransKwOut[i] = min((curMaxMcKwOut[i]+curMaxFcKwOut[i]-auxInKw[i])*veh['transEff'],curMaxTracKw[i]/veh['transEff'])
                debug_flag[i] = 3

            else:
                curMaxTransKwOut[i] = min((curMaxMcKwOut[i]+curMaxFcKwOut[i]-min(curMaxElecKw[i],0))*veh['transEff'],curMaxTracKw[i]/veh['transEff'])
                debug_flag[i] = 4

        ### Cycle Power
        cycDragKw[i] = 0.5*airDensityKgPerM3*veh['dragCoef']*veh['frontalAreaM2']*(((mpsAch[i-1]+cycMps[i])/2.0)**3)/1000.0
        cycAccelKw[i] = (veh['vehKg']/(2.0*(secs[i])))*((cycMps[i]**2)-(mpsAch[i-1]**2))/1000.0
        cycAscentKw[i] = gravityMPerSec2*np.sin(np.arctan(cycGrade[i]))*veh['vehKg']*((mpsAch[i-1]+cycMps[i])/2.0)/1000.0
        cycTracKwReq[i] = cycDragKw[i] + cycAccelKw[i] + cycAscentKw[i]
        spareTracKw[i] = curMaxTracKw[i]-cycTracKwReq[i]
        cycRrKw[i] = gravityMPerSec2*veh['wheelRrCoef']*veh['vehKg']*((mpsAch[i-1]+cycMps[i])/2.0)/1000.0
        cycWheelRadPerSec[i] = cycMps[i]/veh['wheelRadiusM']
        cycTireInertiaKw[i] = (((0.5)*veh['wheelInertiaKgM2']*(veh['numWheels']*(cycWheelRadPerSec[i]**2.0))/secs[i])-((0.5)*veh['wheelInertiaKgM2']*(veh['numWheels']*((mpsAch[i-1]/veh['wheelRadiusM'])**2.0))/secs[i]))/1000.0

        cycWheelKwReq[i] = cycTracKwReq[i] + cycRrKw[i] + cycTireInertiaKw[i]
        regenContrLimKwPerc[i] = veh['maxRegen']/(1+veh['regenA']*np.exp(-veh['regenB']*((cycMph[i]+mpsAch[i-1]*mphPerMps)/2.0+1-0)))
        cycRegenBrakeKw[i] = max(min(curMaxMechMcKwIn[i]*veh['transEff'],regenContrLimKwPerc[i]*-cycWheelKwReq[i]),0)
        cycFricBrakeKw[i] = -min(cycRegenBrakeKw[i]+cycWheelKwReq[i],0)
        cycTransKwOutReq[i] = cycWheelKwReq[i] + cycFricBrakeKw[i]
#         print(" cycWheelKwReq[i] + cycFricBrakeKw[i], cycTracKwReq[i] + cycRrKw[i] + cycTireInertiaKw[i]:", cycWheelKwReq[i] , cycFricBrakeKw[i],cycTracKwReq[i] , cycRrKw[i] ,cycTireInertiaKw[i])

        if cycTransKwOutReq[i]<=curMaxTransKwOut[i]:
            cycMet[i] = 1
            transKwOutAch[i] = cycTransKwOutReq[i]
        else:
            cycMet[i] = -1
            transKwOutAch[i] = curMaxTransKwOut[i]
        
        ################################
        ###   Speed/Distance Calcs   ###
        ################################

        #Cycle is met
        if cycMet[i]==1:
            mpsAch[i] = cycMps[i]

        #Cycle is not met
        else:
            Drag3 = (1.0/16.0)*airDensityKgPerM3*veh['dragCoef']*veh['frontalAreaM2']
            Accel2 = veh['vehKg']/(2.0*(secs[i]))
            Drag2 = (3.0/16.0)*airDensityKgPerM3*veh['dragCoef']*veh['frontalAreaM2']*mpsAch[i-1]
            Wheel2 = 0.5*veh['wheelInertiaKgM2']*veh['numWheels']/(secs[i]*(veh['wheelRadiusM']**2))
            Drag1 = (3.0/16.0)*airDensityKgPerM3*veh['dragCoef']*veh['frontalAreaM2']*((mpsAch[i-1])**2)
            Roll1 = (gravityMPerSec2*veh['wheelRrCoef']*veh['vehKg']/2.0)
            Ascent1 = (gravityMPerSec2*np.sin(np.arctan(cycGrade[i]))*veh['vehKg']/2.0)
            Accel0 = -(veh['vehKg']*((mpsAch[i-1])**2))/(2.0*(secs[i]))
            Drag0 = (1.0/16.0)*airDensityKgPerM3*veh['dragCoef']*veh['frontalAreaM2']*((mpsAch[i-1])**3)
            Roll0 = (gravityMPerSec2*veh['wheelRrCoef']*veh['vehKg']*mpsAch[i-1]/2.0)
            Ascent0 = (gravityMPerSec2*np.sin(np.arctan(cycGrade[i]))*veh['vehKg']*mpsAch[i-1]/2.0)
            Wheel0 = -((0.5*veh['wheelInertiaKgM2']*veh['numWheels']*(mpsAch[i-1]**2))/(secs[i]*(veh['wheelRadiusM']**2)))

            Total3 = Drag3/1e3
            Total2 = (Accel2+Drag2+Wheel2)/1e3
            Total1 = (Drag1+Roll1+Ascent1)/1e3
            Total0 = (Accel0+Drag0+Roll0+Ascent0+Wheel0)/1e3 - curMaxTransKwOut[i]

            Total = [Total3,Total2,Total1,Total0]
            Total_roots = np.roots(Total)
            ind = np.argmin( abs(cycMps[i] - Total_roots) )
            mpsAch[i] = Total_roots[ind]

        mphAch[i] = mpsAch[i]*mphPerMps
        distMeters[i] = mpsAch[i]*secs[i]
        distMiles[i] = distMeters[i]*(1.0/metersPerMile)

        ### Drive Train
        #print('transKwOutAch[i]:',transKwOutAch[i])
        if transKwOutAch[i]>0:
            transKwInAch[i] = transKwOutAch[i]/veh['transEff']
        else:
            transKwInAch[i] = transKwOutAch[i]*veh['transEff']

        if cycMet[i]==1:

            if veh['fcEffType']==4:
                minMcKw2HelpFc[i] = max(transKwInAch[i], -curMaxMechMcKwIn[i])

            else:
                minMcKw2HelpFc[i] = max(transKwInAch[i] - curMaxFcKwOut[i], -curMaxMechMcKwIn[i])
        else:
            minMcKw2HelpFc[i] = max(curMaxMcKwOut[i], -curMaxMechMcKwIn[i])

        if veh['noElecSys'] == 'TRUE':
            regenBufferSoc[i] = 0

        elif veh['chargingOn']:
            regenBufferSoc[i] = max(veh['maxSoc'] - (maxRegenKwh/veh['maxEssKwh']), (veh['maxSoc']+veh['minSoc'])/2)

        else:
            regenBufferSoc[i] = max(((veh['maxEssKwh']*veh['maxSoc'])-(0.5*veh['vehKg']*(cycMps[i]**2)*(1.0/1000)*(1.0/3600)*veh['motorPeakEff']*veh['maxRegen']))/veh['maxEssKwh'],veh['minSoc'])

        essRegenBufferDischgKw[i] = min(curMaxEssKwOut[i], max(0,(soc[i-1]-regenBufferSoc[i])*veh['maxEssKwh']*3600/secs[i]))

        maxEssRegenBufferChgKw[i] = min(max(0,(regenBufferSoc[i]-soc[i-1])*veh['maxEssKwh']*3600.0/secs[i]),curMaxEssChgKw[i])

        if veh['noElecSys']=='TRUE':
            accelBufferSoc[i] = 0

        else:
            accelBufferSoc[i] = min(max((((((((veh['maxAccelBufferMph']*(1/mphPerMps))**2))-((cycMps[i]**2)))/(((veh['maxAccelBufferMph']*(1/mphPerMps))**2)))*(min(veh['maxAccelBufferPercOfUseableSoc']*(veh['maxSoc']-veh['minSoc']),maxRegenKwh/veh['maxEssKwh'])*veh['maxEssKwh']))/veh['maxEssKwh'])+veh['minSoc'],veh['minSoc']), veh['maxSoc'])

        essAccelBufferChgKw[i] = max(0,(accelBufferSoc[i] - soc[i-1])*veh['maxEssKwh']*3600.0/secs[i])
        maxEssAccelBufferDischgKw[i] = min(max(0, (soc[i-1]-accelBufferSoc[i])*veh['maxEssKwh']*3600/secs[i]),curMaxEssKwOut[i])

        if regenBufferSoc[i] < accelBufferSoc[i]:
            essAccelRegenDischgKw[i] = max(min(((soc[i-1]-(regenBufferSoc[i]+accelBufferSoc[i])/2)*veh['maxEssKwh']*3600.0)/secs[i],curMaxEssKwOut[i]),-curMaxEssChgKw[i])

        elif soc[i-1]>regenBufferSoc[i]:
            essAccelRegenDischgKw[i] = max(min(essRegenBufferDischgKw[i],curMaxEssKwOut[i]),-curMaxEssChgKw[i])

        elif soc[i-1]<accelBufferSoc[i]:
            essAccelRegenDischgKw[i] = max(min(-1.0*essAccelBufferChgKw[i],curMaxEssKwOut[i]),-curMaxEssChgKw[i])

        else:
            essAccelRegenDischgKw[i] = max(min(0,curMaxEssKwOut[i]),-curMaxEssChgKw[i])

        fcKwGapFrEff[i] = abs(transKwOutAch[i]-veh['maxFcEffKw'])

        if veh['noElecSys']=='TRUE':
            mcElectInKwForMaxFcEff[i] = 0

        elif transKwOutAch[i]<veh['maxFcEffKw']:

            if fcKwGapFrEff[i] == veh['maxMotorKw']:
                mcElectInKwForMaxFcEff[i] = fcKwGapFrEff[i]/veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]*-1
            else:
                mcElectInKwForMaxFcEff[i] = fcKwGapFrEff[i]/veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,fcKwGapFrEff[i]))-1)]*-1

        else:

            if fcKwGapFrEff[i] == veh['maxMotorKw']:
                mcElectInKwForMaxFcEff[i] = veh['mcKwInArray'][len(veh['mcKwInArray'])-1]
            else:
                mcElectInKwForMaxFcEff[i] = veh['mcKwInArray'][np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,fcKwGapFrEff[i]))-1]

        if veh['noElecSys']=='TRUE':
            electKwReq4AE[i] = 0

        elif transKwInAch[i] > 0:
            if transKwInAch[i] == veh['maxMotorKw']:
        
                electKwReq4AE[i] = transKwInAch[i]/veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]+auxInKw[i]
            else:
                electKwReq4AE[i] = transKwInAch[i]/veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,transKwInAch[i]))-1)]+auxInKw[i]

        else:
            electKwReq4AE[i] = 0

        prevfcTimeOn[i] = fcTimeOn[i-1]

        if veh['maxFuelConvKw']==0:
            canPowerAllElectrically[i] = accelBufferSoc[i]<soc[i-1] and transKwInAch[i]<=curMaxMcKwOut[i] and (electKwReq4AE[i]<curMaxElecKw[i] or veh['maxFuelConvKw']==0)

        else:
            canPowerAllElectrically[i] = accelBufferSoc[i]<soc[i-1] and transKwInAch[i]<=curMaxMcKwOut[i] and (electKwReq4AE[i]<curMaxElecKw[i] or veh['maxFuelConvKw']==0) and (cycMph[i]-0.00001<=veh['mphFcOn'] or veh['chargingOn']) and electKwReq4AE[i]<=veh['kwDemandFcOn']
        if canPowerAllElectrically[i]:

            if transKwInAch[i]<+auxInKw[i]:
                desiredEssKwOutForAE[i] = auxInKw[i]+transKwInAch[i]

            elif regenBufferSoc[i]<accelBufferSoc[i]:
                desiredEssKwOutForAE[i] = essAccelRegenDischgKw[i]

            elif soc[i-1]>regenBufferSoc[i]:
                desiredEssKwOutForAE[i] = essRegenBufferDischgKw[i]

            elif soc[i-1]<accelBufferSoc[i]:
                desiredEssKwOutForAE[i] = -essAccelBufferChgKw[i]

            else:
                desiredEssKwOutForAE[i] = transKwInAch[i]+auxInKw[i]-curMaxRoadwayChgKw[i]

        else:
            desiredEssKwOutForAE[i] = 0

        if canPowerAllElectrically[i]:
            essAEKwOut[i] = max(-curMaxEssChgKw[i],-maxEssRegenBufferChgKw[i],min(0,curMaxRoadwayChgKw[i]-(transKwInAch[i]+auxInKw[i])),min(curMaxEssKwOut[i],desiredEssKwOutForAE[i]))

        else:
            essAEKwOut[i] = 0

        erAEKwOut[i] = min(max(0,transKwInAch[i]+auxInKw[i]-essAEKwOut[i]),curMaxRoadwayChgKw[i])

        if prevfcTimeOn[i]>0 and prevfcTimeOn[i]<veh['minFcTimeOn']-secs[i]:
            fcForcedOn[i] = True
        else:
            fcForcedOn[i] = False

        ##################################################################################################
        ##################################################################################################
        ##################################################################################################
#         EMS = 'ECMS'
#         EMS = 'rule'
        mcMechKw4ForcedFc[i]=transKwInAch[i]    
            
        ##################################################################################################
        ##################################################################################################
        ##################################################################################################   

        if fcForcedOn[i]==False or canPowerAllElectrically[i]==False:            
            fcForcedState[i] = 1
#             mcMechKw4ForcedFc[i] = (soc[i-1]-soc_ref)*soc_slop
#             mcMechKw4ForcedFc[i] = 0            

        elif transKwInAch[i]<0:            
            fcForcedState[i] = 2
#             mcMechKw4ForcedFc[i] = (soc[i-1]-soc_ref)*soc_slop
#             mcMechKw4ForcedFc[i] = transKwInAch[i]            

        elif veh['maxFcEffKw']==transKwInAch[i]:
            fcForcedState[i] = 3
#             mcMechKw4ForcedFc[i] = (soc[i-1]-soc_ref)*soc_slop
#             mcMechKw4ForcedFc[i] = 0
            
        elif veh['idleFcKw'] > transKwInAch[i] and cycAccelKw[i] >=0:
            fcForcedState[i] = 4
#             mcMechKw4ForcedFc[i] = (soc[i-1]-soc_ref)*soc_slop
#             mcMechKw4ForcedFc[i] = transKwInAch[i] - veh['idleFcKw']
            
        elif veh['maxFcEffKw']>transKwInAch[i]:
            fcForcedState[i] = 5
#             mcMechKw4ForcedFc[i] = (soc[i-1]-soc_ref)*soc_slop
#             mcMechKw4ForcedFc[i] = 0
            
        else:
            ##################################################################################################        
            ## max ICE efficiency, the rest is filled in by EM
            fcForcedState[i] = 6            
#             mcMechKw4ForcedFc[i] = (soc[i-1]-soc_ref)*soc_slop          
#             mcMechKw4ForcedFc[i] = transKwInAch[i] - veh['maxFcEffKw']            
        
        if (-mcElectInKwForMaxFcEff[i]-curMaxRoadwayChgKw[i])>0:
            essDesiredKw4FcEff[i] = (-mcElectInKwForMaxFcEff[i]-curMaxRoadwayChgKw[i]) * veh['essDischgToFcMaxEffPerc']

        else:
            essDesiredKw4FcEff[i] = (-mcElectInKwForMaxFcEff[i]-curMaxRoadwayChgKw[i]) * veh['essChgToFcMaxEffPerc']

        if accelBufferSoc[i]>regenBufferSoc[i]:
            essKwIfFcIsReq[i] = min(curMaxEssKwOut[i],veh['mcMaxElecInKw']+auxInKw[i],curMaxMcElecKwIn[i]+auxInKw[i], max(-curMaxEssChgKw[i], essAccelRegenDischgKw[i]))

        elif essRegenBufferDischgKw[i]>0:
            essKwIfFcIsReq[i] = min(curMaxEssKwOut[i],veh['mcMaxElecInKw']+auxInKw[i],curMaxMcElecKwIn[i]+auxInKw[i], max(-curMaxEssChgKw[i], min(essAccelRegenDischgKw[i],mcElecInLimKw[i]+auxInKw[i], max(essRegenBufferDischgKw[i],essDesiredKw4FcEff[i]))))

        elif essAccelBufferChgKw[i]>0:
            essKwIfFcIsReq[i] = min(curMaxEssKwOut[i],veh['mcMaxElecInKw']+auxInKw[i],curMaxMcElecKwIn[i]+auxInKw[i], max(-curMaxEssChgKw[i], max(-1*maxEssRegenBufferChgKw[i], min(-essAccelBufferChgKw[i],essDesiredKw4FcEff[i]))))


        elif essDesiredKw4FcEff[i]>0:
            essKwIfFcIsReq[i] = min(curMaxEssKwOut[i],veh['mcMaxElecInKw']+auxInKw[i],curMaxMcElecKwIn[i]+auxInKw[i], max(-curMaxEssChgKw[i], min(essDesiredKw4FcEff[i],maxEssAccelBufferDischgKw[i])))

        else:
            essKwIfFcIsReq[i] = min(curMaxEssKwOut[i],veh['mcMaxElecInKw']+auxInKw[i],curMaxMcElecKwIn[i]+auxInKw[i], max(-curMaxEssChgKw[i], max(essDesiredKw4FcEff[i],-maxEssRegenBufferChgKw[i])))

        erKwIfFcIsReq[i] = max(0,min(curMaxRoadwayChgKw[i],curMaxMechMcKwIn[i],essKwIfFcIsReq[i]-mcElecInLimKw[i]+auxInKw[i]))

        mcElecKwInIfFcIsReq[i] = essKwIfFcIsReq[i]+erKwIfFcIsReq[i]-auxInKw[i]

        if veh['noElecSys']=='TRUE':
            mcKwIfFcIsReq[i] = 0

        elif  mcElecKwInIfFcIsReq[i] == 0:
            mcKwIfFcIsReq[i] = 0

        elif mcElecKwInIfFcIsReq[i] > 0:

            if mcElecKwInIfFcIsReq[i] == max(veh['mcKwInArray']):
                 mcKwIfFcIsReq[i] = mcElecKwInIfFcIsReq[i]*veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]
            else:
                 mcKwIfFcIsReq[i] = mcElecKwInIfFcIsReq[i]*veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwInArray']>min(max(veh['mcKwInArray'])-0.01,mcElecKwInIfFcIsReq[i]))-1)]

        else:
            if mcElecKwInIfFcIsReq[i]*-1 == max(veh['mcKwInArray']):
                mcKwIfFcIsReq[i] = mcElecKwInIfFcIsReq[i]/veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]
            else:
                mcKwIfFcIsReq[i] = mcElecKwInIfFcIsReq[i]/(veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwInArray']>min(max(veh['mcKwInArray'])-0.01,mcElecKwInIfFcIsReq[i]*-1))-1)])

        if veh['maxMotorKw']==0:
            mcMechKwOutAch[i] = 0

        elif fcForcedOn[i]==True and canPowerAllElectrically[i]==True and (veh['vehPtType']==2.0 or veh['vehPtType']==3.0) and veh['fcEffType']!=4:
            mcMechKwOutAch[i] =  mcMechKw4ForcedFc[i]

        elif transKwInAch[i]<=0:
            if veh['fcEffType']!=4 and veh['maxFuelConvKw']> 0:
                if canPowerAllElectrically[i] == 1:
                    mcMechKwOutAch[i] = -min(curMaxMechMcKwIn[i],-transKwInAch[i])
                else:
                    mcMechKwOutAch[i] = min(-min(curMaxMechMcKwIn[i], -transKwInAch[i]),max(-curMaxFcKwOut[i], mcKwIfFcIsReq[i]))
            else:
                mcMechKwOutAch[i] = min(-min(curMaxMechMcKwIn[i],-transKwInAch[i]),-transKwInAch[i])

        elif canPowerAllElectrically[i] == 1:
            mcMechKwOutAch[i] = transKwInAch[i]

        else:
            ################################################################################################
            ################################################################################################
            ################################################################################################
            #  
            mcMechKwOutAch[i] = mcMechKw4ForcedFc[i]
#             mcMechKwOutAch[i] = max(minMcKw2HelpFc[i],mcKwIfFcIsReq[i])
        #print('mcMechKwOutAch[i]:',mcMechKwOutAch[i])
        if mcMechKwOutAch[i]==0:
            mcElecKwInAch[i] = 0.0
            motor_index_debug[i] = 0

        elif mcMechKwOutAch[i]<0:

            if mcMechKwOutAch[i]*-1 == max(veh['mcKwInArray']):
                mcElecKwInAch[i] = mcMechKwOutAch[i]*veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]
            else:
                mcElecKwInAch[i] = mcMechKwOutAch[i]*veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwInArray']>min(max(veh['mcKwInArray'])-0.01,mcMechKwOutAch[i]*-1))-1)]

        else:
            if veh['maxMotorKw'] == mcMechKwOutAch[i]:
                mcElecKwInAch[i] = mcMechKwOutAch[i]/veh['mcFullEffArray'][len(veh['mcFullEffArray'])-1]
            else:                
                mcElecKwInAch[i] = mcMechKwOutAch[i]/veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,mcMechKwOutAch[i]))-1)]
#                 print("i mcElecKwInAch[i],veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,mcMechKwOutAch[i]))-1)]:",i-2,mcMechKwOutAch[i],veh['mcFullEffArray'][max(1,np.argmax(veh['mcKwOutArray']>min(veh['maxMotorKw']-0.01,mcMechKwOutAch[i]))-1)],mcElecKwInAch[i])

        if curMaxRoadwayChgKw[i] == 0:
            roadwayChgKwOutAch[i] = 0

        elif veh['fcEffType']==4:
            roadwayChgKwOutAch[i] = max(0, mcElecKwInAch[i], maxEssRegenBufferChgKw[i], essRegenBufferDischgKw[i], curMaxRoadwayChgKw[i])

        elif canPowerAllElectrically[i] == 1:
            roadwayChgKwOutAch[i] = erAEKwOut[i]

        else:
            roadwayChgKwOutAch[i] = erKwIfFcIsReq[i]

        minEssKw2HelpFc[i] = mcElecKwInAch[i] + auxInKw[i] - curMaxFcKwOut[i] - roadwayChgKwOutAch[i]

        if veh['maxEssKw'] == 0 or veh['maxEssKwh'] == 0:
            essKwOutAch[i]  = 0

        elif veh['fcEffType']==4:

            if transKwOutAch[i]>=0:
                essKwOutAch[i] = min(max(minEssKw2HelpFc[i],essDesiredKw4FcEff[i],essAccelRegenDischgKw[i]),curMaxEssKwOut[i],mcElecKwInAch[i]+auxInKw[i]-roadwayChgKwOutAch[i])

            else:
                essKwOutAch[i] = mcElecKwInAch[i]+auxInKw[i]-roadwayChgKwOutAch[i]

        elif highAccFcOnTag[i]==1 or veh['noElecAux']=='TRUE':
            essKwOutAch[i] = mcElecKwInAch[i]-roadwayChgKwOutAch[i]
        else:
            #print('mcElecKwInAch[i],auxInKw[i],roadwayChgKwOutAch[i]',mcElecKwInAch[i],auxInKw[i],roadwayChgKwOutAch[i])
            essKwOutAch[i] = mcElecKwInAch[i]+auxInKw[i]-roadwayChgKwOutAch[i]

        if veh['maxFuelConvKw']==0:
            fcKwOutAch[i] = 0

        elif veh['fcEffType']==4:
            fcKwOutAch[i] = min(curMaxFcKwOut[i], max(0, mcElecKwInAch[i]+auxInKw[i]-essKwOutAch[i]-roadwayChgKwOutAch[i]))

        elif veh['noElecSys']=='TRUE' or veh['noElecAux']=='TRUE' or highAccFcOnTag[i]==1:
            fcKwOutAch[i] = min(curMaxFcKwOut[i], max(0, transKwInAch[i]-mcMechKwOutAch[i]+auxInKw[i]))

        else:
            fcKwOutAch[i] = min(curMaxFcKwOut[i], max(0, transKwInAch[i]-mcMechKwOutAch[i]))

        if fcKwOutAch[i]==0:
            fcKwInAch[i] = 0.0
            fcKwOutAch_pct[i] = 0

        if veh['maxFuelConvKw'] == 0:
            fcKwOutAch_pct[i] = 0
        else:
            fcKwOutAch_pct[i] = fcKwOutAch[i] / veh['maxFuelConvKw']

        if fcKwOutAch[i] == 0:
            fcKwInAch[i] = 0
        else:
            if fcKwOutAch[i] == veh['fcMaxOutkW']:
                fcKwInAch[i] = fcKwOutAch[i]/f(fcKwOutAch[i]/veh['maxFuelConvKw'])
#                 fcKwInAch[i] = fcKwOutAch[i]/veh['fcEffArray'][len(veh['fcEffArray'])-1]
            else:
                fcKwInAch[i] = fcKwOutAch[i]/f(fcKwOutAch[i]/veh['maxFuelConvKw'])
#                 fcKwInAch[i] = fcKwOutAch[i]/(veh['fcEffArray'][max(1,np.argmax(veh['fcKwOutArray']>min(fcKwOutAch[i],veh['fcMaxOutkW']-0.001))-1)])

        fsKwOutAch[i] = np.copy( fcKwInAch[i] )

        fsKwhOutAch[i] = fsKwOutAch[i]*secs[i]*(1/3600.0)


        if veh['noElecSys']=='TRUE':
            essCurKwh[i] = 0

        elif essKwOutAch[i]<0:
            essCurKwh[i] = essCurKwh[i-1]-essKwOutAch[i]*(secs[i]/3600.0)*np.sqrt(veh['essRoundTripEff'])

        else:
            essCurKwh[i] = essCurKwh[i-1]-essKwOutAch[i]*(secs[i]/3600.0)*(1/np.sqrt(veh['essRoundTripEff']))

        if veh['maxEssKwh']==0:
            soc[i] = 0.0

        else:
            soc[i] = essCurKwh[i]/veh['maxEssKwh']

        if canPowerAllElectrically[i]==True and fcForcedOn[i]==False and fcKwOutAch[i]==0:
            fcTimeOn[i] = 0
        else:
            fcTimeOn[i] = fcTimeOn[i-1] + secs[i]

        ### Battery wear calcs

        if veh['noElecSys']!='TRUE':

            if essCurKwh[i]>essCurKwh[i-1]:
                addKwh[i] = (essCurKwh[i]-essCurKwh[i-1]) + addKwh[i-1]
            else:
                addKwh[i] = 0

            if addKwh[i]==0:
                dodCycs[i] = addKwh[i-1]/veh['maxEssKwh']
            else:
                dodCycs[i] = 0

            if dodCycs[i]!=0:
                essPercDeadArray[i] = np.power(veh['essLifeCoefA'],1.0/veh['essLifeCoefB']) / np.power(dodCycs[i],1.0/veh['essLifeCoefB'])
            else:
                essPercDeadArray[i] = 0

        ### Energy Audit Calculations
        dragKw[i] = 0.5*airDensityKgPerM3*veh['dragCoef']*veh['frontalAreaM2']*(((mpsAch[i-1]+mpsAch[i])/2.0)**3)/1000.0
        if veh['maxEssKw'] == 0 or veh['maxEssKwh']==0:
            essLossKw[i] = 0
        elif essKwOutAch[i]<0:
            essLossKw[i] = -essKwOutAch[i] - (-essKwOutAch[i]*np.sqrt(veh['essRoundTripEff']))
        else:
            essLossKw[i] = essKwOutAch[i]*(1.0/np.sqrt(veh['essRoundTripEff']))-essKwOutAch[i]
        accelKw[i] = (veh['vehKg']/(2.0*(secs[i])))*((mpsAch[i]**2)-(mpsAch[i-1]**2))/1000.0
        ascentKw[i] = gravityMPerSec2*np.sin(np.arctan(cycGrade[i]))*veh['vehKg']*((mpsAch[i-1]+mpsAch[i])/2.0)/1000.0
        rrKw[i] = gravityMPerSec2*veh['wheelRrCoef']*veh['vehKg']*((mpsAch[i-1]+mpsAch[i])/2.0)/1000.0
        
        ### TMS
        ii=i
        # battery current
        I = np.abs(essKwOutAch[ii])/N_cell*1e3/V_nominal
        for i in range((i-1)*int(1/dt)+1,(i-1)*int(1/dt)+1+int(1/dt)):                       
            ### control - essT
#             if i<int(100/dt):
#                 essT_ref           = 25+273
#             else:
# #                 if i<int(1200/dt):
#                 essT_ref           = 24+273
# #                 else:
# #                     essT_ref           = 25+273 
            essT_ref            = T_ess_ref+273
#             essT_ref            = 25+273    
            e_essT[i]           = essT[i-1]-essT_ref
            p_term_essT[i]      = kp_essT*e_essT[i]
            i_term_essT[i]      = i_term_essT[i-1]+ki_essT*e_essT[i]*dt
            PID_u_essT          = p_term_essT[i]+i_term_essT[i]
            essmdotcool[i]      = max(0.0001,200/3600*PID_u_essT)

            ### control - chill superheat
            T_superheat_ref     = 3
            e_chillTrfg[i]      = chillTrfgout[i-1]-chillTrfg[i-1]-T_superheat_ref
            p_term_chillTrfg[i] = kp_chillTrfg*e_chillTrfg[i]
            i_term_chillTrfg[i] = i_term_chillTrfg[i-1]+ki_chillTrfg*e_chillTrfg[i]*dt
            PID_u_chillTrfg     = p_term_chillTrfg[i]+i_term_chillTrfg[i]
            rfgNcomp[i]         =  max(0.05,PID_u_chillTrfg)

            ### control - evap superheat
            e_evapTrfg[i]       = evapTrfgout[i-1]-evapTrfg[i-1]-T_superheat_ref
            p_term_evapTrfg[i]  = kp_evapTrfg*e_evapTrfg[i]
            i_term_evapTrfg[i]  = i_term_evapTrfg[i-1]+ki_evapTrfg*e_evapTrfg[i]*dt
            PID_u_evapTrfg      = p_term_evapTrfg[i]+i_term_evapTrfg[i]
            if AC_switch=='on':
                evapmdotrfg[i]  = max(0.002,PID_u_evapTrfg)
            else:
                evapmdotrfg[i]  = 0

            ### control - cond subcool
            T_subcool_ref       = -5
            e_condTrfg[i]       = condTrfgout[i-1]-condTrfg[i-1]-T_subcool_ref
            p_term_condTrfg[i]  = kp_condTrfg*e_condTrfg[i]
            i_term_condTrfg[i]  = i_term_condTrfg[i-1]+ki_condTrfg*e_condTrfg[i]*dt
            PID_u_condTrfg      = p_term_condTrfg[i]+i_term_condTrfg[i]    
            condmdotair[i]      =  max(0.001,PID_u_condTrfg)

            ### control - Tcab
#             if i<int(1500/dt):
#                 T_cab_ref           = 25+273
#             else:
#                 if i<int(1650/dt):
#                     T_cab_ref           = 22+273
#                 else:
#                     T_cab_ref           = 25+273      
            T_cab_ref           = T_cab_REF+273
            e_Tcab[i]           = cabTair[i-1]-T_cab_ref
            p_term_Tcab[i]      = kp_Tcab*e_Tcab[i]
            i_term_Tcab[i]      = i_term_Tcab[i-1]+ki_Tcab*e_Tcab[i]*dt
            PID_u_Tcab          = p_term_Tcab[i]+i_term_Tcab[i]    
            if AC_switch=='on':
                evapmdotair[i]  = max(0.005,PID_u_Tcab)
            else:
                evapmdotair[i]  = 0

            ### ess and cool loop
            essQheat[i-1]  = N_cell*I**2*R
            esshcool       = 140   # [w/m^2K]
            essAcool       = N_cell*0.01*0.065   # [m^2]
            essTcoolin[i]  = chillTcoolout[i-1]
            essT[i]        = essT[i-1] + (essQheat[i-1] - essHXeff*esshcool*essAcool*(essT[i-1]-essTcool[i-1]))*dt/essm/essc
            essTcool[i]    = essTcool[i-1] + (essmdotcool[i]*essccool*(essTcoolin[i]-essTcoolout[i-1]) + \
                                         esshcool*essAcool*(essT[i-1]-essTcool[i-1]))*dt/essmcool/essccool
            essTcoolout[i] = essTcool[i]

            ## compressor
            rfgdensincomp  = CP.PropsSI('D','P',plowrfg[i],'Q',1,fluid)
            coolPpump[i]   = essmdotcool[i]*(coolpdownpump-coolpuppump)/cooldens/cooliseffpump
            rp             = phighrfg[i]/plowrfg[i]
            rfgiseffcomp   = min(0.8,max(0.5,(a1+a2*rp)+(a3+a4*rp)*rfgNcomp[i]+(a5+a6*rp)*rfgNcomp[i]**2))
            rfgvoleffcomp  = min(0.9,max(0.75,(b1+b2*rp)+(b3+b4*rp)*rfgNcomp[i]+(b5+b6*rp)*rfgNcomp[i]**2))
            m_dot          = max(evapmdotrfg[i-1]+0.00001,rfgvoleffcomp*(rfgNcomp[i]*rfgVcomp*1e3*rfgdensincomp*60*1e-6)/3600)
            chillmdotrfg[i]= m_dot-evapmdotrfg[i]
            radmdotcool[i] = essmdotcool[i]
            coolmdot[i]    = essmdotcool[i]
            
            ## radiator of the liquid cooling loop
            # radiator is turned off in this study due to the high ambient temperature
            radmdotair[i]  = 0.0001 # [kg/s]  fan power - 2014 JPS-experimental investigation on heat pipe cooling for HEV and EV lithium-ion battery
            radTcoolin[i]  = essTcoolout[i]
            radTcool[i]    = radTcool[i-1] + (radmdotcool[i]*radccool* (radTcoolin[i]-radTcoolout[i-1]) \
                                           + radhcoolw*radAcoolw*(radTw[i-1]-radTcool[i-1]))*dt/radmcool/radccool
            radTcoolout[i] = radTcool[i]
            radTw[i]       = radTw[i-1] + (-radeff*radhcoolw*radAcoolw* (radTw[i-1]-radTcool[i-1]) \
                                           - radhairw*radAairw*(radTw[i-1]-radTair[i-1]))*dt/radmw/radcw
            radTairout[i]  = (radhairw*radAairw*radTw[i]+(-0.5*radhairw*radAairw+radmdotair[i]*radcair)*Tamb)\
                                           /(0.5*radhairw*radAairw+radmdotair[i]*radcair)
            radTair[i]     = 0.5*radTairout[i]+0.5*Tamb
            
            ## chiller of the liquid cooling loop
            chillTcoolin[i]= radTcoolout[i]
            chillTcool[i]  = chillTcool[i-1] + (coolmdot[i]*chillccool*    (chillTcoolin[i]-chillTcoolout[i-1]) \
                                           + chillhcoolw*chillAcoolw*(chillTw[i-1]-chillTcool[i-1]))*dt/chillmcool/chillccool    
            chillTcoolout[i]=chillTcool[i]
            chillTw[i]     =chillTw[i-1] + (-chillhcoolw*chillAcoolw*(chillTw[i-1]-chillTcool[i-1])  \
                                       -chilleff*chillhrfgw*chillArfgw*  (chillTw[i-1]-chillTrfg[i-1]))*dt/chillmw/chillcw   

            ## refrigerant loop
            rfgNpump       = 1000  # [rpm]
            rfgppumpout    = plowrfg[0]     # [bar]
            rfgppumpin     = 1.8*1e5   # [bar]
            rfgpcompout    = phighrfg[0]     # [bar]
            rfgpcompin     = plowrfg[0]   # [bar]
            rfgdens        = 1000  # [kg/m^3]
            rfgiseffpump   = 0.7
            rfgiseffexp    = 0.4 # 2021 ECM A novel electric vehicle thermal management system based on cooling and heating of batteries by refrigerant
            evapPfanair[i] = evapmdotair[i]*0.05/0.03*1e3 # [W] # 0.05kW 0.03kg/s air - fan power - 2002 - Thermal Evaluation of Toyota Prius Battery Pack
            condswitchfanair=0
            condPfanair    =    condswitchfanair*0.15 # [kW] 2011 - Evaluation of Impact of Active Grille Shutter on Vehicle Thermal Management
            cabmdotair[i]  = evapmdotair[i]    # [kg/s]

            hmin           = 100e3
            hmax           = 700e3

            evapTairout[i] = (evaphairw*evapAairw*evapTw[i-1]+evapmdotair[i]*evapcair*cabTair[i-1])\
                            /(evaphairw*evapAairw            +evapmdotair[i]*evapcair)
            evapTair[i]    = evapTairout[i]
            q_cab      =  R_sun*A_sun + K_amb*(Tamb - cabTair[i-1]) + M_driver*A_driver
            cabTair[i]     = cabTair[i-1] + (cabmdotair[i]*cabcair* (evapTairout[i]-cabTair[i-1]) \
                                           + q_cab)*dt/cabmair/cabcair
            evapTw[i]      = evapTw[i-1] + (-evapeff*evaphrfgw*evapArfgw* (evapTw[i-1]-evapTrfg[i-1]) \
                                           - evaphairw*evapAairw*(evapTw[i-1]-evapTair[i]))*dt/evapmw/evapcw

            reshrfg[i]     = min(hmax,reshrfg[i-1] + evapmdotrfg[i]*(condhrfgout[i-1]-reshrfg[i-1])*dt/resmrfg)
            resTrfg[i]     = CP.PropsSI('T','P',plowrfg[i],'H',reshrfg[i],fluid)

            ressrfgout     = CP.PropsSI('S','P',phighrfg[i],'H',reshrfg[i-1],fluid)

            evaphisrfgin   = CP.PropsSI('H','P',plowrfg[i],'S',ressrfgout,fluid)

            evaphrfgin[i]  = reshrfg[i-1] - 1/rfgiseffexp*(reshrfg[i-1]-evaphisrfgin)
            evaphrfg[i]    = min(hmax,evaphrfg[i-1] + (evapmdotrfg[i]*(evaphrfgin[i]-evaphrfgout[i-1]) \
                                           + evaphrfgw*evapArfgw*(evapTw[i]-evapTrfg[i-1]))*dt/evapmrfg)   # 2021 ECM A novel electric vehicle thermal management system based on cooling and heating of batteries by refrigerant
            evaphrfgout[i] = min(hmax,max(hmin,2*evaphrfg[i]-evaphrfgin[i]))
            evapTrfgin[i]  = CP.PropsSI('T','P',plowrfg[i],'H',evaphrfgin[i],fluid)
            evapTrfg[i]    = CP.PropsSI('T','P',plowrfg[i],'H',evaphrfg[i],fluid)
            evapTrfgout[i] = CP.PropsSI('T','P',plowrfg[i],'H',evaphrfgout[i],fluid)

            chillhrfg[i]   = min(hmax,chillhrfg[i-1] + (chillmdotrfg[i]*(evaphrfgin[i]-chillhrfgout[i-1]) \
                                           + chillhrfgw*chillArfgw*(chillTw[i]-chillTrfg[i-1]))*dt/chillmrfg)
            chillhrfgout[i]= min(hmax,max(hmin,2*chillhrfg[i]-evaphrfgin[i]))
            chillTrfg[i]   = CP.PropsSI('T','P',plowrfg[i],'H',chillhrfg[i],fluid)
            chillTrfgout[i]= CP.PropsSI('T','P',plowrfg[i],'H',chillhrfgout[i],fluid)

            condTairout[i] = (condhairw*condAairw*condTw[i-1]+(-0.5*condhairw*condAairw+condmdotair[i]*condcair)*Tamb)\
                                           /(0.5*condhairw*condAairw+condmdotair[i]*condcair)
            condTair[i]    = 0.5*condTairout[i]+0.5*Tamb
            condTw[i]      = condTw[i-1] + (-condeff*condhrfgw*condArfgw* (condTw[i-1]-condTrfg[i-1]) \
                                           - condhairw*condAairw*(condTw[i-1]-condTair[i]))*dt/condmw/condcw    

            mixhrfg[i]     = (evapmdotrfg[i]*evaphrfgout[i]+chillmdotrfg[i]*chillhrfgout[i])/(evapmdotrfg[i]+chillmdotrfg[i])
            mixsrfgout     = CP.PropsSI('S','P',plowrfg[i],'H',mixhrfg[i],fluid)
            mixhisrfgout   = CP.PropsSI('H','P',phighrfg[i],'S',mixsrfgout,fluid)

            condhrfgin[i]  = mixhrfg[i] - rfgiseffcomp*(mixhrfg[i]-mixhisrfgout)
            condTrfgin[i]  = CP.PropsSI('T','P',phighrfg[i],'H',condhrfgin[i],fluid)
            condhrfg[i]    = min(hmax,condhrfg[i-1] + ((evapmdotrfg[i]+chillmdotrfg[i])* (condhrfgin[i]-condhrfgout[i-1]) \
                                           + condhrfgw*condArfgw*(condTw[i]-condTrfg[i-1]))*dt/condmrfg)  # 2021 ECM A novel electric vehicle thermal management system based on cooling and heating of batteries by refrigerant
            condhrfgout[i] = max(hmin,min(hmax,2*condhrfg[i]-condhrfgin[i]))
            condTrfg[i]    = CP.PropsSI('T','P',phighrfg[i],'H',condhrfg[i],fluid)
            condTrfgout[i] = CP.PropsSI('T','P',phighrfg[i],'H',condhrfgout[i],fluid)

            rfgPcomp[i]    = (evapmdotrfg[i]+chillmdotrfg[i])*(phighrfg[i]-plowrfg[i])/rfgiseffcomp/rfgdensincomp
            sumPTMS[i]     = (rfgPcomp[i] + evapPfanair[i] + coolPpump[i])/1e3
            
            ## real-time status display
            if i%(300/dt)==0:
                print('!!Comp/fan/colpum/N','{:<4}'.format(round(rfgPcomp[i]/1e3,3)),'{:<4}'.format(round(evapPfanair[i]/1e3,3)),\
                      '{:<4}'.format(round(coolPpump[i]/1e3,3)),'{:>6}'.format(round(rfgNcomp[i],1)),'mdot/eff:',\
                      int((evapmdotrfg[i])*3600),'{:>2}'.format(int((chillmdotrfg[i])*3600)),'{:<4}'.format(round(rfgiseffcomp,2)),
                      '{:<4}'.format(round(rfgvoleffcomp,2)),'CPU(s):','{:>2}'.format(int(time.time()-t0)),'Sim(s):',\
                      '{:>4}'.format(int(i*dt)),'Veh P','{:>6}'.format(round(essKwOutAch[ii],2)))
        i=ii   

        
    print("Time(s):",round(time.time()-t0,2))
    #  plot
    from plotly.subplots import make_subplots
    import plotly.express as px
    import plotly.graph_objects as go
    fig1 = make_subplots(rows=5, cols=5,subplot_titles=\
    ("essT","chillTrfgout","evapTrfgout","condTrfgout","cabTair",\
     "e_essT","e_chillTrfg","e_evapTrfg","e_condTrfg","e_Tcab",\
     "essmdotcool","chillmdotrfg","evapmdotrfg","condmdotair","evapmdotair",\
     "p_term_essT",'p_term_chillTrfg','p_term_evapTrfg','p_term_condTrfg','p_term_Tcab',\
     'i_term_essT','i_term_chillTrfg','i_term_evapTrfg',"i_term_condTrfg","i_term_Tcab"))
    T=i
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=essT-273)         ,row=1, col=1)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=chillTrfgout-273) ,row=1, col=2)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=evapTrfgout-273)  ,row=1, col=3)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=condTrfgout-273)  ,row=1, col=4)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=cabTair-273)      ,row=1, col=5)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=e_essT)           ,row=2, col=1)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=e_chillTrfg)      ,row=2, col=2)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=e_evapTrfg)       ,row=2, col=3)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=e_condTrfg)       ,row=2, col=4)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=e_Tcab)           ,row=2, col=5)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=essmdotcool)      ,row=3, col=1)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=chillmdotrfg)     ,row=3, col=2)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=evapmdotrfg)      ,row=3, col=3)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=condmdotair)      ,row=3, col=4)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=evapmdotair)      ,row=3, col=5)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=p_term_essT)      ,row=4, col=1)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=p_term_chillTrfg) ,row=4, col=2)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=p_term_evapTrfg)  ,row=4, col=3)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=p_term_condTrfg)  ,row=4, col=4)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=p_term_Tcab)      ,row=4, col=5)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=i_term_essT)      ,row=5, col=1)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=i_term_chillTrfg) ,row=5, col=2)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=i_term_evapTrfg)  ,row=5, col=3)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=i_term_condTrfg)  ,row=5, col=4)
    fig1.add_trace(go.Scatter(x=np.linspace(0,T,len(essT)), y=i_term_Tcab)      ,row=5, col=5)

    fig1.update_yaxes(range=[20, 40],  row=1, col=1)
    fig1.update_yaxes(range=[-10, 0],  row=1, col=2)
    fig1.update_yaxes(range=[-10, 0],  row=1, col=3)
    fig1.update_yaxes(range=[30, 60],  row=1, col=4)
    fig1.update_yaxes(range=[20, 30],  row=1, col=5)
    fig1.update_yaxes(range=[0, 0.05], row=3, col=1)
    fig1.update_yaxes(range=[0, 0.01], row=3, col=2)
    fig1.update_yaxes(range=[0, 0.05], row=3, col=3)
    fig1.update_yaxes(range=[0, 1.5],  row=3, col=4)
    fig1.update_yaxes(range=[0, 0.1],  row=3, col=5)

    fig1.update_layout(height=700, width=1000, title_text="Side By Side Subplots")
    fig1.show()

    import pygame
    pygame.init()
    pygame.mixer.init()
    sounda= pygame.mixer.Sound("beep.wav")
    sounda.play()

    ############################################
    ### Calculate Results and Assign Outputs ###
    ############################################

    #########################################################
    essKwh = essCurKwh[0]-essCurKwh[-1]
    print('sum(fsKwhOutAch),essKwh',sum(fsKwhOutAch),essKwh)
    #########################################################
    
    output = dict()

    if sum(fsKwhOutAch) == 0:
        output['mpgge'] = 0

    else:
        output['mpgge'] = sum(distMiles)/(sum(fsKwhOutAch)*(1/kWhPerGGE))
        output['mpge_cor']=sum(distMiles)/((sum(fsKwhOutAch)+essKwh)*(1/kWhPerGGE))

    roadwayChgKj = sum(roadwayChgKwOutAch*secs)
    essDischKj = -(soc[-1]-initSoc)*veh['maxEssKwh']*3600.0
    output['battery_kWh_per_mi'] = (essDischKj/3600.0) / sum(distMiles)
    output['electric_kWh_per_mi'] = ((roadwayChgKj+essDischKj)/3600.0) / sum(distMiles)
    output['maxTraceMissMph'] = mphPerMps*max(abs(cycMps-mpsAch))
    fuelKj = sum(np.asarray(fsKwOutAch)*np.asarray(secs))
    roadwayChgKj = sum(np.asarray(roadwayChgKwOutAch)*np.asarray(secs))
    essDischgKj = -(soc[-1]-initSoc)*veh['maxEssKwh']*3600.0

    if (fuelKj+roadwayChgKj)==0:
        output['ess2fuelKwh'] = 1.0

    else:
        output['ess2fuelKwh'] = essDischgKj/(fuelKj+roadwayChgKj)

    fuelKg=np.asarray(fsKwhOutAch)/veh['fuelKwhPerKg']
    fuelKgAch=np.zeros(len(fuelKg))
    fuelKgAch[0]=fuelKg[0]
    for qw1 in range(1,len(fuelKg)):
        fuelKgAch[qw1]=fuelKgAch[qw1-1]+fuelKg[qw1]
    output['initial_soc'] = soc[0]
    output['final_soc'] = soc[-1]


    if output['mpgge'] == 0:
        Gallons_gas_equivalent_per_mile = output['electric_kWh_per_mi']/33.7

    else:
         Gallons_gas_equivalent_per_mile = 1/output['mpgge'] + output['electric_kWh_per_mi']/33.7

    output['mpgge_elec'] = 1/Gallons_gas_equivalent_per_mile
    output['soc'] = np.asarray(soc)
    output['distance_mi'] = sum(distMiles)
    duration_sec = cycSecs[-1]-cycSecs[0]
    output['avg_speed_mph'] = sum(distMiles) / (duration_sec/3600.0)
    accel = np.diff(mphAch) / np.diff(cycSecs)
    output['avg_accel_mphps'] = np.mean(accel[accel>0])

    if max(mphAch)>60:
        output['ZeroToSixtyTime_secs'] = np.interp(60,mphAch,cycSecs)
    else:
        output['ZeroToSixtyTime_secs'] = 0.0

    #######################################################################
    ####  Time series information for additional analysis / debugging. ####
    ####             Add parameters of interest as needed.             ####
    #######################################################################

    output['fcKwOutAch'] = np.asarray(fcKwOutAch)
    output['fsKwhOutAch'] = np.asarray(fsKwhOutAch)
    output['fcKwInAch'] = np.asarray(fcKwInAch)
    output['essKwOutAch'] = np.asarray(essKwOutAch)
    output['time'] = np.asarray(cycSecs)
    output['fcForcedState'] = np.asarray(fcForcedState)
    output['transKwInAch'] = np.asarray(transKwInAch)
    output['mcMechKwOutAch'] = np.asarray(mcMechKwOutAch)    
    output['auxInKw'] = np.asarray(auxInKw)
    output['mcElecKwInAch'] = np.asarray(mcElecKwInAch)
    output['mcMechKw4ForcedFc'] = np.asarray(mcMechKw4ForcedFc)
    output['canPowerAllElectrically'] = np.asarray(canPowerAllElectrically)
    output['fuelKg']=np.array(fuelKg)
    output['fuelKgAch']=np.array(fuelKgAch)
    output['mphAch']=np.array(mphAch)
    output['distMiles']=np.array(distMiles)
    
#     # aerodynamic dragg power
#     cycDragKw[i]
#     # rolling resistance power
#     cycRrKw[i]
#     # braking power
#     cycRegenBrakeKw[i]
#     # acceleration power
#     cycAccelKw[i]
#     # A/C compressor/fan air
#     rfgPcomp[i]/1e3
#     evapPfanair[i]/1e3
    output['cycDragKw']=np.array(cycDragKw)
    output['cycRrKw']=np.array(cycRrKw)
    output['cycRegenBrakeKw']=np.array(cycRegenBrakeKw)
    output['cycAccelKw']=np.array(cycAccelKw)
    c= np.linspace(0,int(1/dt)*len(cycRrKw),len(cycRrKw),dtype=int)
    output['rfgPcompKw']=np.array(rfgPcomp[[c]]/1e3)
    output['evapPfanairKw']=np.array(evapPfanair[c]/1e3)

    return output
