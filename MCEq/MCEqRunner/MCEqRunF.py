import MCEq.density_profiles as den_p
import numpy as np
import CRFluxModels as pm
from MCEq.core import MCEqRun
from mceq_config import config, mceq_config_without





#This function calculates a cord across the Earth.  LLo and LLhi are the minimum and maximum values
# of the parameter along the cord.  LToH converts values of L to atmospheric heights H
def LToHConv(b,Rs,H):
    LToH=0
    if b>1.:
        LToH=lambda(L):(pow(L*L + b*b*Rs*Rs, 0.5)-Rs)
        LLo=-Rs*pow(pow(1+H/Rs,2)-b*b,0.5)
        LHi=Rs*pow(pow(1+H/Rs,2)-b*b,0.5)
    else :
        LToH=lambda(L):(Rs* (pow(b*b+pow(pow(1-b*b,0.5)+L/Rs,2),0.5)-1) )
        LLo=0
        LHi=pow(pow(H+Rs,2)-b*b*Rs*Rs,0.5)-Rs*pow(1-b*b,0.5)
    return [LToH, LLo, LHi]


#Setup the MCEQ module with some position and density spectrum
def SetupMCEq(rawpositions,rawdensities,Settings={}):

    
    ad = den_p.GeneralizedTarget()
    densities=[]
    positions=[]
    if(rawpositions[0]>0):
        positions.append(0)
        densities.append(rawdensities[0])
    for i in range(0,len(rawdensities)):
        if not np.isnan(rawdensities[i]):
            densities.append(rawdensities[i])
            positions.append(rawpositions[i])
 #   ad.len_target=positions[-1]
    ad.set_length(positions[-1]-positions[0])
    for i in range(0, len(densities)-1):
        ad.add_material(positions[i],densities[i],'Layer'+str(i))
    ad.print_table()

    if Settings['PrimaryModel']=='HillasGaisser_H3a':
        PrimaryObject=pm.HillasGaisser2012
        PrimaryString='H3a'

    if Settings['PrimaryModel']=='HillasGaisser_H4a':
        PrimaryObject=pm.HillasGaisser2012
        PrimaryString='H4a'

    if Settings['PrimaryModel']=='GaisserHonda':
        PrimaryObject=pm.GaisserHonda
        PrimaryString=None

    if Settings['PrimaryModel']=='CombinedGHAndHG_H3a':
        PrimaryObject=pm.CombinedGHandHG
        PrimaryString="H3a"

    if Settings['PrimaryModel']=='CombinedGHAndHG_H4a':
        PrimaryObject=pm.CombinedGHandHG
        PrimaryString="H4a"

    if Settings['PrimaryModel']=='PolyGonato':
        PrimaryObject=pm.PolyGonato
        PrimaryString=None

    if Settings['PrimaryModel']=='Thunman':
        PrimaryObject=pm.Thunman
        PrimaryString=None

    if Settings['PrimaryModel']=='ZatsepinSokolskaya':
        PrimaryObject=pm.ZatsepinSokolskaya
        PrimaryString=None
 






    mceq_run = MCEqRun(
        interaction_model=Settings['HadronicModel'],

        primary_model=(PrimaryObject, PrimaryString),
        #Do not provide any default values to avoid unnecessary initilizations
        theta_deg=0.,
        density_model=('GeneralizedTarget', None),
        # Expand the rest of parameters
        **mceq_config_without(['density_model','integrator'])
        )
#    config['integrator'] = 'odepack'

    mceq_run.density_model=ad
    res_group=["mu+", "mu-", "pi_mu+", "pi_mu-", 
               "k_mu+", "k_mu-", "pr_mu+", "pr_mu-"]
    mceq_run.set_obs_particles(res_group)
    return [mceq_run, ad]


#Given an MCEq solution, get the array of position dependent fluxes to store
def GetPosArray(TheMceq, ListOfStrings, posgrid,Positions):
    ReturnArray=[]
    nparticles = lambda pname, g_idx: np.sum(
        TheMceq.get_solution(pname,mag=0,grid_idx=g_idx)*TheMceq.e_widths)
    ReturnArray.append(posgrid)
    ReturnArray.append(Positions)
    for part in ListOfStrings:
        ReturnArray.append([nparticles(part,i) for i in range(len(posgrid))])
    return ReturnArray


#Given an MCEq solution, get the array of energy dependent fluxes to store
def GetEnergyArray(TheMceq, ListOfStrings,mag=3):
    ReturnArray=[]
    ReturnArray.append(TheMceq.e_grid)
    for i in ListOfStrings:
        ReturnArray.append(TheMceq.get_solution(i,mag=mag))
    return ReturnArray

#Get IDs for each type of MCEq calculated flux
def MakeNameDicts():
    IDToName=dict()
    NameToID=dict()
    Parts=['conv_','pr_','total_','pi_','k_','obs_']
    Nus=['numu','antinumu','nutau','antinutau','nue','antinue','mu+','mu-']
    count=0
    for j in range(0, len(Nus)):
        for i in range(0, len(Parts)):
            Name=Parts[i]+Nus[j]
            NameToID[count]=Name
            IDToName[Name]=count
            count=count+1
    return [NameToID,IDToName]


#Setup density profiles to test
def SetupTestPoints(bvals,DensFunction,SolarRadius,ChromosHeight,Divs=25,BinSchemeIndex=1):
    TestPoints=[]
    eps=0.00001
    for b in bvals:
        if(b<1):
            LToH,LLo,LHi = LToHConv(b,SolarRadius,ChromosHeight)
            vars=[]
            vars=(pow(np.linspace(1-eps,eps,Divs),BinSchemeIndex))*(LHi-LLo) + LLo
            densities=DensFunction(LToH(vars))
            TestPoints.append([b,LHi-vars,densities])        
        else:
            LToH,LLo,LHi =LToHConv(b,SolarRadius,ChromosHeight)
            vars=np.linspace(LHi*(1-eps),LLo*(1-eps),Divs)
            densities=DensFunction(LToH(vars))
            TestPoints.append([b,LHi-vars,densities])
    return TestPoints

# Runs MCEq for a provided density profile and makes output files
def RunMCEq(Positions,Densities,Settings={}):
    #Run the simulation
    MakeFiles=Settings['MakeFiles']
    FileName=Settings['FileName']
    
    [mceq_run,ad]=SetupMCEq(Positions,Densities,Settings)
    x_grid = ad.s_h2X(Positions)
    
    mceq_run.solve(int_grid=x_grid)
    
    #Save the output arrays
    IDToName,NameToID = MakeNameDicts()
    ReturnPosArray=GetPosArray(mceq_run, IDToName.values(),x_grid,Positions)
    ReturnEnergyArray=GetEnergyArray(mceq_run, IDToName.values(),2)
    if(MakeFiles):
        np.savetxt(FileName+"_E.txt",ReturnEnergyArray)
        np.savetxt(FileName+"_X.txt",ReturnPosArray)
    return [mceq_run,ad]
