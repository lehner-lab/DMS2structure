
xplor.requireVersion("2.43.4")


#
# slow cooling protocol in torsion angle space for protein G. Uses 
# NOE, J-coupling restraints.
#
# this script performs annealing from an extended structure
#
# CDS 2009/07/24
#

# this checks for typos on the command-line. User-customized arguments can
# also be specified.
#
xplor.parseArguments()


# filename for output structures. This string must contain the STRUCTURE
# literal so that each calculated structure has a unique name. The SCRIPT
# literal is replaced by this filename (or stdin if redirected using <),
# but it is optional.
#
outFilename = "XPLOR_output/SCRIPT_STRUCTURE.pdb"
numberOfStructures=100   #usually you want to create at least 20 

# protocol module has many high-level helper functions.
#
import protocol

protocol.initRandomSeed(1603)   #set random seed - by time

command = xplor.command

# generate PSF data from sequence and initialize the correct parameters.
#
from psfGen import seqToPSF
seqToPSF('protG.seq') ## modify ##


# generate random extended initial structure with correct covalent geometry
#
protocol.genExtendedStructure()

#
# a PotList contains a list of potential terms. This is used to specify which
# terms are active during refinement.
#
from potList import PotList
potList = PotList()

# parameters to ramp up during the simulated annealing protocol
#
from simulationTools import MultRamp, StaticRamp, InitialParams, IVMAction

rampedParams=[]
highTempParams=[]

# IVM setup
#   the IVM is used for performing dynamics and minimization in torsion-angle
#   space, and in Cartesian space.
#
from ivm import IVM
dyn  = IVM()
minc = IVM() # minc used for final cartesian minimization

# initialize ivm topology for torsion-angle dynamics


#
# 
#


# compare atomic Cartesian rmsd with a reference structure
#  backbone and heavy atom RMSDs will be printed in the output
#  structure files
#
# from posDiffPotTools import create_PosDiffPot
# refRMSD = create_PosDiffPot("refRMSD","name CA or name C or name N",
#                            pdbFile='g_xray.pdb', ## modify ##
#                            cmpSel="not name H*")

# set up NOE potential
noe=PotList('noe')
potList.append(noe)
from noePotTools import create_NOEPot
for (name,scale,file) in [('all',1,"NOE_contacts55rand_perfect.tbl"), ## modify ##
                          #add entries for additional tables
                          ]:
    pot = create_NOEPot(name,file)
    ## use soft PotType##
    # pot.setPotType("soft") # if you think there may be bad NOEs 
    pot.setAveType("shortest") ## use shortest distance between atom selections (side chain heavy atoms)
    pot.setScale(scale)
    noe.append(pot)
rampedParams.append( MultRamp(2,30, "noe.setScale( VALUE )") )


## no J couplings ##
# set up J coupling - with Karplus coefficients
##from jCoupPotTools import create_JCoupPot
##jCoup = create_JCoupPot("jcoup","jna_coup.tbl",
##                        A=6.98,B=-1.38,C=1.72,phase=-60.0)
##potList.append(jCoup)

# Set up dihedral angles
from xplorPot import XplorPot
dihedralRestraintFilename="dihe_ss_1pga.tbl" ## modify, and check that if works if dihedral angles are incomplete ##
protocol.initDihedrals(dihedralRestraintFilename,
                       #useDefaults=False  # by default, symmetric sidechain
                                           # restraints are included
                       )
potList.append( XplorPot('CDIH') )
highTempParams.append( StaticRamp("potList['CDIH'].setScale(10)") )
rampedParams.append( StaticRamp("potList['CDIH'].setScale(200)") )
# set custom values of threshold values for violation calculation
#
potList['CDIH'].setThreshold( 5 )



# gyration volume term 
#
from gyrPotTools import create_GyrPot
gyr = create_GyrPot("Vgyr",
                    "resid 1:56") # selection should exclude disordered tails
potList.append(gyr)
rampedParams.append( MultRamp(.002,1,"gyr.setScale(VALUE)") )


## what is this exactly doing?? ##
# hbdb - hbond database-based term
#
protocol.initHBDB()
potList.append( XplorPot('HBDB') )

#New torsion angle database potential
#
from torsionDBPotTools import create_TorsionDBPot
torsionDB = create_TorsionDBPot('torsionDB')
potList.append( torsionDB )
rampedParams.append( MultRamp(.002,2,"torsionDB.setScale(VALUE)") )

#
# setup parameters for atom-atom repulsive term. (van der Waals-like term)
#
from repelPotTools import create_RepelPot,initRepel
repel = create_RepelPot('repel')
potList.append(repel)
rampedParams.append( StaticRamp("initRepel(repel,use14=False)") )
rampedParams.append( MultRamp(.004,4,  "repel.setScale( VALUE)") )
# nonbonded interaction only between CA atoms
highTempParams.append( StaticRamp("""initRepel(repel,
                                               use14=True,
                                               scale=0.004,
                                               repel=1.2,
                                               moveTol=45,
                                               interactingAtoms='name CA'
                                               )""") )

# Selected 1-4 interactions.
import torsionDBPotTools
repel14 = torsionDBPotTools.create_Terminal14Pot('repel14')
potList.append(repel14)
highTempParams.append(StaticRamp("repel14.setScale(0)"))
rampedParams.append(MultRamp(0.004, 4, "repel14.setScale(VALUE)"))


potList.append( XplorPot("BOND") )
potList.append( XplorPot("ANGL") )
potList['ANGL'].setThreshold( 5 )
rampedParams.append( MultRamp(0.4,1,"potList['ANGL'].setScale(VALUE)") )
potList.append( XplorPot("IMPR") )
potList['IMPR'].setThreshold( 5 )
rampedParams.append( MultRamp(0.1,1,"potList['IMPR'].setScale(VALUE)") )
      


# Give atoms uniform weights, configure bath/molecule friction coeff.
#
protocol.massSetup()


# IVM setup
#   the IVM is used for performing dynamics and minimization in torsion-angle
#   space, and in Cartesian space.
#
from ivm import IVM
dyn = IVM()

# initialize ivm topology for torsion-angle dynamics
#
protocol.torsionTopology(dyn)

# minc used for final cartesian minimization
#
minc = IVM()
protocol.initMinimize(minc)

protocol.cartesianTopology(minc)



# object which performs simulated annealing
#
from simulationTools import AnnealIVM
init_t  = 3500.     # Need high temp and slow annealing to converge
cool = AnnealIVM(initTemp =init_t,
                 finalTemp=25,
                 tempStep =12.5,
                 ivm=dyn,
                 rampedParams = rampedParams)


def calcOneStructure(loopInfo):
    """ this function calculates a single structure, performs analysis on the
    structure, and then writes out a pdb file, with remarks.
    """

    # generate a new structure with randomized torsion angles
    #
    from monteCarlo import randomizeTorsions
    randomizeTorsions(dyn)

    # set torsion angles from restraints
    #
    from torsionTools import setTorsionsFromTable
    setTorsionsFromTable(dihedralRestraintFilename) ## does this work if incomplete? ##
    protocol.fixupCovalentGeom(maxIters=100,useVDW=1)
    protocol.writePDB(loopInfo.filename()+".init")


    # initialize parameters for high temp dynamics.
    InitialParams( rampedParams )
    # high-temp dynamics setup - only need to specify parameters which
    #   differfrom initial values in rampedParams
    InitialParams( highTempParams )

    # high temp dynamics
    #
    protocol.initDynamics(dyn,
                          potList=potList, # potential terms to use
                          bathTemp=init_t,
                          initVelocities=1,
                          finalTime=100,   # stops at 100ps or 1000 steps
                          numSteps=1000,   # whichever comes first
                          printInterval=100)

    dyn.setETolerance( init_t/100 )  #used to det. stepsize. default: t/1000 
    dyn.run()

    # initialize integrator for simulated annealing
    #
    protocol.initDynamics(dyn,
                          potList=potList,
                          numSteps=100,       #at each temp: 100 steps or
                          finalTime=.2 ,       # .2ps, whichever is less
                          printInterval=100)

    # perform simulated annealing
    #
    cool.run()
              
              
    # final torsion angle minimization
    #
    protocol.initMinimize(dyn,
                          printInterval=50)
    dyn.run()

    # final all- atomic degrees of freedom minimization
    #
    protocol.initMinimize(minc,
                          potList=potList,
                          dEPred=10)
    minc.run()

    #do analysis and write structure when function returns
    pass


## modify ##
from simulationTools import StructureLoop, FinalParams
StructureLoop(numStructures=numberOfStructures,
              doWriteStructures=True,
              pdbTemplate=outFilename,
              structLoopAction=calcOneStructure,
              genViolationStats=True,
              averageTopFraction=0.25, #report stats on best 25% of structs
              averageSortPots=[potList['torsionDB']],
              averageContext=FinalParams(rampedParams),
              #averageCrossTerms=refRMSD,
              averageFilename="XPLOR_output/SCRIPT_ave.pdb",
              averagePotList=potList).run()
