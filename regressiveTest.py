import sys
import os.path
import numpy as np
from rmgpy.chemkin import loadChemkinFile
from rmgpy.tools.plot import *
from rmgpy.tools.cantera import *

def getNearestTime(timepoint, timeArray):
    """
    `timePoint`: the desired time point
    `timeArray`: the array of times in which to search for the nearest time to the one selected

    Returns a tuple containing (index, value of the time point closest to the timepoint desired)
    within the time array.
    """
    index = findNearest(timepoint, timeArray)
    return (index, timeArray(index))


def curvesSimilar(t1, y1, t2, y2, tol):
    """
    This function returns True if the two given curves are similar enough within tol. Otherwise returns False.

    t1: time/domain of standard curve we assume to be correct
    y1: values of standard curve, usually either temperature in (K) or log of a mol fraction
    t2: time/domain of test curve
    y2: values of test curve, usually either temperature in (K) or log of a mol fraction

    The test curve is first synchronized to the standard curve using geatNearestTime function. We then calculate the value of
    (y1-y2')^2/y1^2, giving us a normalized difference for every point. If the average value of these differences is less
    than tol, we say the curves are similar.

    We choose this criteria because it is compatible with step functions we expect to see in ignition systems.
    """
    # Make synchornized version of t2,y2 called t2sync,y2sync.
    t2sync=[]
    y2sync=[]
    for timepoint1 in t1:
        (index, timepoint2)=getNearestTime(timepoint1, t2sync)
        t2sync.append(timepoint2)
        y2sync.append(y2[index])

    # Get R^2 value equivalent:
    normalizedError=[(y1[x]-y2sync[x])**2/y1[x]**2 for x in range(len(y1))]/len(y1)

    if normalizedError > tol:
        return False
    else: 
        return True


class ObservablesTestCase:
    """
    We use this class to run regressive tests

    ======================= ==============================================================================================
    Attribute               Description
    ======================= ==============================================================================================
    Inputted attributs:
    'title'                 A string describing the test. For regressive tests, should be same as example's name.
    `oldDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the old model
    `newDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the new model
    `conditions`            A list of the :class: 'Condition' objects describing reaction conditions
    'exptData'              An array of GenericData objects

    Generated Attributes:
    'results'
    'oldDict'
    'newDict'
    ======================= ==============================================================================================


    """
    def __init__(self, title='', oldDir='', newDir='', exptData = []):
        self.title=title
        self.newDir=newDir
        self.oldDir=oldDir
        self.conditions=None
        self.exptData=exptData

        # load the species and reactions from each model
        oldSpeciesList, oldReactionList = loadChemkinFile(os.path.join(oldDir,'chem_annotated.inp'),
                                                          os.path.join(oldDir,'species_dictionary.txt'))

        newSpeciesList, newReactionList = loadChemkinFile(os.path.join(newDir,'chem_annotated.inp'),
                                                          os.path.join(newDir,'species_dictionary.txt'))
        self.oldSim = Cantera(speciesList = oldSpeciesList,
                              reactionList = oldReactionList,
                              outputDirectory = oldDir)
        self.newSim = Cantera(speciesList = newSpeciesList,
                              reactionList = newReactionList
                              outputDirectory = newDir)

        # load each chemkin file into the cantera model
        self.oldSim.loadChemkinModel(os.path.join(oldDir,'chem_annotated.inp'))
        self.newSim.loadChemkinModel(os.path.join(newDir,'chem_annotated.inp'))

    def __str__(self):
        """
        Return a string representation of this test case, using its title'.
        """
        return 'Observables Test Case: {0}'.format(self.title)

    def generateConditions(reactorType, reactionTime, molFracList, Tlist, Plist):
        """
        Creates a list of conditions from from the lists provided. 
        
        `reactorType`: a string indicating the Cantera reactor type
        `reactionTime`: ScalarQuantity object for time
        `molFracList`: a list of dictionaries containing species smiles and their mole fraction values
        `Tlist`: ArrayQuantity object of temperatures
        `Plist`: ArrayQuantity object of pressures
        `Vlist`: ArrayQuantity object of volumes
        
        This saves all the reaction conditions into both the old and new cantera jobs.
        """
        # Store the conditions in the observables test case, for bookkeeping
        self.conditions = generateConditions(reactorType, reactionTime, molFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)

        # Map the mole fractions dictionaries to species objects from the old and new models
        oldMolFracList = []
        newMolFracList = []

        for molFracCondition in molFracList:
            oldCondition = {}
            newCondition = {} 
            oldSpeciesDict = getRMGSpeciesFromSMILES(molFracCondition.keys(), self.oldSim.speciesList)
            newSpeciesDict = getRMGSpeciesFromSMILES(molFracCondition.keys(), self.newSim.speciesList)
            for smiles, molfrac in enumerate(molFracCondition):
                if oldSpeciesDict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the old model!'.format(smiles))
                if newSpeciesDict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the new model!'.format(smiles))

                oldCondition[oldSpeciesDict[smiles]] = molfrac
                newCondition[newSpeciesDict[smiles]] = molfrac
            oldMolFracList.append(oldCondition)
            oldMolFracList.append(newCondition)
        
        # Generate the conditions in each simulation
        self.oldSim.generateConditions(reactorType, reactionTime, oldMolFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)
        self.newSim.generateConditions(reactorType, reactionTime, newMolFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)

    def compare(self, plot=False):
        """
        Compare the old and new model
        `plot`: if set to True, it will comparison plots of the two models comparing their species.

        Returns a list of variables failed in a list of tuples in the format:
        
        (CanteraCondition, variable label, variableOld, variableNew)

        """
        oldConditionData, newConditionData = self.runSimulations()
        
        #For now make print statements about each none matching data
        conditionsBroken=[]
        variablesFailed=[]
        for i in range(len(oldConditionData)):
            timeOld, dataListOld = oldConditionData[i]
            timeNew, dataListNew = newConditionData[i]

            for variableOld, variableNew in zip(dataListOld, dataListNew):
                if not curvesSimilar(timeOld.data, variableOld.data, timeNew.data, variableNew.data, 0.05):
                    if i not in conditionsBroken: conditionsBroken.append(i)

                    if variableOld.species:
                        label = variableOld.species.molecule[0].toSMILES()
                        print "Species profile for {0} does not match between the old model ({1}) and \
                         the new model ({2}) in condition {3:d}.".format(variableOld.species.molecule[0].toSMILES(),
                                                                         variableOld.label, 
                                                                         variableNew.label,
                                                                         i+1)
                    else:
                        label = variableOld.label
                        print "{0} does not match between the old model and \
                         the new model in condition {1:d}.".format(variableOld.label, i+1)

                    variablesFailed.append((self.conditions[i], label, variableOld, variableNew))
        print ''
        print 'The following reaction conditions were broken:'
        print ''
        for index in conditionsBroken:
            print "Condition {0:d}:"
            print str(self.conditions[index])
            print ''

        if plot:
            # Ignore Inerts
            inertList = ['[Ar]','[He]','[N#N]','[Ne]']
            for i in range(len(oldConditionData)):
                time, dataList = oldConditionData[i]
                speciesData = [data for data in dataList if data.species.molecule[0].toSMILES() not in inertList]
                oldSpeciesPlot = SimulationPlot(xVar=time, yVar=speciesData, ylabel='Mole Fraction')

                time, dataList = newConditionData[i]
                speciesData = [data for data in dataList if data.species.molecule[0].toSMILES() not in inertList]
                newSpeciesPlot = SimulationPlot(xVar=time, yVar=speciesData, ylabel='Mole Fraction')

                # Name after the index of the condition
                # though it may be better to name it after the actual conditions in T, P, etc
                oldSpeciesPlot.comparePlot(newSpeciesPlot,filename='simulation_condition_{0}.png'.format(i+1))

        return variablesFailed

    def runSimulations(self):
        """
        Run a selection of conditions in Cantera and return
        generic data objects containing the time, pressure, temperature,
        and mole fractions from the simulations.

        Returns (oldConditionData, newConditionData)
        where conditionData is a list of of tuples: (time, dataList) for each condition in the same order as conditions
        time is a GenericData object which gives the time domain for each profile
        dataList is a list of GenericData objects for the temperature, profile, and mole fraction of major species
        """
        oldConditionData = self.oldSim.simulate()
        newConditionData = self.newSim.simulate()
        return (oldConditionData, newConditionData)
