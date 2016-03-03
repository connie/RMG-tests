import sys
import os.path
from rmgpy.tools.regressiveTest import ObservablesTestCase
from rmgpy.quantity import Quantity


# Conditions
reactorType = 'IdealGasReactor'
molFracList=[{'CC': 0.05, '[Ar]': 0.95}]
Plist=Quantity([278643.8],'Pa')
Tlist=Quantity([1100,1200,1300],'K')
terminationTime = Quantity(5e-5,'s')

# Create observables test case and compare the old and new models

minimal = ObservablesTestCase(title = 'Ethane Pyrolysis',
                              oldDir = 'old',
                              newDir = 'new')

minimal.generateConditions(reactorType = reactorType,
                           reactionTime = terminationTime,
                           molFracList = molFracList,
                           Tlist = Tlist,
                           Plist = Plist)

print minimal
minimal.compare(plot = True)
