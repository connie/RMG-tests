import sys
import os.path
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from regressiveTest import ObservablesTestCase

molFracList=[{'CC': 0.05, '[Ar]': 0.95}]
Plist=[278643.8]
Tlist=range(1100,1300,100)
terminationTime = 5e-5
conditions=rt.generateAllConditions("IdealGasReactor", terminationTime, molFracList, Tlist, Plist)

majorSpeciesSmiles=['CC', 'C=C', 'C#C','C', '[CH3]', 'C[CH2]', '[H]', '[H][H]', 'C=[CH]', 'CC[CH2]', 'Ar']

A=rt.ObservablesTestCase("Ethane Pyrolysis (Minimal)", "regression/egA",
                    "regression/egB", conditions, exptData=None)
A.generateConditions("IdealGasReactor", terminationTime, molFracList, Tlist=Tlist, Plist=Plist)


print A
A.compare(True)
