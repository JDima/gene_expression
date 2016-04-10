__author__ = 'JDima'

import stochpy

smod = stochpy.SSA()
smod.Model("expression2.psc", dir="")
smod.DoStochSim(end = 50,mode = 'time')
smod.PlotSpeciesTimeSeries()