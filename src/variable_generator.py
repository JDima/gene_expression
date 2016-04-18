__author__ = 'JDima'


def getSite(core, tf, site):
    return "count_repr".format(site, core, tf)

def getCounts():
    return "crepres", "cactiv"

def getBinded(core, tf, site=""):
    return "Nbinded_core_{0}_{1}".format(str(core), tf)


def getFree(core, tf, site=""):
    return "Nfree_core_{0}_{1}".format(str(core), tf)


def getOnState(core, tf, site):
    return "OnState_core_{0}_site_{2}_{1}".format(str(core), tf, str(site))


def getOffState(core, tf, site):
    return "OffState_core_{0}_site_{2}_{1}".format(str(core), tf, str(site))

def getVariables(core, tf, site):
    return getFree(core, tf), \
           getBinded(core, tf), \
           getOnState(core, tf, site),\
           getOffState(core, tf, site)
