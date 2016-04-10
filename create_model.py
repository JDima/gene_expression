__author__ = 'JDima'

tfs = ["hb", "Kr", "gt", "kni", "bcd", "cad", "tll", "hkb"]

from operator import or_

def getSite(core, tf, site):
    return "site_{0}_core_{1}_tf_{2}".format(site, core, tf)

def getBinded(core, tf, site):
    return "N_binded_core_{0}_tf_{1}_site_{2}".format(str(core), tf, str(site))

def getFree(core, tf, site):
    return "N_free_core_{0}_tf_{1}_site_{2}".format(str(core), tf, str(site))

def getOnState(core, tf, site):
    return "OnState_core_{0}_tf_{1}_site_{2}".format(str(core), tf, str(site))

def getOffState(core, tf, site):
    return "OffState_core_{0}_tf_{1}_site_{2}".format(str(core), tf, str(site))

class Model:
    def __init__(self, count_сore, transciption_factors, tf_probs):
        self.count_core = count_сore
        self.transcipton_factors = transciption_factors
        self.tf_probs = tf_probs

    def create_model(self):
        reactions = "# Reactions\n"
        species = set()

        br = BindingReaction(self.count_core, self.transcipton_factors, self.tf_probs)
        reacs, binded_species, sites, onStates = br.getReactions()
        reactions += reacs

        unbr = UnbindingReaction(self.count_core, self.transcipton_factors, self.tf_probs)
        reacs, unbinded_species, sites, offStates = unbr.getReactions()
        reactions += reacs

        trans_start = TranscriptionStartReaction(sites, 0.05)
        reactions += trans_start.getReactions()

        tr = TranslationReaction()
        reactions += tr.getReactions()

        degrad_mRNA = DegradationRNAReaction()
        reactions += degrad_mRNA.getReactions()

        degrad_prot = DegradationProteinReaction()
        reactions += degrad_prot.getReactions()

        species = unbinded_species

        reactions += "# InitVar\n"

        for site in sites:
            reactions += site + " = 1.0\n"

        for state in onStates:
            reactions += state + " = 0.0\n"

        for state in offStates:
            reactions += state + " = 1.0\n"

        for specie in species:
            if "free" in specie:
                reactions += specie + " = 67.975\n"
            else:
                reactions += specie + " = 0.0\n"

        reactions += "mRNA = 0.0\n"
        reactions += "Nprotein = 0.0\n\n"

        reactions += "# Parameter initialisations\nTranslationConst = 0.6\nRdegradation = 0.0002\nProtDegradation = 0.0002"

        return reactions

class BindingReaction:
    def __init__(self, count_сore, transciption_factors, tf_probs):
        self.count_core = count_сore
        self.transcipton_factors = transciption_factors
        self.tf_probs = tf_probs
        self.tf_power = 0.3

    def getBindingReaction(self, core, site, tf_power, tf, prob):
        freeN = getFree(core, tf, site)
        bindedN = getBinded(core, tf, site)
        siteVar = getSite(core, tf, site)
        onState = getOnState(core, tf, site)
        offState = getOffState(core, tf, site)
        format_str = ""
        if tf == "bcd" or tf == "cad":
            format_str = "binding_core_{0}_tf_{1}_site_{2}:\n\t{3} + {5}{6} + {10} > {4} + {9}\n\t{7} * {8} * {10}\n\n"
        else:
            format_str = "binding_core_{0}_tf_{1}_site_{2}:\n\t{3} + {10}> {4} + {9} + {5}{6}\n\t{7} * {8} * {10}\n\n"
        return format_str.format(
                    str(core), tf,  str(site), freeN, bindedN,
                    "{" + str(tf_power) + "}", siteVar, freeN, str(prob), onState, offState), freeN, bindedN, siteVar, onState

    def getReactions(self):
        reactions = ""
        species = []
        sites = []
        states = []

        for core in range(self.count_core):
            for site, tf_prob in enumerate(self.tf_probs):
                reaction, freeN, bindedN, siteVar, onState = self.getBindingReaction(core, site, self.tf_power, tf_prob[0], tf_prob[1])
                reactions += reaction

                species.append(freeN)
                species.append(bindedN)
                sites.append(siteVar)
                states.append(onState)
            reactions += "\n"
        return reactions, species, sites, states


class UnbindingReaction:
    def __init__(self, count_сore, transciption_factors, tf_probs):
        self.count_core = count_сore
        self.transcipton_factors = transciption_factors
        self.tf_probs = tf_probs
        self.tf_power = 0.1

    def getUnbindingReaction(self, core, site, tf_power, tf, prob):
        freeN = getFree(core, tf, site)
        bindedN = getBinded(core, tf, site)
        siteVar = getSite(core, tf, site)
        onState = getOnState(core, tf, site)
        offState = getOffState(core, tf, site)
        format_str = ""
        if tf == "bcd" or tf == "cad":
            format_str = "unbinding_core_{0}_tf_{1}_site_{2}:\n\t{3} + {10} > {4} + {5}{6} + {9}\n\t{7} * {8} * {10}\n\n"
        else:
            format_str = "unbinding_core_{0}_tf_{1}_site_{2}:\n\t{3} + {5}{6} + {10} > {4} + {9}\n\t{7} * {8} * {10}\n\n"

        return format_str.format(
                    str(core), tf,  str(site), bindedN, freeN,
                    "{" + str(tf_power) + "}", siteVar, bindedN, str(1 - prob), offState, onState), freeN, bindedN, siteVar, offState


    def getReactions(self):
        reactions = ""
        species = []
        sites = []
        states = []

        for core in range(self.count_core):
            for site, tf_prob in enumerate(self.tf_probs):
                reaction, freeN, bindedN, siteVar, offState = self.getUnbindingReaction(core, site, self.tf_power, tf_prob[0], tf_prob[1])
                reactions += reaction

                species.append(freeN)
                species.append(bindedN)
                sites.append(siteVar)
                states.append(offState)
            reactions += "\n"
        return reactions, species, sites, states

class TranscriptionStartReaction:
    def __init__(self, sites, defParam):
        self.sites = sites
        self.defParam = defParam

    def getReactions(self):
        prob = "*".join(self.sites)
        reaction =  "transcription_start:\n\t$pool > mRNA\n\t{0}*{1}\n\n".format(self.defParam, prob)
        return reaction

class TranslationReaction:
    def getReactions(self):
        return "translation:\n\tmRNA > Nprotein\n\tmRNA*TranslationConst\n\n"

class DegradationRNAReaction:
    def getReactions(self):
        return "degradation_mRNA:\n\tmRNA > $pool\n\tmRNA*Rdegradation\n\n"

class DegradationProteinReaction:
    def getReactions(self):
        return "degradation_protein:\n\tNprotein > $pool\n\tNprotein*ProtDegradation\n\n"

class DiffusionReaction:
    def getReactions(self):
        return "degradation_protein:\n\tNprotein > $pool\n\tNprotein*ProtDegradation\n\n"

def prob_tf(file):
    tf_prob = []
    for line in open(file):
        _, _, tf, _, _, prob,  = line.split()
        tf_prob.append((tfs[int(tf)], float(prob)))
    return tf_prob

tf_prob = prob_tf("sitesDR.ann")
tf_prob = tf_prob[:2]
model = Model(2, tfs, tf_prob)
print(model.create_model())

with open("expression2.psc","w") as out:
    print(model.create_model(),file=out)


import stochpy
import matplotlib.gridspec as gridspec

sim_end = 50

smod = stochpy.SSA()
smod.Model("expression2.psc", dir="")



smod.DoStochSim(end = sim_end, mode = 'time', method='TauLeaping')

smod.PlotSpeciesTimeSeries()

gs = gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[0.3,1,0.3,1])

ax1 = stochpy.plt.subplot(gs[0])

smod.PlotSpeciesTimeSeries(species2plot ='Nprotein',xlabel ='',ylabel ='')
smod.plot.ResetPlotnum()
smod.PlotSpeciesTimeSeries(species2plot ='Nprotein',colors = ['r'])
stochpy.plt.xlim([0,sim_end])
stochpy.plt.legend('',frameon=False)
stochpy.plt.xticks([])
stochpy.plt.title('')
stochpy.plt.ylabel('Nprotein')
stochpy.plt.yticks([0,15,30])
stochpy.plt.text(101,27,'D',fontsize = 14)


ax2 = stochpy.plt.subplot(gs[1])
smod.plot.ResetPlotnum()
smod.PlotSpeciesTimeSeries(species2plot ='mRNA',colors = ['#32CD32'])
stochpy.plt.xlim([0, sim_end])
stochpy.plt.legend('',frameon=False)
stochpy.plt.xticks([])
stochpy.plt.title('')
stochpy.plt.xlabel('')
stochpy.plt.ylabel('mRNA')
stochpy.plt.yticks([0,100,200,300,400])
stochpy.plt.text(101,27,'B',fontsize = 14)
stochpy.plt.savefig('stochpy_test.pdf')