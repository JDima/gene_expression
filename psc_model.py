__author__ = 'JDima'

from lxml import etree

tfs = ["hb", "Kr", "gt", "kni", "bcd", "cad", "tll", "hkb"]

def getSite(core, tf, site):
    return "site{0}core{1}tf{2}".format(site, core, tf)

def getBinded(core, tf, site):
    return "NbindedCore{0}tf{1}site{2}".format(str(core), tf, str(site))

def getFree(core, tf, site):
    return "NfreeCore{0}tf{1}site{2}".format(str(core), tf, str(site))

def getOnState(core, tf, site):
    return "OnStateCore{0}tf{1}site{2}".format(str(core), tf, str(site))

def getOffState(core, tf, site):
    return "OffStateCore{0}tf{1}site{2}".format(str(core), tf, str(site))


class ModelPSC:

    def add_child(self, root, child_name, text = None):
        child = etree.Element(child_name)
        if text != None:
            child.text = str(text)
        root.append(child)
        return child


    def __init__(self, count_сore, transciption_factors, tf_probs):
        self.count_core = count_сore
        self.transcipton_factors = transciption_factors
        self.tf_probs = tf_probs

        self.model = etree.Element('Model')
        self.add_child(self.model, 'Description', 'Gene expression')
        self.reaction_id = 1

    def add_reaction(self, root, id, description, rate, reactants, products, type = 'customized'):
        reaction = self.add_child(root, 'Reaction')
        self.add_child(reaction, 'Id', id)
        self.add_child(reaction, 'Description', description)
        self.add_child(reaction, 'Type', type)
        self.add_child(reaction, 'PropensityFunction', rate)

        if reactants:
            ereactants = self.add_child(reaction, 'Reactants')
            ereactants.extend([etree.Element('SpeciesReference', {'id': reactant, 'stoichiometry' : reactants[reactant]}) for reactant in reactants.keys()])
            reaction.append(ereactants)

        if products:
            eproducts = self.add_child(reaction, 'Products')
            eproducts.extend([etree.Element('SpeciesReference', {'id': product, 'stoichiometry' : products[product]}) for product in products.keys()])
            reaction.append(eproducts)

        self.reaction_id += 1

        return reaction

    def add_binding_reactions(self, root):
        for core in range(self.count_core):
            for site, tf_prob in enumerate(self.tf_probs):
                freeN = getFree(core, tf_prob[0], site)
                bindedN = getBinded(core, tf_prob[0], site)
                siteVar = getSite(core, tf_prob[0], site)
                onState = getOnState(core, tf_prob[0], site)
                offState = getOffState(core, tf_prob[0], site)

                desc = "binding_core_{0}_tf_{1}_site_{2}".format(str(core), tf_prob[0],  str(site))
                rate = "{0} * {1} * {2}".format(str(tf_prob[1]), offState, freeN)

                if tf_prob[0] == "bcd" or tf_prob[0] == "cad":
                    reactants = { freeN : '1', siteVar :  '1', offState : '1'}
                    products = { bindedN :  '1', onState : '1'}
                else:
                    reactants = { freeN : '1', offState : '1'}
                    products = { bindedN :  '1', onState : '1', siteVar :  '1'}

                self.add_reaction(root, self.reaction_id, desc, rate, reactants, products)

        return freeN, bindedN, siteVar, onState, offState


    def add_unbinding_reactions(self, root):
        species = []
        sites = []
        states = []
        for core in range(self.count_core):
            for site, tf_prob in enumerate(self.tf_probs):
                freeN = getFree(core, tf_prob[0], site)
                bindedN = getBinded(core, tf_prob[0], site)
                siteVar = getSite(core, tf_prob[0], site)
                onState = getOnState(core, tf_prob[0], site)
                offState = getOffState(core, tf_prob[0], site)

                desc = "unbinding_core_{0}_tf_{1}_site_{2}".format(str(core), tf_prob[0],  str(site))
                rate = "{0} * {1} * {2}".format(str(tf_prob[1]), onState, bindedN)

                if tf_prob[0] == "bcd" or tf_prob[0] == "cad":
                    reactants = { bindedN : '1', onState : '1'}
                    products = { freeN :  '1', offState : '1', siteVar :  '1'}
                else:
                    reactants = { bindedN : '1', onState : '1', siteVar :  '1'}
                    products = { freeN :  '1', offState : '1'}

                self.add_reaction(root, self.reaction_id, desc, rate, reactants, products)

                species.append(freeN)
                species.append(bindedN)
                sites.append(siteVar)
                states.append(offState)
                states.append(onState)
        return species, sites, states

    def add_transription_start(self, root, sites):
        rate = " * ".join(sites)
        desc = "Transcription start"

        self.add_reaction(root, self.reaction_id, desc, rate, {}, {'mRNA' : '1'})

    def add_translation_reaction(self, root):
        desc = "Translation"
        self.add_reaction(root, self.reaction_id, desc, 'mRNA * translation', {'mRNA' : '1'}, {'Nprotein' : '1'})

    def add_degradation_mRNA_reaction(self, root):
        desc = "Degradation mRNA"
        self.add_reaction(root, self.reaction_id, desc, 'mRNA * mRNA_degrad', {'mRNA' : '1'}, {})

    def add_degradation_protein_reaction(self, root):
        desc = "Degradation protein"
        self.add_reaction(root, self.reaction_id, desc, 'Nprotein * Protein_degrad', {'Nprotein' : '1'}, {})

    def add_reactions(self):
        pl_root = etree.Element('ReactionsList')
        self.add_binding_reactions(pl_root)
        species, sites, states = self.add_unbinding_reactions(pl_root)
        self.add_transription_start(pl_root, sites)
        self.add_translation_reaction(pl_root)
        self.add_degradation_mRNA_reaction(pl_root)
        self.add_degradation_protein_reaction(pl_root)
        return species, sites, states, pl_root

    def add_parametr(self, pl_root, key, value):
        parameter = self.add_child(pl_root, 'Parameter')
        id = etree.Element('id')
        id.text = key
        expression = etree.Element('Expression')
        expression.text = str(value)
        parameter.extend([id, expression])

    def add_parameters_list(self, species, sites, states):
        pl_root = etree.Element('ParametersList')

        parametrs = {'translation' : 0.2, 'mRNA_degrad' : 0.2, 'Protein_degrad' : 0.2}
        for key in parametrs.keys():
            self.add_parametr(pl_root, key, parametrs[key])

        return pl_root

    def add_specie(self, root, id, desc, init):
        sroot = self.add_child(root, 'Species')
        self.add_child(sroot, 'Id', id)
        self.add_child(sroot, 'Description', desc)
        self.add_child(sroot, 'InitialPopulation', init)

    def add_species_list(self, species, sites, states):
        sl_root = etree.Element('SpeciesList')

        for specie in species:
            self.add_specie(sl_root, specie, specie, 68 if 'free' in specie else 0.0)

        for state in states:
            self.add_specie(sl_root, state, state, 0.0 if 'On' in state else 1)

        for site in sites:
            self.add_specie(sl_root, site, site, 2)

        self.add_specie(sl_root, 'mRNA', 'mRNA', 0.0)
        self.add_specie(sl_root, 'Nprotein', 'Nprotein', 0.0)

        return len(states) + len(species) + len(sites) + 2, sl_root


    def create_model(self):
        species, sites, states, reac_root = self.add_reactions()
        pl_root = self.add_parameters_list(species, sites, states)
        countSpecies, s_root = self.add_species_list(species, sites, states)

        self.add_child(self.model, 'NumberOfReactions', self.reaction_id - 1)
        self.add_child(self.model, 'NumberOfSpecies', countSpecies)

        self.model.append(pl_root)
        self.model.append(reac_root)
        self.model.append(s_root)

        return self.model


def prob_tf(file):
    tf_prob = []
    for line in open(file):
        _, _, tf, _, _, prob,  = line.split()
        tf_prob.append((tfs[int(tf)], float(prob)))
    return tf_prob

tf_prob = prob_tf("sitesDR.ann")

tf_prob = tf_prob[:2]

model = ModelPSC(2, tfs, tf_prob)
doc = model.create_model()


print(etree.tostring(doc, pretty_print=True))

f = open('model.xml','wb')
f.write(etree.tostring(doc, pretty_print=True))
f.close()