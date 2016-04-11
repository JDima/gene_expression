__author__ = 'JDima'

from lxml import etree

tfs = ["hb", "Kr", "gt", "kni", "bcd", "cad", "tll", "hkb"]

def getSite(core, tf, site):
    # return "site{0}core{1}".format(site, core)
    return "site{0}core{1}tf{2}".format(site, core, tf)

def getBinded(core, tf, site = ""):
    return "NbindedCore{0}tf{1}".format(str(core), tf)
    # return "NbindedCore{0}tf{1}site{2}".format(str(core), tf, str(site))

def getFree(core, tf, site = ""):
    return "NfreeCore{0}tf{1}".format(str(core), tf)
    # return "NfreeCore{0}tf{1}site{2}".format(str(core), tf, str(site))

def getOnState(core, tf, site):
    return "OnStateCore{0}tf{1}site{2}".format(str(core), tf, str(site))
    # return "OnStateCore{0}site{1}".format(str(core), str(site))

def getOffState(core, tf, site):
    return "OffStateCore{0}tf{1}site{2}".format(str(core), tf, str(site))
    # return "OffStateCore{0}site{1}".format(str(core), str(site))


class ModelPSC:

    def add_child(self, root, child_name, text = None):
        child = etree.Element(child_name)
        if text != None:
            child.text = str(text)
        root.append(child)
        return child


    def __init__(self, count_сore, transciption_factors, tf_probs, diffuse, h_dist):
        self.count_core = count_сore
        self.transcipton_factors = transciption_factors
        self.tf_probs = tf_probs

        self.model = etree.Element('Model')
        self.add_child(self.model, 'Description', 'Gene expression')
        self.reaction_id = 1

        self.diffuse = diffuse
        self.h_dist = h_dist
        self.trans_start = 0.1

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
        species = set()
        sites = set()
        states = set()
        for core in range(self.count_core):
            for site, tf_prob in enumerate(self.tf_probs):
                freeN = getFree(core, tf_prob[0], site)
                bindedN = getBinded(core, tf_prob[0], site)
                siteVar = getSite(core, tf_prob[0], site)
                onState = getOnState(core, tf_prob[0], site)
                offState = getOffState(core, tf_prob[0], site)

                desc = "unbinding_core_{0}_tf_{1}_site_{2}".format(str(core), tf_prob[0],  str(site))
                rate = "{0} * {1} * {2}".format(str(1 - tf_prob[1]), onState, bindedN)

                if tf_prob[0] == "bcd" or tf_prob[0] == "cad":
                    reactants = { bindedN : '1', onState : '1'}
                    products = { freeN :  '1', offState : '1', siteVar :  '1'}
                else:
                    reactants = { bindedN : '1', onState : '1', siteVar :  '1'}
                    products = { freeN :  '1', offState : '1'}

                self.add_reaction(root, self.reaction_id, desc, rate, reactants, products)

                species.add(freeN)
                species.add(bindedN)
                sites.add(siteVar)
                states.add(offState)
                states.add(onState)
        return species, sites, states

    def add_transription_start(self, root, sites):
        desc = "Transcription start"

        for core in range(self.count_core):
            site_by_core = []
            for site in sites:
                if 'core' + str(core) in site:
                    site_by_core.append(site)
            site_by_core.append(str(self.trans_start))
            rate = " * ".join(site_by_core)
            self.add_reaction(root, self.reaction_id, desc, rate, {}, {'mRNA' + str(core) : '1'})

    def add_translation_reaction(self, root):
        desc = "Translation"
        for core in range(self.count_core):
            self.add_reaction(root, self.reaction_id, desc, 'mRNA' + str(core) + ' * translation', {'mRNA' + str(core) : '1'}, {'Nprotein' + str(core) : '1'})

    def add_degradation_mRNA_reaction(self, root):
        desc = "Degradation mRNA"
        for core in range(self.count_core):
            self.add_reaction(root, self.reaction_id, desc, 'mRNA' + str(core) + ' * mRNA_degrad', {'mRNA' + str(core) : '1'}, {})

    def add_degradation_protein_reaction(self, root):
        desc = "Degradation protein"
        for core in range(self.count_core):
            self.add_reaction(root, self.reaction_id, desc, 'Nprotein' + str(core) + ' * Protein_degrad', {'Nprotein' + str(core) : '1'}, {})




    def add_diff_reaction(self, root, desc, fromP, toP):
        dconst = str(self.diffuse / self.h_dist)
        self.add_reaction(root, self.reaction_id,
                          desc, dconst + ' * (' + fromP + ' - ' + toP + ')' + ' * (' + fromP + ' - ' + toP + ')',
                          {fromP : '1'}, {toP : '1'})

    def add_diffusion_reaction(self, root):
        if self.count_core == 1:
          return

        dconst = self.diffuse / self.h_dist

        tfs = { tf for tf, prob in self.tf_probs}

        for tf in tfs:
            freeLeft = getFree(0, tf)
            freeLeftNeighbour = getFree(1, tf)
            self.add_diff_reaction(root,  "diffusionToRight" + tf + "core" + str(0), freeLeft, freeLeftNeighbour)

            freeRight = getFree(self.count_core - 1, tf)
            freeRightNeighbour = getFree(self.count_core - 2, tf)
            self.add_diff_reaction(root,  "diffusionToLeft" + tf + "core" + str(self.count_core - 1), freeRight, freeRightNeighbour)

        if self.count_core == 2:
            return

        for core in range(1, self.count_core - 1):
            for tf in tfs:
                freeN = getFree(core, tf)
                freeNLeft = getFree(core - 1, tf)
                freeNRight = getFree(core + 1, tf)

                self.add_diff_reaction(root,  "diffusionToRight" + tf + "core" + str(core), freeN, freeNRight)
                self.add_diff_reaction(root,  "diffusionToLeft" + tf + "core" + str(core), freeN, freeNLeft)



    def add_reactions(self):
        pl_root = etree.Element('ReactionsList')
        self.add_binding_reactions(pl_root)
        species, sites, states = self.add_unbinding_reactions(pl_root)
        self.add_transription_start(pl_root, sites)
        self.add_translation_reaction(pl_root)
        self.add_degradation_mRNA_reaction(pl_root)
        self.add_degradation_protein_reaction(pl_root)
        self.add_diffusion_reaction(pl_root)
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

        parametrs = {'translation' : 0.2, 'mRNA_degrad' : 0.2, 'Protein_degrad' : 0.2, 'Transcription_start' : self.trans_start}
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

        for core in range(self.count_core):
            self.add_specie(sl_root, 'mRNA' + str(core), 'mRNA' + str(core), 0.0)
            self.add_specie(sl_root, 'Nprotein' + str(core), 'Nprotein' + str(core), 0.0)

        return len(states) + len(species) + len(sites) + 2 * self.count_core, sl_root


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

tf_prob = tf_prob[:1]

model = ModelPSC(3, tfs, tf_prob, 0.002, 10000)
doc = model.create_model()


print(etree.tostring(doc, pretty_print=True))

f = open('model.xml','wb')
f.write(etree.tostring(doc, pretty_print=True))
f.close()