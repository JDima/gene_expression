import math
from src.adder import adder
from src.variable_generator import *
from lxml import etree

tf_count = {
    "bcd": [ 27.407, 21.598, 16.408, 12.328, 9.036, 7.160, 5.754, 3.495, 1.596, 0.160],
    "cad": [ 22.903, 29.267, 40.315, 49.829, 52.643, 53.950, 56.071, 66.993, 79.425, 56.242],
    "tll": [1.738, 0.669, 1.113, 1.000, 0.944, 1.356, 1.251, 2.976, 22.632, 99.587],
    "hkb": [7.476, 7.567, 7.735, 8.101, 7.519, 6.138, 5.543, 4.688, 6.486, 16.806],
    "hb": [164.428, 154.157, 25.905, 3.606, 0.785, 0.009, -0.002, 30.725, 111.775, 74.768],
    "Kr": [16.052, 131.507, 200.316, 177.865, 77.429, 10.447, 1.579, 0.394, 0.483, 0.000],
    "gt": [50.155, 2.280, 1.393, 0.660, 0.028, 18.754, 156.540, 75.497, 8.498, 0.880]
}


class ModelStochKit:
    def __init__(self, count_сore, transciption_factors, tf_probs,
                 diffuse, trans_start, protein_degrad,
                 mRNA_degrad, translation, bases):
        self.count_core = count_сore
        self.transcipton_factors = transciption_factors
        self.tf_probs = tf_probs

        self.model = etree.Element('Model')
        self.adder = adder()
        self.adder.add_child(self.model, 'Description', 'Gene expression')

        self.diffuse = diffuse
        self.trans_start = trans_start
        self.protein_degrad = protein_degrad
        self.mRNA_degrad = mRNA_degrad
        self.translation = translation
        self.bases = bases


    def add_binding_reactions(self, root):
        for core in range(self.count_core):
            repres, activ = getCounts(core)
            for site, tf_prob in enumerate(self.tf_probs):
                freeN, bindedN, onState, offState = getVariables(core, tf_prob[0], tf_prob[2])
                desc = "binding_core_{0}_tf_{1}_site_{2}".format(str(core), tf_prob[0], str(tf_prob[2]))
                rate = "{0} * {1} * {2}".format(str(tf_prob[1]), offState, freeN)

                if tf_prob[0] == "bcd" or tf_prob[0] == "cad":
                    reactants = {freeN: '1', offState: '1'}
                    products = {bindedN: '1', activ: '1', onState: '1'}
                else:
                    reactants = {freeN: '1', offState: '1'}
                    products = {bindedN: '1', repres: '1', onState: '1'}


                for isite in range (max (0, site - 2), min(site + 10, len(self.tf_probs))):
                    if isite == site:
                        continue

                    if abs(int(self.tf_probs[site][2]) - int(self.tf_probs[isite][2])) > 10:
                        continue

                    _, _, onState, offState = getVariables(core, self.tf_probs[isite][0], self.tf_probs[isite][2])
                    rate += '* {0}'.format(offState)
                    reactants[offState] = '1'

                self.adder.add_reaction(root, desc, rate, reactants, products)

        return freeN, bindedN, onState, offState

    def add_unbinding_reactions(self, root):
        species, states = [], []

        # add Protein
        for core in range(self.count_core):
            freeN, _, _, _ = getVariables(core, 'kni', 0)
            species.append(freeN)

        for core in range(self.count_core):
            repres, activ = getCounts(core)
            for site, tf_prob in enumerate(self.tf_probs):
                freeN, bindedN, onState, offState = getVariables(core, tf_prob[0], tf_prob[2])

                desc = "unbinding_core_{0}_tf_{1}_site_{2}".format(str(core), tf_prob[0], str(tf_prob[2]))
                rate = "{0} * {1} * {2}".format(str(1 - tf_prob[1]), onState, bindedN)

                if tf_prob[0] == "bcd" or tf_prob[0] == "cad":
                    reactants = {bindedN: '1', activ: '1', onState: '1'}
                    products = {freeN: '1', offState: '1'}
                else:
                    reactants = {bindedN: '1', repres: '1', onState: '1'}
                    products = {freeN: '1', offState: '1'}

                if offState not in states:
                    states.append(offState)
                if onState not in states:
                    states.append(onState)

                for isite in range (max (0, site - 10), min(site + 10, len(self.tf_probs))):
                    if isite == site:
                        continue

                    if abs(int(self.tf_probs[site][2]) - int(self.tf_probs[isite][2])) > 10:
                        continue

                    _, _, onState, offState = getVariables(core, self.tf_probs[isite][0], self.tf_probs[isite][2])
                    products[offState] = '1'

                self.adder.add_reaction(root, desc, rate, reactants, products)

                if freeN not in species:
                    species.append(freeN)
                if bindedN not in species:
                    species.append(bindedN)
        return species, states

    def add_transription_start(self, root):
        desc = "Transcription start"

        for core in range(self.count_core):
            repres, activ = getCounts(core)
            abase = self.bases[2 * core]
            rbase = self.bases[2 * core + 1]
            # rate = "(max ({0} , {1}) / max ({0} + 1, {1} + 1)) * pow ( 1.3, {0} ) * pow ( 0.7, {1} ) * {2}".format(activ, repres, self.trans_start)
            rate = "{2} * pow ( {3}, {0} ) * pow ( {4}, {1})".format(activ, repres, self.trans_start, abase, rbase)
            self.adder.add_reaction(root, desc, rate, {}, {'mRNA' + str(core): '1'})

    def add_translation_reaction(self, root):
        desc = "Translation"
        for core in range(self.count_core):
            freeN, _, _, _ = getVariables(core, 'kni', 0)
            self.adder.add_reaction(root, desc, 'mRNA' + str(core) + ' * translation',
                              {'mRNA' + str(core): '1'}, {freeN: '1'})

    def add_degradation_mRNA_reaction(self, root):
        desc = "Degradation mRNA"
        for core in range(self.count_core):
            self.adder.add_reaction(root, desc, 'mRNA' + str(core) + ' * mRNA_degrad',
                              {'mRNA' + str(core): '1'}, {})

    def add_degradation_protein_reaction(self, root):
        desc = "Degradation protein"
        for core in range(self.count_core):
            freeN, _, _, _ = getVariables(core, 'kni', 0)
            self.adder.add_reaction(root, desc, freeN + ' * protein_degrad',
                              {freeN : '1'}, {})

    def add_diff_reaction(self, root, desc, fromP, toP):
        dconst = str(self.diffuse)
        self.adder.add_reaction(root,
                          desc, dconst + ' * max(' + fromP + ' - ' + toP + ', 0.0 )',
                          {fromP: '1'}, {toP: '1'})

    def add_diffusion_reactions(self, root):
        if self.count_core == 1:
            return

        tfs = {tf for tf, prob, site in self.tf_probs}

        for tf in tfs:
            freeLeft = getFree(0, tf)
            freeLeftNeighbour = getFree(1, tf)
            self.add_diff_reaction(root, "diffusionToRight" + tf + "core" + str(0), freeLeft, freeLeftNeighbour)

            freeRight = getFree(self.count_core - 1, tf)
            freeRightNeighbour = getFree(self.count_core - 2, tf)
            self.add_diff_reaction(root, "diffusionToLeft" + tf + "core" + str(self.count_core - 1), freeRight,
                                   freeRightNeighbour)

        if self.count_core == 2:
            return

        for core in range(1, self.count_core - 1):
            for tf in tfs:
                freeN = getFree(core, tf)
                freeNLeft = getFree(core - 1, tf)
                freeNRight = getFree(core + 1, tf)

                self.add_diff_reaction(root, "diffusionToRight" + tf + "core" + str(core), freeN, freeNRight)
                self.add_diff_reaction(root, "diffusionToLeft" + tf + "core" + str(core), freeN, freeNLeft)

    def add_reactions(self):
        pl_root = etree.Element('ReactionsList')
        self.add_binding_reactions(pl_root)
        species, states = self.add_unbinding_reactions(pl_root)
        self.add_transription_start(pl_root)
        self.add_translation_reaction(pl_root)
        self.add_degradation_mRNA_reaction(pl_root)
        self.add_degradation_protein_reaction(pl_root)
        self.add_diffusion_reactions(pl_root)
        return species, states, pl_root

    def add_parameters_list(self):
        pl_root = etree.Element('ParametersList')

        parametrs = {'translation': self.translation,
                     'mRNA_degrad': self.mRNA_degrad,
                     'protein_degrad': self.protein_degrad,
                     'Transcription_start': self.trans_start}
        for key in parametrs.keys():
            self.adder.add_parametr(pl_root, key, parametrs[key])

        return pl_root

    def getIndex(self, patterns, titles):
        inds = []
        for ind, title in enumerate(titles):
            for pattern in patterns:
                if pattern in title:
                    continue
                else:
                    break
            else:
                inds.append(ind)
        return inds

    def head(self, species, states, rna_protein):
        g = open('means_head.txt', 'w')
        means_head = ["iter"] + species

        means_head += states + rna_protein
        g.write('\t'.join(means_head) + '\n')
        g.close()

        kni_indexes = self.getIndex(["Nfree", "kni"], means_head)
        self.kni_indexes = kni_indexes

        rna_indexes = self.getIndex(["mRNA"], means_head)
        self.rna_indexes = rna_indexes


    def add_species_list(self, species, states):
        sl_root = etree.Element('SpeciesList')

        for specie in species:
            _, _, core, tf = specie.split("_")
            if tf == 'kni':
                self.adder.add_specie(sl_root, specie, specie, 0.0)
            else:
                _, _, core, tf = specie.split("_")
                self.adder.add_specie(sl_root, specie, specie, int(tf_count[tf][int(core)]) if 'free' in specie else 0.0)


        for state in states:
            self.adder.add_specie(sl_root, state, state, 0.0 if 'On' in state else 1)

        rna_protein = []

        for core in range(self.count_core):
            repres, activ = getCounts(core)
            self.adder.add_specie(sl_root, repres, repres, 0)
            self.adder.add_specie(sl_root, activ, activ, 0)

            self.adder.add_specie(sl_root, 'mRNA' + str(core), 'mRNA' + str(core), 0.0)
            repres, activ = getCounts(core)

            rna_protein.append(repres)
            rna_protein.append(activ)
            rna_protein.append('mRNA' + str(core))

        self.head(species, states, rna_protein)

        return len(states) + len(species) + 2 * self.count_core + self.count_core, sl_root

    def create_model(self):
        species, states, reac_root = self.add_reactions()
        pl_root = self.add_parameters_list()
        countSpecies, s_root = self.add_species_list(species, states)

        self.adder.add_child(self.model, 'NumberOfReactions', self.adder.reaction_id - 1)
        self.adder.add_child(self.model, 'NumberOfSpecies', countSpecies)

        self.model.append(pl_root)
        self.model.append(reac_root)
        self.model.append(s_root)

        return self.model
