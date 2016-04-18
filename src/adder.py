__author__ = 'JDima'
from lxml import etree


class adder:
    def __init__(self):
        self.reaction_id = 1

    def add_child(self, root, child_name, text=None):
        child = etree.Element(child_name)
        if text != None:
            child.text = str(text)
        root.append(child)
        return child

    def add_reaction(self, root, description, rate, reactants, products, type='customized'):
        reaction = self.add_child(root, 'Reaction')
        self.add_child(reaction, 'Id', self.reaction_id)
        self.add_child(reaction, 'Description', description)
        self.add_child(reaction, 'Type', type)
        self.add_child(reaction, 'PropensityFunction', rate)

        if reactants:
            ereactants = self.add_child(reaction, 'Reactants')
            ereactants.extend(
                [etree.Element('SpeciesReference', {'id': reactant, 'stoichiometry': reactants[reactant]}) for reactant
                 in reactants.keys()])
            reaction.append(ereactants)

        if products:
            eproducts = self.add_child(reaction, 'Products')
            eproducts.extend(
                [etree.Element('SpeciesReference', {'id': product, 'stoichiometry': products[product]}) for product in
                 products.keys()])
            reaction.append(eproducts)

        self.reaction_id += 1

        return reaction

    def add_parametr(self, pl_root, key, value):
        parameter = self.add_child(pl_root, 'Parameter')
        id = etree.Element('id')
        id.text = key
        expression = etree.Element('Expression')
        expression.text = str(value)
        parameter.extend([id, expression])

    def add_specie(self, root, id, desc, init):
        sroot = self.add_child(root, 'Species')
        self.add_child(sroot, 'Id', id)
        self.add_child(sroot, 'Description', desc)
        self.add_child(sroot, 'InitialPopulation', init)
