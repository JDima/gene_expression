from src.common import read_tf_prob, tfs, save_model
from src.parser import initParserCommandLine
from src.stochkit_model import ModelStochKit

import os


def runStoch(ssa, model, t = 1000, r = 1, i = 1000):
    os.system("{0} -m {1} -t {2} -r {3} -i {4} --force --seed 19". format(ssa, model, t, r, i))

if __name__ == "__main__":

    (options, args) = initParserCommandLine()

    tf_prob = read_tf_prob("../sitesDR.ann")
    # tf_prob = tf_prob[:options.sites]
    tf_prob = tf_prob[:10]

    model = ModelStochKit(count_сore=options.cores, transciption_factors=tfs, tf_probs=tf_prob,
                          diffuse = options.diffuse, trans_start = options.trans_start,
                          protein_degrad = options.kni_degrad, mRNA_degrad = options.mrna_degrad,
                          translation = options.translation, start_cout_tf = 68)

    doc = model.create_model()

    save_model(doc, 'model' + str(options.sites) +'.xml')
    save_model(doc, 'model.xml')

    runStoch(options.tau_path, options.model, options.time, options.realizations, options.intervals)