from src.common import read_tf_prob, tfs, save_model
from src.parser import initParserCommandLine
from src.stochkit_model import ModelStochKit

import os
from time import gmtime, strftime

def runStoch(ssa, model, t = 1000, r = 1, i = 1000):
    outdir = strftime("%Y-%m-%d_%H:%M:%S", gmtime())
    os.makedirs(outdir)
    os.system("{0} -m {1} -t {2} -r {3} -i {4} --force --seed 19 --out-dir {5} >> out". format(ssa, model, t, r, i, outdir))
    return outdir


def analyze(outdir, kni_indexes, rna_indexes):

    rna_res = [0.493, 5.751, 10.861, 50.735, 190.953, 201.132, 84.133, 38.528, 16.438, 0.066]
    kni_res = [0.826, 1.102, 1.879, 8.801, 27.885, 25.785, 10.496, 4.076, 2.344, 0.000]

    file_means = "{0}/stats/means.txt".format(outdir)
    with open(file_means, 'r', encoding='utf-8') as file:
        lines = file.readlines()
        strs = lines[-1].split("\t")

        err = 0
        for ind, kni_ind in enumerate(kni_indexes):
            err += (float(strs[kni_ind]) - kni_res[ind]) ** 2

        err /= len(kni_indexes)

        for ind, rna_ind in enumerate(rna_indexes):
            err += (float(strs[rna_ind]) - rna_res[ind]) ** 2
        err /= (len(kni_indexes) + len(rna_indexes))

    return err

if __name__ == "__main__":

    (options, args) = initParserCommandLine()

    tf_prob = read_tf_prob("../sitesDR.ann")
    # tf_prob = tf_prob[:options.sites]
    tf_prob = tf_prob[:10]

    model = ModelStochKit(count_сore=options.cores, transciption_factors=tfs, tf_probs=tf_prob,
                          diffuse = options.diffuse, trans_start = options.trans_start,
                          protein_degrad = options.kni_degrad, mRNA_degrad = options.mrna_degrad,
                          translation = options.translation, bases = options.bases)

    doc = model.create_model()

    save_model(doc, 'model' + str(options.sites) +'.xml')
    save_model(doc, 'model.xml')

    outdir = runStoch(options.tau_path, options.model, options.time, options.realizations, options.intervals)
    res = analyze(outdir, model.kni_indexes, model.rna_indexes)
    print(str(res))
