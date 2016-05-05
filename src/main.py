from src.common import read_tf_prob, tfs, save_model
from src.stochkit_model import ModelStochKit
from lxml import etree
import sys


if __name__ == "__main__":
    tf_prob = read_tf_prob("../sitesDR.ann")
    # ctfs = len(tf_prob)
    ctfs = 10
    if len(sys.argv) > 1:
        ctfs = sys.argv[1]

    tf_prob = tf_prob[:ctfs]

    # model = ModelStochKit(10, tfs, tf_prob, 0.002, 10, 0.1, 0.2, 0.2, 0.2, 68)
    model = ModelStochKit(count_сore=10, transciption_factors=tfs, tf_probs=tf_prob,
                          diffuse = 0.002, h_dist = 10, trans_start = 0.002069, Protein_degrad = 0.000009627,
                            mRNA_degrad = 0.000019254, translation = 0.005172, start_cout_tf = 68)
    doc = model.create_model()

    save_model(doc, 'model' + str(ctfs) +'.xml')
