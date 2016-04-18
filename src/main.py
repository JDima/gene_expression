from src.stochkit_model import ModelStochKit

__author__ = 'JDima'

from lxml import etree

tfs = ["hb", "Kr", "gt", "kni", "bcd", "cad", "tll", "hkb"]


# def read_tf_prob(file):
#     tf_prob = []
#     for line in open(file):
#         _, _, tf, _, _, prob, = line.split()
#         tf_prob.append((tfs[int(tf)], float(prob)))
#     return tf_prob

def read_tf_prob(file):
    tf_prob = {}
    for line in open(file):
        site, _, tf, _, _, prob, = line.split()
        tf_prob[site]  = ( tfs[int(tf)], float(prob))
    return tf_prob


def save_model(doc, out='/../model.xml'):
    f = open(out, 'wb')
    f.write(etree.tostring(doc, pretty_print=True))
    f.close()


if __name__ == "__main__":
    tf_prob = read_tf_prob("../sitesDR.ann")

    tf_prob = tf_prob[:2]

    model = ModelStochKit(2, tfs, tf_prob, 0.002, 10)
    doc = model.create_model()

    save_model(doc, 'model.xml')
