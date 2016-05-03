from src.stochkit_model import ModelStochKit
from lxml import etree
import sys


tfs = ["hb", "Kr", "gt", "kni", "bcd", "cad", "tll", "hkb"]


def read_tf_prob(file):
    tf_prob = []
    for line in open(file):
        site, _, tf, _, _, prob, = line.split()
        tf_prob.append((tfs[int(tf)], float(prob), site))
    return sorted(tf_prob, key=lambda tup: int(tup[2]))


def save_model(doc, out='/../model.xml'):
    f = open(out, 'wb')
    f.write(etree.tostring(doc, pretty_print=True))
    f.close()


if __name__ == "__main__":
    tf_prob = read_tf_prob("../sitesDR.ann")
    ctfs = 10
    if len(sys.argv) > 1:
        ctfs = sys.argv[1]

    tf_prob = tf_prob[:ctfs]

    model = ModelStochKit(10, tfs, tf_prob, 0.002, 10, 0.1, 0.2, 0.2, 0.2, 68)
    doc = model.create_model()

    save_model(doc, 'model.xml')
