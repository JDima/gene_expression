__author__ = 'JDima'
from lxml import etree

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
