from src.common import read_tf_prob, tfs, save_model
from src.stochkit_model import ModelStochKit
from lxml import etree
import sys


if __name__ == "__main__":
    tf_prob = read_tf_prob("../sitesDR.ann")
    ctfs = 10
    if len(sys.argv) > 1:
        ctfs = sys.argv[1]

    tf_prob = tf_prob[:ctfs]

    model = ModelStochKit(10, tfs, tf_prob, 0.002, 10, 0.1, 0.2, 0.2, 0.2, 68)
    doc = model.create_model()

    save_model(doc, 'model.xml')
