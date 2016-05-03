import os
import sys
from src.common import read_tf_prob, tfs, save_model
from src.stochkit_model import ModelStochKit

__author__ = 'JDima'

import pytest

class TestBoundValues:

    # No reactions
    def test_tf_no_tf(self):
        tf_prob = read_tf_prob("../sitesDR.ann")
        ctfs = 4

        tf_prob = tf_prob[:ctfs]

        model = ModelStochKit(3, tfs, tf_prob, 0.002, 10, 0.1, 0.2, 0.2, 0.2, 0)
        doc = model.create_model()

        save_model(doc, 'result//test_tf_no_tf.xml')

    # All reasctions, but kni is 0, because does not have tfs in first 4 in file(sitesDR.ann)
    def test_tf_only_bind_unbind(self):
        tf_prob = read_tf_prob("../sitesDR.ann")
        ctfs = 4

        tf_prob = tf_prob[:ctfs]

        model = ModelStochKit(3, tfs, tf_prob, 0.002, 10, 0.0, 0.0, 0.0, 0.0, 68)
        doc = model.create_model()

        save_model(doc, 'result//test_tf_only_bind_unbind.xml')


    def test_good_high_trans_no_degra(self):
        tf_prob = read_tf_prob("../sitesDR.ann")
        ctfs = 4

        tf_prob = tf_prob[:ctfs]

        model = ModelStochKit(3, tfs, tf_prob, 0.002, 10, 1.0, 0.0, 0.0, 1.0, 68)
        doc = model.create_model()

        save_model(doc, 'result//test_good_high_trans_no_degra.xml')

    def test_good_low_trans_high_degra(self):
        tf_prob = read_tf_prob("../sitesDR.ann")
        ctfs = 4

        tf_prob = tf_prob[:ctfs]

        model = ModelStochKit(3, tfs, tf_prob, 0.002, 10, 0.1, 0.9, 0.9, 0.1, 68)
        doc = model.create_model()

        save_model(doc, 'result//test_good_low_trans_high_degra.xml')

    def test_good_no_trans(self):
        tf_prob = read_tf_prob("../sitesDR.ann")
        ctfs = 4

        tf_prob = tf_prob[:ctfs]

        model = ModelStochKit(3, tfs, tf_prob, 0.002, 10, 1.0, 0.0, 0.0, 0.0, 68)
        doc = model.create_model()

        save_model(doc, 'result//test_good_no_trans.xml')