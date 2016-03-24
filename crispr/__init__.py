from __future__ import absolute_import, division, print_function   # , unicode_literals
from . import config, main, analyse, cut, digest, blast, zhang, score, stretch, amplicon, prime, pcr

__all__ = [config, main] \
          + [analyse, cut, digest] \
          + [blast, zhang, score, stretch] \
          + [amplicon, prime, pcr]

__version__ = '0.0.1'

__author__ = 'Alain Domissy <alaindomissy@gmail.com>'
