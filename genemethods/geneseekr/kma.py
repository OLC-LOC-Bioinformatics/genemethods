#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, MetadataObject
from genemethods.MLSTsippr.mlst import GeneSippr as MLSTSippr
from genemethods.sipprverse_reporter.reports import Reports
from genemethods.geneseekr.geneseekr import GeneSeekr
from genemethods.geneseekr.parser import Parser
import multiprocessing
from glob import glob
import logging
import os

__author__ = 'adamkoziol'


class KMA(object):

    def seekr(self):
        pass

    def __init__(self, args, analysistype='geneseekr', cutoff=70, program='kma', genus_specific=False, unique=False,
                 pipeline=True):
        pass
