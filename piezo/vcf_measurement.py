#! /usr/bin/env python

from datreant import Treant

class VCFMeasurement(Treant):

    _treanttype='VCFMeasurement'

    def __init__(self, vcf_path, new=False, categories=None, tags=["vcf",'genetics'], species='M.tuberculosis'):

        Treant.__init__(self, vcf_path, categories=categories, tags=tags)
