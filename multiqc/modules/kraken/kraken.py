#!/usr/bin/env python

""" MultiQC module to parse results from kraken  """

from __future__ import print_function

import logging
from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Kraken', anchor='kraken',
        href="http://ccb.jhu.edu/software/kraken/",
        info="a system for assigning taxonomic labels to short DNA sequences.")

        self.kraken_report_data  = dict()
        
        for f in self.find_log_files('kraken/report', filehandles=True):
            kraken_report_data[f.name] = parse_kraken_report(f)

        if len(self.kraken_report_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.kraken_report_data)))

        self.plot()


    def parse_kraken_report(self, f):
        """ Go through the report file and parse it """
        data = {}
        for line in f['f']:
            line = line.rstrip('\n').split("\t")
            percentage_reads_clade = float(line[0])
            number_reads_clade = int(line[1])
            number_reads_taxon = int(line[2])
            rank_code = line[3]
            ncbi_taxonomy_id = line[4]
            species_name = line[5].lstrip(' ')
            data[species_name] = percentage_reads_clade

        # sanity check
        # TODO
        return data

    def plot(self):
        """ Generate the plot """

        helptext = '''
        '''

        pconfig = {
            'id': 'kraken',
            'title': 'Kraken: Report Table',
        }

        self.add_section(
            anchor = 'kraken',
            description = '',
            helptext = helptext,
            plot = table.plot(self.kraken_report_data, pconfig)
        )
