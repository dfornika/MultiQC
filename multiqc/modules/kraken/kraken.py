#!/usr/bin/env python

""" MultiQC module to parse results from kraken  """

from __future__ import print_function
from collections import OrderedDict
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

        self.kraken  = {}

        for f in self.find_log_files('kraken/report', filehandles=True):
            self.parse_kraken(f)
            
        if len(self.kraken) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.kraken)))

        self.add_section( plot = self.kraken_table() )


    def parse_kraken(self, f):
        """ Go through the report file and parse it """
        s_name = f['s_name']
        if getattr(config, 'kraken_fn_snames', False):
            s_name = f['s_name']
        if f['s_name'] in self.kraken:
            log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))

        self.kraken[s_name] = {}
        for line in f['f']:
            line = line.rstrip('\n').split("\t")
            percentage_reads_clade = float(line[0])
            number_reads_clade = int(line[1])
            number_reads_taxon = int(line[2])
            rank_code = line[3]
            ncbi_taxonomy_id = line[4]
            clade_name = line[5].lstrip(' ')
            if clade_name == "unclassified":
                self.kraken[s_name]['Percent Unclassified'] = percentage_reads_clade
        self.add_data_source(f, s_name)

    def kraken_table(self):
        """ Generate the table """
        headers = OrderedDict()

        table_config = {
            'namespace': 'kraken'
        }

        return table.plot(self.kraken, headers, table_config)
