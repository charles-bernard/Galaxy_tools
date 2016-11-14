#!/usr/bin/env python

import logging
import os
import shutil
import subprocess
import sys
import tempfile
from glob import glob
from optparse import OptionParser

def __main__():
    parser = OptionParser()
    parser.add_option("--tool_dir", dest = 'tool_dir', action = 'store', nargs = 1, metavar = 'tool_dir', type = "str")
    parser.add_option("--cytosine_report", dest='cytosine_report', action='store', nargs=1, metavar='cytosine_report', type="str")
    parser.add_option("--context", action="store_true")
    parser.add_option("--CpG_bedgraph", dest='CpG_bedgraph', action="store", nargs = 1, metavar = 'CpG_bedgraph', type='str')
    parser.add_option("--CHH_bedgraph", dest='CHH_bedgraph', action="store", nargs=1, metavar='CHH_bedgraph', type='str')
    parser.add_option("--CHG_bedgraph", dest='CHG_bedgraph', action="store", nargs=1, metavar='CHG_bedgraph', type='str')
    parser.add_option("--CXX_bedgraph", dest='CXX_bedgraph', action="store", nargs=1, metavar='CXX_bedgraph', type='str')
    parser.add_option("--coverage_bedgraph", dest='coverage_bedgraph', action="store", nargs=1, metavar='coverage_bedgraph', type='str')
    parser.add_option("--tdf", action="store_true")
    parser.add_option("--igv_genome", dest='igv_genome', action='store', nargs=1, metavar='igv_genome', type='str')
    parser.add_option("--CpG_tdf", dest='CpG_tdf', action="store", nargs=1, metavar='CpG_tdf', type='str')
    parser.add_option("--CHH_tdf", dest='CHH_tdf', action="store", nargs=1, metavar='CHH_tdf', type='str')
    parser.add_option("--CHG_tdf", dest='CHG_tdf', action="store", nargs=1, metavar='CHG_tdf', type='str')
    parser.add_option("--CXX_tdf", dest='CXX_tdf', action="store", nargs=1, metavar='CXX_tdf', type='str')
    parser.add_option("--coverage_tdf", dest='coverage_tdf', action="store", nargs=1, metavar='coverage_tdf', type='str')
    parser.add_option("--log", dest='log_report', action="store", nargs=1, metavar='log_report')

    options, args = parser.parse_args()

    if options.log_report:
    	logging.basicConfig(level=logging.INFO, filename=options.log_report, filemode="a+", format='%(message)s')
    else:
    	logging.basicConfig(level=logging.INFO, filename='galaxy_log_report.log', filemode="a+", format='%(message)s')
    logging.info('_______________________________________________________________\n')
    logging.info('List of options: %s\n' % options)

    script_path=os.path.join(options.tool_dir, 'bismark2bedgraph.sh')
    tmp_dir = tempfile.mkdtemp(prefix='tmp')

    if options.context:
        context_option = '-c'
    else:
        context_option = ''

    if options.tdf:
        tdf_option = '--tdf'
        genome_option = '--igv_genome'
        genome = options.igv_genome
    else:
        tdf_option = ''
        genome_option = ''
        genome = ''

    logging.info('_______________________________________________________________\n')
    logging.info('COMMAND:\n')
    logging.info("bash %s -i %s %s %s %s %s -e current_job -o %s\n" % (script_path, options.cytosine_report, context_option, tdf_option, genome_option, genome, tmp_dir))

    proc = subprocess.Popen(['bash', script_path, '--tool_dir', options.tool_dir, '-i', options.cytosine_report,\
                             context_option, tdf_option, genome_option, genome, '-e', 'current_job', '-o', tmp_dir], \
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    proc_out, proc_err = proc.communicate()

    if proc_out:
    	logging.info('_______________________________________________________________\n')
    	logging.info('Tool STDOUT:\n')
    	logging.info(proc_out)
    if proc_err:
    	logging.info('_______________________________________________________________\n')
    	logging.info('ERROR:\n')
    	msg = 'The tool crashed with the following error:\n%s\n' % proc_err
    	logging.critical(msg)
    	sys.exit(msg)

    if options.context:
        shutil.move( glob(os.path.join(tmp_dir, '*_CpG.bedgraph'))[0], options.CpG_bedgraph )
        shutil.move( glob(os.path.join(tmp_dir, '*_CHG.bedgraph'))[0], options.CHG_bedgraph )
        shutil.move( glob(os.path.join(tmp_dir, '*_CHH.bedgraph'))[0], options.CHH_bedgraph )
        if options.tdf:
            shutil.move(glob(os.path.join(tmp_dir, '*_CpG.tdf'))[0], options.CpG_tdf)
            shutil.move(glob(os.path.join(tmp_dir, '*_CHG.tdf'))[0], options.CHG_tdf)
            shutil.move(glob(os.path.join(tmp_dir, '*_CHH.tdf'))[0], options.CHH_tdf)
    else:
        shutil.move(glob(os.path.join(tmp_dir, '*_CXX.bedgraph'))[0], options.CXX_bedgraph)
        if options.tdf:
            shutil.move(glob(os.path.join(tmp_dir, '*_CXX.tdf'))[0], options.CXX_tdf)
    shutil.move(glob(os.path.join(tmp_dir, '*_coverage.bedgraph'))[0], options.coverage_bedgraph)
    if options.tdf:
        shutil.move(glob(os.path.join(tmp_dir, '*_coverage.tdf'))[0], options.coverage_tdf)

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

if __name__=="__main__": __main__()