#!/usr/bin/env python

import logging
from optparse import OptionParser
import os
import subprocess
import sys
import tempfile


def determine_context(analysis, plot_CG, plot_CHG, plot_CHH, plot_allC,
                      hypo_dmr_CG, hypo_dmr_CHG, hypo_dmr_CHH, hypo_dmr_C, hyper_dmr_CG, hyper_dmr_CHG, hyper_dmr_CHH, hyper_dmr_C):
    if analysis == 'context_dependent':
        context = ['CG', 'CHG', 'CHH']
        infile = [[hypo_dmr_CG, hypo_dmr_CHG, hypo_dmr_CHH], [hyper_dmr_CG, hyper_dmr_CHG, hyper_dmr_CHH]]
        outfile = [plot_CG, plot_CHG, plot_CHH]
    else:
        context = ['all_C']
        infile = [[hypo_dmr_C], [hyper_dmr_C]]
        outfile = [plot_allC]
    return context, infile, outfile


def __main__():
    """
    READ THE OPTIONS
    """
    parser = OptionParser()
    parser.add_option("--tool_dir", dest='tool_dir', action='store', nargs=1, metavar='tool_dir', type="str")
    parser.add_option("--name", dest="name", action='store', nargs=1, metavar="name", type="str")
    parser.add_option("--analysis_type", dest='analysis_type', action='store', nargs=1, metavar='analysis_type', type='str')
    # REPORTS
    parser.add_option("--cmd_report", dest='cmd_report', action='store', nargs=1, metavar='cmd_report', type='str')
    ## INPUTS
    parser.add_option("--hypo_dmr_allC", dest='hypo_dmr_C', action='store', nargs=1, metavar='hypo_dmr_C', type='str')
    parser.add_option("--hyper_dmr_allC", dest='hyper_dmr_C', action='store', nargs=1, metavar='hyper_dmr_C', type='str')
    parser.add_option("--hypo_dmr_CG", dest='hypo_dmr_CG', action='store', nargs=1, metavar='hypo_dmr_CG', type='str')
    parser.add_option("--hyper_dmr_CG", dest='hyper_dmr_CG', action='store', nargs=1, metavar='hyper_dmr_CG', type='str')
    parser.add_option("--hypo_dmr_CHG", dest='hypo_dmr_CHG', action='store', nargs=1, metavar='hypo_dmr_CHG', type='str')
    parser.add_option("--hyper_dmr_CHG", dest='hyper_dmr_CHG', action='store', nargs=1, metavar='hyper_dmr_CHG', type='str')
    parser.add_option("--hypo_dmr_CHH", dest='hypo_dmr_CHH', action='store', nargs=1, metavar='hypo_dmr_CHH', type='str')
    parser.add_option("--hyper_dmr_CHH", dest='hyper_dmr_CHH', action='store', nargs=1, metavar='hyper_dmr_CHH', type='str')
    ## CHR LEN
    parser.add_option("--chr_len", dest="chr_len", action="store", nargs=1, type='str')
    ## OUTPUTS
    parser.add_option("--plot_allC", dest='plot_allC', action='store', nargs=1, metavar='plot_allC', type='str')
    parser.add_option("--plot_CG", dest='plot_CG', action='store', nargs=1, metavar='plot_CG', type='str')
    parser.add_option("--plot_CHG", dest='plot_CHG', action='store', nargs=1, metavar='plot_CHG', type='str')
    parser.add_option("--plot_CHH", dest='plot_CHH', action='store', nargs=1, metavar='plot_CHH', type='str')

    """
    EXTRACT THE OPTIONS
    """
    opt, arg = parser.parse_args()

    """
    CREATE TEMPORARY WORKING DIRECTORY
    """
    working_dir = tempfile.mkdtemp(prefix='tmp')
    print("WORKING DIR: %s" % working_dir)

    """
    INITIALIZE COMMAND REPORT
    """
    logging.basicConfig(level=logging.INFO, filename=opt.cmd_report, filemode="a+", format='%(message)s')

    """
    DETERMINE THE LISTS OF CONTEXTS
    """
    list_context, list_infiles, list_outfiles = determine_context(opt.analysis_type,
        opt.plot_CG, opt.plot_CHG, opt.plot_CHH, opt.plot_allC,
        opt.hypo_dmr_CG, opt.hypo_dmr_CHG, opt.hypo_dmr_CHH, opt.hypo_dmr_C, opt.hyper_dmr_CG, opt.hyper_dmr_CHG, opt.hyper_dmr_CHH, opt.hyper_dmr_C)
    nContext = len(list_context)

    """
    PLOT
    """
    script_path = os.path.join(opt.tool_dir, "plot_DMRs.R")
    for i in range (0, nContext):
        cmd = "Rscript %s --context %s --chr_len %s --hypo %s --hyper %s --output %s" % (script_path, list_context[i], opt.chr_len,
                                list_infiles[0][i], list_infiles[1][i], list_outfiles[i])
        logging.info("COMMAND:\n%s\n" % cmd)
        proc = subprocess.Popen(['Rscript', script_path, '--context', list_context[i], '--chr_len', opt.chr_len,
                                 '--hypo', list_infiles[0][i], '--hyper', list_infiles[1][i], '--output', list_outfiles[i]], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        proc_out, proc_err = proc.communicate()
        if proc_err:
            sys.exit(proc_err)

if __name__=="__main__": __main__()