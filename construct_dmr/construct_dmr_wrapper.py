#!/usr/bin/env python

import logging
import multiprocessing
from optparse import OptionParser
import os
import re
import shutil
import subprocess
import sys
import tempfile


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    return

def create_logger(name, filename):
    # Initialize the logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    # Create a file handler
    handler = logging.FileHandler(filename, mode="a+")
    handler.setLevel(logging.INFO)
    # Create a logging format
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    # Add the handlers to the logger
    logger.addHandler(handler)
    return logger

def determine_context(analysis, min_C_CG, min_C_CHG, min_C_CHH, min_C_allC, \
    infile_CG, infile_CHG, infile_CHH, infile_allC, \
    hypo_dmr_CG, hypo_dmr_CHG, hypo_dmr_CHH, hypo_dmr_C, hyper_dmr_CG, hyper_dmr_CHG, hyper_dmr_CHH, hyper_dmr_C, \
    hypo_dmr_CG_fig, hypo_dmr_CHG_fig, hypo_dmr_CHH_fig, hypo_dmr_C_fig, hyper_dmr_CG_fig, hyper_dmr_CHG_fig, hyper_dmr_CHH_fig, hyper_dmr_C_fig, \
    hypo_dmr_CG_gff, hypo_dmr_CHG_gff, hypo_dmr_CHH_gff, hypo_dmr_C_gff, hyper_dmr_CG_gff, hyper_dmr_CHG_gff, hyper_dmr_CHH_gff, hyper_dmr_C_gff):
    if analysis == 'context_dependent':
        context = ['CG', 'CHG', 'CHH']
        min_C = [min_C_CG, min_C_CHG, min_C_CHH]
        infile = [infile_CG, infile_CHG, infile_CHH, infile_allC]
        out_file = [[hypo_dmr_CG, hypo_dmr_CHG, hypo_dmr_CHH], [hyper_dmr_CG, hyper_dmr_CHG, hyper_dmr_CHH]]
        out_file_fig = [[hypo_dmr_CG_fig, hypo_dmr_CHG_fig, hypo_dmr_CHH_fig], [hyper_dmr_CG_fig, hyper_dmr_CHG_fig, hyper_dmr_CHH_fig]]
        out_file_gff = [[hypo_dmr_CG_gff, hypo_dmr_CHG_gff, hypo_dmr_CHH_gff], [hyper_dmr_CG_gff, hyper_dmr_CHG_gff, hyper_dmr_CHH_gff]]
    else:
        context = ['all_C']
        min_C = [min_C_allC]
        infile = [infile_allC]
        out_file = [[hypo_dmr_C], [hyper_dmr_C]]
        out_file_fig = [[hypo_dmr_C_fig], [hyper_dmr_C_fig]]
        out_file_gff = [[hypo_dmr_C_gff], [hyper_dmr_C_gff]]
    return context, min_C, infile, out_file, out_file_fig, out_file_gff

def produce_stats(R_cmd, context, logger):
    produce_stats_result = subprocess.Popen(args=R_cmd, shell=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
    produce_stats_out, produce_stats_err = produce_stats_result.communicate()
    if produce_stats_out:
        logger.info("Running generate_stats_for_filtering on %s context produced the additional output:\n"
                    "STDOUT: \n%s\n" % (context, produce_stats_out))
    if produce_stats_err:
        msg = "The tool failed to generate the stats recquired for filtering the windows for %s context\n\
        ERROR: \n%s\n" % (context, produce_stats_err)
        logger.critical(msg)
        sys.exit(msg)
    return

def construct_dmr(awk_cmd, context, logger):
    construct_dmr_result = subprocess.Popen(args=awk_cmd, shell=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    construct_dmr_out, construct_dmr_err = construct_dmr_result.communicate()
    if construct_dmr_out:
        logger.info("Running construct_dmr on %s context produced the additional output:\n"
                    "STDOUT: \n%s\n" % (context, construct_dmr_out))
    if construct_dmr_err:
        msg = "The tool failed to construct the dmrs for %s context\n\
        ERROR: \n%s\n" % (context, construct_dmr_err)
        logger.critical(msg)
        sys.exit(msg)
    return

def __main__():

    parallel = True

    """
    READ THE OPTIONS
    """
    parser = OptionParser()
    parser.add_option("--tool_dir", dest='tool_dir', action='store', nargs=1, metavar='tool_dir', type="str")
    parser.add_option("--name", dest="name", action='store', nargs=1, metavar="name", type="str")
    parser.add_option("--analysis_type", dest='analysis_type', action='store', nargs=1, metavar='analysis_type', type='str')
    parser.add_option("--infile_CG", dest="infile_CG", action='store', nargs=1, metavar="infile_CG", type="str")
    parser.add_option("--infile_CHG", dest="infile_CHG", action='store', nargs=1, metavar="infile_CHG", type="str")
    parser.add_option("--infile_CHH", dest="infile_CHH", action='store', nargs=1, metavar="infile_CHH", type="str")
    parser.add_option("--infile_allC", dest="infile_allC", action='store', nargs=1, metavar="infile_allC", type="str")
    # SETTINGS
    parser.add_option("--min_C_CG", dest="min_C_CG", action='store', nargs=1, metavar="min_C_CG", type="int")
    parser.add_option("--min_C_CHG", dest="min_C_CHG", action='store', nargs=1, metavar="min_C_CHG", type="int")
    parser.add_option("--min_C_CHH", dest="min_C_CHH", action='store', nargs=1, metavar="min_C_CHH", type="int")
    parser.add_option("--min_C_allC", dest="min_C_allC", action='store', nargs=1, metavar="min_C_allC", type="int")
    parser.add_option("--min_cov", dest="min_cov", action='store', nargs=1, metavar="min_cov", type="int")
    parser.add_option("--common_C_ratio", dest="min_ratio", action='store', nargs=1, metavar="min_ratio", type="float")
    parser.add_option("--n_gap", dest="n_gap", action='store', nargs=1, metavar="nb_gap", type="int")
    # REPORTS
    parser.add_option("--cmd_report", dest='cmd_report', action='store', nargs=1, metavar='cmd_report', type='str')
    parser.add_option("--log_report", dest='log_report', action='store', nargs=1, metavar='log_report', type='str')
    ## OUPUTS
    parser.add_option("--hypo_dmr_C", dest='hypo_dmr_C', action='store', nargs=1, metavar='hypo_dmr_C', type='str')
    parser.add_option("--hypo_dmr_C_fig", dest='hypo_dmr_C_fig', action='store', nargs=1, metavar='hypo_dmr_C_fig', type='str')
    parser.add_option("--hypo_dmr_C_gff", dest='hypo_dmr_C_gff', action='store', nargs=1, metavar='hypo_dmr_C_gff', type='str')
    parser.add_option("--hyper_dmr_C", dest='hyper_dmr_C', action='store', nargs=1, metavar='hyper_dmr_C', type='str')
    parser.add_option("--hyper_dmr_C_fig", dest='hyper_dmr_C_fig', action='store', nargs=1, metavar='hyper_dmr_C_fig', type='str')
    parser.add_option("--hyper_dmr_C_gff", dest='hyper_dmr_C_gff', action='store', nargs=1, metavar='hyper_dmr_C_gff', type='str')
    parser.add_option("--hypo_dmr_CG", dest='hypo_dmr_CG', action='store', nargs=1, metavar='hypo_dmr_CG', type='str')
    parser.add_option("--hypo_dmr_CG_fig", dest='hypo_dmr_CG_fig', action='store', nargs=1, metavar='hypo_dmr_CG_fig', type='str')
    parser.add_option("--hypo_dmr_CG_gff", dest='hypo_dmr_CG_gff', action='store', nargs=1, metavar='hypo_dmr_CG_gff', type='str')
    parser.add_option("--hyper_dmr_CG", dest='hyper_dmr_CG', action='store', nargs=1, metavar='hyper_dmr_CG', type='str')
    parser.add_option("--hyper_dmr_CG_fig", dest='hyper_dmr_CG_fig', action='store', nargs=1, metavar='hyper_dmr_CG_fig', type='str')
    parser.add_option("--hyper_dmr_CG_gff", dest='hyper_dmr_CG_gff', action='store', nargs=1, metavar='hyper_dmr_CG_gff', type='str')
    parser.add_option("--hypo_dmr_CHG", dest='hypo_dmr_CHG', action='store', nargs=1, metavar='hypo_dmr_CHG', type='str')
    parser.add_option("--hypo_dmr_CHG_fig", dest='hypo_dmr_CHG_fig', action='store', nargs=1, metavar='hypo_dmr_CHG_fig', type='str')
    parser.add_option("--hypo_dmr_CHG_gff", dest='hypo_dmr_CHG_gff', action='store', nargs=1, metavar='hypo_dmr_CHG_gff', type='str')
    parser.add_option("--hyper_dmr_CHG", dest='hyper_dmr_CHG', action='store', nargs=1, metavar='hyper_dmr_CHG', type='str')
    parser.add_option("--hyper_dmr_CHG_fig", dest='hyper_dmr_CHG_fig', action='store', nargs=1, metavar='hyper_dmr_CHG_fig', type='str')
    parser.add_option("--hyper_dmr_CHG_gff", dest='hyper_dmr_CHG_gff', action='store', nargs=1, metavar='hyper_dmr_CHG_gff', type='str')
    parser.add_option("--hypo_dmr_CHH", dest='hypo_dmr_CHH', action='store', nargs=1, metavar='hypo_dmr_CHH', type='str')
    parser.add_option("--hypo_dmr_CHH_fig", dest='hypo_dmr_CHH_fig', action='store', nargs=1, metavar='hypo_dmr_CHH_fig', type='str')
    parser.add_option("--hypo_dmr_CHH_gff", dest='hypo_dmr_CHH_gff', action='store', nargs=1, metavar='hypo_dmr_CHH_gff', type='str')
    parser.add_option("--hyper_dmr_CHH", dest='hyper_dmr_CHH', action='store', nargs=1, metavar='hyper_dmr_CHH', type='str')
    parser.add_option("--hyper_dmr_CHH_fig", dest='hyper_dmr_CHH_fig', action='store', nargs=1, metavar='hyper_dmr_CHH_fig', type='str')
    parser.add_option("--hyper_dmr_CHH_gff", dest='hyper_dmr_CHH_gff', action='store', nargs=1, metavar='hyper_dmr_CHH_gff', type='str')

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
    INITIALIZE COMMAND & LOG REPORTS
    """
    if opt.cmd_report:
        cmd_report_name = opt.cmd_report
    else:
        cmd_report_name = os.path.join(working_dir, 'tmp_cmd_report.txt')
    if opt.log_report:
        log_report_name = opt.log_report
    else:
        log_report_name = os.path.join(working_dir, 'tmp_log_report.txt')
    cmd_report = create_logger('cmd_report', cmd_report_name)
    log_report = create_logger('log_report', log_report_name)

    """
    DETERMINE THE LISTS OF CONTEXTS
    """
    list_context, list_min_C, list_infiles, list_outfiles, list_outfiles_fig, list_outfiles_gff = determine_context(opt.analysis_type,
        opt.min_C_CG, opt.min_C_CHG, opt.min_C_CHH, opt.min_C_allC,
        opt.infile_CG, opt.infile_CHG, opt.infile_CHH, opt.infile_allC,
        opt.hypo_dmr_CG, opt.hypo_dmr_CHG, opt.hypo_dmr_CHH, opt.hypo_dmr_C, opt.hyper_dmr_CG, opt.hyper_dmr_CHG, opt.hyper_dmr_CHH, opt.hyper_dmr_C,
        opt.hypo_dmr_CG_fig, opt.hypo_dmr_CHG_fig, opt.hypo_dmr_CHH_fig, opt.hypo_dmr_C_fig, opt.hyper_dmr_CG_fig, opt.hyper_dmr_CHG_fig, opt.hyper_dmr_CHH_fig, opt.hyper_dmr_C_fig,
        opt.hypo_dmr_CG_gff, opt.hypo_dmr_CHG_gff, opt.hypo_dmr_CHH_gff, opt.hypo_dmr_C_gff, opt.hyper_dmr_CG_gff, opt.hyper_dmr_CHG_gff, opt.hyper_dmr_CHH_gff, opt.hyper_dmr_C_gff)
    nContext = len(list_context)

    """
    CREATE AN OUTPUT DIR PER CONTEXT
    """
    context_dir = [''] * nContext
    for i in range(0, nContext):
        context_dir[i] = os.path.join(working_dir, list_context[i])
        os.makedirs(context_dir[i])

    """
    DETERMINE THE STRUCTURE OF THE METHYLKIT OUTPUT REPORT
    """
    nA = [''] * nContext
    nB = [''] * nContext
    is_pool = [''] * nContext
    for i in range(0, nContext):
        with open(list_infiles[i], 'r') as f: first_line = f.readline().strip()
        nA[i] = len(re.findall(r'nb_Cs_A', first_line))
        nB[i] = len(re.findall(r'nb_Cs_A', first_line))
        nCov = len(re.findall(r'tot_cov_(A|B)', first_line))
        if(nCov > 2):
            is_pool[i] = False
        else:
            is_pool[i] = True

    """
    DETERMINE THE LIST OF R OUTPUT FILES
    """
    R_output_files = [''] * nContext
    for i in range(0, nContext):
        R_output_files[i] = os.path.join(context_dir[i], "filtering_stats_table_%s.txt" % list_context[i])

    """
    PREPARE nContext R COMMAND:
    """
    R_script = os.path.join(opt.tool_dir, 'generate_data_for_filtering.R')
    R_cmd = [''] * nContext

    for i in range (0, nContext):
        R_cmd[i] = ("Rscript \"%s\" "
                    "--infile \"%s\" --outfile \"%s\" "
                    "--nA %d --nB %d --pool %s "
                    % (R_script,
                    list_infiles[i], R_output_files[i],
                    int(nA[i]), int(nB[i]), is_pool[i]))
        cmd_report.info('COMMAND:\n%s\n\n' % R_cmd[i])

    log_report.info("1. GENERATE STATISTICS IN ORDER TO FILTER IRRELEVANT GENOMIC WINDOWS\n\n")

    """
    LAUNCH EACH R COMMAND
    """
    produce_stats_tasks = [['', '', '']] * nContext
    for i in range(0, nContext):
        produce_stats_tasks[i] = [R_cmd[i], list_context[i], log_report]
        jobs = []
        for t in produce_stats_tasks:
            if parallel is True:
                p = multiprocessing.Process(target=produce_stats, args=(t))
                jobs.append(p)
                p.start()
            else:
                produce_stats(*t)
        # Wait until all the jobs are finished
        for j in jobs: j.join()

    """
    PREPARE nContext awk COMMAND:
    """
    awk_script = os.path.join(opt.tool_dir, 'join_dmw_into_dmr.awk')
    awk_cmd = [''] * nContext

    for i in range(0, nContext):
        awk_cmd[i] = ("awk -v outdir=\"%s\" -v name=\"%s\" -v context=%s "
                    "-v min_C=%d -v min_cov=%d -v min_ratio=%.3f -v n_gap=%d "
                    "-f \"%s\" \"%s\""
                    % (context_dir[i], opt.name, list_context[i],
                    list_min_C[i], opt.min_cov, opt.min_ratio, opt.n_gap,
                    awk_script, R_output_files[i]))
        cmd_report.info('COMMAND:\n%s\n\n' % awk_cmd[i])

    log_report.info("2. CONSTRUCT THE DMRS\n\n")

    """
    LAUNCH EACH AWK COMMAND
    """
    construct_dmr_tasks = [['', '', '']] * nContext
    for i in range(0, nContext):
        construct_dmr_tasks[i] = [awk_cmd[i], list_context[i], log_report]
    jobs = []
    for t in construct_dmr_tasks:
        if parallel is True:
            p = multiprocessing.Process(target=construct_dmr, args=(t))
            jobs.append(p)
            p.start()
        else:
            construct_dmr(*t)
    # Wait until all the jobs are finished
    for j in jobs: j.join()

    """
    REDIRECT OUTPUT
    """
    dmr_type = ['hypo', 'hyper']
    for d in range(0, 2):
        if opt.analysis_type == 'context_dependent':
            # the filenames are determined in the awk script
            for i in range(0, nContext):
                shutil.move(os.path.join(context_dir[i], dmr_type[d] + "_differenceDMRs.txt"), (list_outfiles[d][i]))
                shutil.move(os.path.join(context_dir[i], dmr_type[d] + "_differenceDMRs_for_figure_and_annotation.txt"), list_outfiles_fig[d][i])
                shutil.move(os.path.join(context_dir[i], dmr_type[d] + "_differenceDMRs.gff"), list_outfiles_gff[d][i])

    """
    REMOVE TEMPORARY WORKING DIR
    """
    #cleanup_before_exit(working_dir)

if __name__=="__main__": __main__()
