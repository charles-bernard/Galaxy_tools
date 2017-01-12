#!/usr/bin/env python

import logging
import multiprocessing
from optparse import OptionParser
import os
from PyPDF2 import PdfFileReader, PdfFileMerger
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

def determine_context(analysis, diff_CG, diff_CHG, diff_CHH, diff_C, report_CG, report_CHG, report_CHH, report_C):
    if analysis == 'context_dependent':
        context = ['CG', 'CHG', 'CHH']
        diff = [diff_CG, diff_CHG, diff_CHH]
        report = [report_CG, report_CHG, report_CHH]
    elif analysis == 'context_independent':
        context = ['all_C']
        diff = [diff_C]
        report = [report_C]
    else:
        context = ['CG', 'CHG', 'CHH', 'all_C']
        diff = [diff_CG, diff_CHG, diff_CHH, diff_C]
        report = [report_CG, report_CHG, report_CHH, report_C]
    return context, diff, report

def convert_Creport_to_methkit(awk_cmd, infilename, logger):
    awk_result = subprocess.Popen(args=awk_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    awk_out, awk_err = awk_result.communicate()
    if awk_out:
        logger.info("Converting \"%s\" into methyl kit input files produced the additional awk output:\n"
                    "STDOUT: \n%s\n" % (infilename, awk_out))
    if awk_err:
        msg = "The tool failed to convert \"%s\" into methyl kit input files\n\
        ERROR: \n%s\n" % (infilename, awk_err)
        logger.critical(msg)
        sys.exit(msg)
    return

def launch_methylkit(methylkit_cmd, context, logger):
    methylkit_result = subprocess.Popen(args=methylkit_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    methylkit_out, methylkit_err = methylkit_result.communicate()
    if methylkit_out:
        logger.info("Running methylkit on %s context produced the additional output:\n"
                    "STDOUT: \n%s\n" % (context, methylkit_out))
    if methylkit_err:
        msg = "The tool failed to run methylkit differential analysis for %s context\n\
        ERROR: \n%s\n" % (context, methylkit_err)
        logger.critical(msg)
        sys.exit(msg)
    return

def merge_pdf(list_pdfs, output_filename):
    merger = PdfFileMerger()
    for pdf in list_pdfs:
        merger.append(PdfFileReader(pdf, "rb"))
    with open(output_filename, 'wb') as output:
        merger.write(output)
    return

def __main__():
    """
    WHAT IT DOES:
    1. Convert each Bismark genome-wide cytosine report into methylKit input format
    2. Run methylKit for each context
    3. Write Tables & Plots
    4. Redirect the chosen output files to Galaxy
    """

    """
    TO DO:
    TEST IF INPUT FILES ARE INDEED IN BISMARK FORMAt
    """

    """
    READ THE OPTIONS
    """
    parser = OptionParser()
    # DIRECTORIES
    parser.add_option("--tool_dir", dest='tool_dir', action='store', nargs=1, metavar='tool_dir', type="str")
    # PROFILE A
    parser.add_option("--name_A", dest="name_A", action='store', nargs=1, metavar="name_A", type="str")
    parser.add_option("--files_A", dest="files_A", action='store', nargs=1, metavar="files_A", type="str")
    # PROFILE B
    parser.add_option("--name_B", dest="name_B", action='store', nargs=1, metavar="name_B", type="str")
    parser.add_option("--files_B", dest="files_B", action='store', nargs=1, metavar="files_B", type="str")
    # SETTINGS
    parser.add_option("--n_threads", dest='n_threads', action='store', nargs=1, metavar='n_threads', type='int')
    parser.add_option("--chr_len", dest='chr_len', action='store', nargs=1, metavar='chr_len', type='str')
    parser.add_option("--win", dest="win", action='store', nargs=1, metavar="win", type="int")
    parser.add_option("--min_cov", dest="min_cov", action='store', nargs=1, metavar="min_cov", type="int")
    parser.add_option("--max_cov", dest="max_cov", action='store', nargs=1, metavar="max_cov", type="int")
    parser.add_option("--norm", action="store_true")
    parser.add_option("--qvalue", dest="qvalue", action='store', nargs=1, metavar="qvalue", type="float")
    parser.add_option("--pool", action="store_true")
    parser.add_option("--analysis_type", dest='analysis_type', action='store', nargs=1, metavar='analysis_type', type='str')
    parser.add_option("--diff_CG", dest="diff_CG", action='store', nargs=1, metavar="diff_CG", type="int")
    parser.add_option("--diff_CHG", dest="diff_CHG", action='store', nargs=1, metavar="diff_CHG", type="int")
    parser.add_option("--diff_CHH", dest="diff_CHH", action='store', nargs=1, metavar="diff_CHH", type="int")
    parser.add_option("--diff_C", dest="diff_C", action='store', nargs=1, metavar="diff_C", type="int")
    # OUTPUTS
    parser.add_option("--win_report_C", dest='win_report_C', action='store', nargs=1, metavar='win_report_C', type="str")
    parser.add_option("--win_report_CG", dest='win_report_CG', action='store', nargs=1, metavar='win_report_CG', type="str")
    parser.add_option("--win_report_CHG", dest='win_report_CHG', action='store', nargs=1, metavar='win_report_CHG', type="str")
    parser.add_option("--win_report_CHH", dest='win_report_CHH', action='store', nargs=1, metavar='win_report_CHH', type="str")
    parser.add_option("--plot", dest='plot', action='store', nargs=1, metavar='plots', type='str')
    parser.add_option("--plot_cov", action="store_true")
    parser.add_option("--plot_meth", action="store_true")
    parser.add_option("--plot_diff", action="store_true")
    parser.add_option("--plot_cor", action="store_true")
    parser.add_option("--plot_pca", action="store_true")
    parser.add_option("--cmd_report", dest='cmd_report', action='store', nargs=1, metavar='cmd_report', type='str')
    parser.add_option("--log_report", dest='log_report', action='store', nargs=1, metavar='log_report', type='str')

    """
    EXTRACT THE OPTIONS
    """
    opt, arg = parser.parse_args()

    """
    EXIT IF THERE IS NOTHING TO RETURN
    """
    if not(opt.win_report_C or opt.win_report_CG or opt.win_report_CHG or opt.win_report_CHH):
        if not opt.plot:
            sys.exit('ERROR:\nNo output files to return !\n')
        is_win_report = False
    else:
        is_win_report = True

    """
    CREATE TEMPORARY WORKING DIRECTORY
    """
    working_dir = tempfile.mkdtemp(prefix='tmp')

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


    log_report.info("1. PREPARE PARAMETERS FOR THE CONVERSION\n"
                    "FROM BISMARK GENOME-WIDE CYTOSINE METHYLATION REPORT\n"
                    "TO METHYLKIT INPUT FILES\n\n")

    """
    DETERMINE THE LISTS OF CONTEXTS
    """
    list_context, list_diff, list_win_report = determine_context(opt.analysis_type, \
        opt.diff_CG, opt.diff_CHG, opt.diff_CHH, opt.diff_C, \
        opt.win_report_CG, opt.win_report_CHG, opt.win_report_CHH, opt.win_report_C)
    nContext = len(list_context)

    """
    INITIALIZE PARALLELISM PARAMETERS
    """
    if opt.n_threads == 1:
        parallel = False
        methkit_threads = 1
    else:
        parallel = True
        n_cpu = multiprocessing.cpu_count()
        if opt.n_threads > n_cpu:
            opt.n_threads = n_cpu
        if opt.n_threads >= 3 * nContext:
            methkit_threads = (opt.n_threads - nContext) / nContext
        else:
            methkit_threads = 1

    """
    CREATE AN OUTPUT DIR PER CONTEXT
    """
    context_dir = [''] * nContext
    for i in range(0, nContext):
        context_dir[i] = os.path.join(working_dir, list_context[i])
        os.makedirs(context_dir[i])

    """
    RETRIEVE EACH INFILE NAME FROM THE STRING "opt.files_(A|B)" OF INFILE NAMES
    """
    # This is achieved by splitting the string with the pattern ";--infile--:"
    list_files_A = opt.files_A.split(";--infile--:")
    list_files_B = opt.files_B.split(";--infile--:")
    list_files = list_files_A + list_files_B
    nA = len(list_files_A)
    nB = len(list_files_B)
    n = nA + nB

    """
    DETERMINE ID OF EACH SAMPLE
    """
    list_id_files = [''] * n
    ## File Basename method:
    # for i in range (0, n):
    #    list_id_files[i] = os.path.basename(list_files[i]).rsplit('.', 1)[0]
    ## Groupname method:
    for i in range(0, nA):
        list_id_files[i]="%s-%d" % (opt.name_A, i+1)
    for i in range(0, nB):
        list_id_files[(i + nA)]="%s-%d" % (opt.name_B, i+1)

    """
    DETERMINE THE FILENAME OF EACH METHYL KIT INPUT FILE
    """
    # This is achieved by concatenating the basename of each infile name
    # with the extension "_<context>.methkit"
    list_methkit_files = [ [''] * nContext for file in range(n)]
    for i in range (0, n):
        for c in range (0, nContext):
            list_methkit_files[i][c] = os.path.join(context_dir[c], '%s_%s.methkit' % (list_id_files[i], list_context[c]))

    """
    JOIN THE METHKIT FILENAMES WHICH HAVE THE SAME BASENAME (awk output argument),
    THEN WHICH BELONG TO THE SAME CONTEXT (methkit input argument)
    """
    # "fileA1_CG.methkit;--outfile--:fileA1_CHG.methkit;--outfile--:fileA1_CHH.methkit"
    list_str_methkit_files_by_name = [''] * n
    for i in range (0, n):
        list_str_methkit_files_by_name[i] = ';--outfile--:'.join(list_methkit_files[i])

    # "fileA1_CG.methkit;--infile--;fileA2_CG.methkit"
    list_str_methkit_files_by_ctx_A = [''] * nContext
    list_str_methkit_files_by_ctx_B = [''] * nContext
    transpose_methkit_files = zip(*list_methkit_files)
    for c in range (0, nContext):
        list_str_methkit_files_by_ctx_A[c] = ';--infile--:'.join(transpose_methkit_files[c][0:nA])
        list_str_methkit_files_by_ctx_B[c] = ';--infile--:'.join(transpose_methkit_files[c][nA:n])

    """
    JOIN THE IDS WHICH BELONG TO THE GROUP A, THEN TO THE GROUP B (methkit input argument)
    """
    str_id_A = ';--id--:'.join(list_id_files[0:nA])
    str_id_B = ';--id--:'.join(list_id_files[nA:n])


    log_report.info("2. RUN THE CONVERSION\n\n")

    """
    CONVERT CYTOSINE REPORT FORMAT TO METHYLKIT INPUTFILE FORMAT WHILE
    FILTERING IRRELEVANT CYTOSINES

    CONVERSION:
    From: <chr> <base> <strand> <countC> <countT> <context> <exact_context>
    To:   <chrBase> <chr> <base> <strand> <coverage> <freqC> <freqT> (for each context)
    FILTRATION:
    retain only cytosines satisfying min & max coverage
    """

    """
    PREPARE n AWK COMMAND:
    """
    awk_script = os.path.join(opt.tool_dir, 'cytosine_report_to_methylkit_input.awk')
    awk_cmd = [''] * n

    for i in range (0, n):
        awk_cmd[i] = 'awk -v win_len=%s -v min_cov=%s -v max_cov=%s -v analysis_type=\"%s\" -v output_files=\"%s\" -f \"%s\" \"%s\"' \
            % (opt.win, opt.min_cov, opt.max_cov, opt.analysis_type, list_str_methkit_files_by_name[i], awk_script, list_files[i])
        cmd_report.info('COMMAND:\n%s\n\n' % awk_cmd[i])

    """
    LAUNCH EACH AWK COMMAND
    """
    # Conversion is heavy, so we proceed to parallel processing
    conversion_tasks = [['', '', '']] * n
    for i in range(0, n):
        conversion_tasks[i] = [awk_cmd[i], list_id_files[i], log_report]
    jobs = []
    for t in conversion_tasks:
        if parallel is True:
            p = multiprocessing.Process(target=convert_Creport_to_methkit, args=(t))
            jobs.append(p)
            p.start()
        else:
            convert_Creport_to_methkit(*t)
    # Wait until all the jobs are finished
    for j in jobs: j.join()

    log_report.info("3. PREPARE PARAMETERS FOR THE DIFFERENTIAL ANALYSES\n\n")

    """
    PREPARE PLOT OUTPUT FILES
    """
    list_pdfs = [''] * nContext
    for c in range(0, nContext):
        list_pdfs[c] = os.path.join(context_dir[c], "plot.pdf")

    """
    PREPARE nContext METHYLKIT LAUNCHER COMMAND
    """
    # Amira settings (uncomment if desired):
    parallel = False
    methkit_threads = opt.n_threads

    methylkit_script = os.path.join(opt.tool_dir, 'methylkit_launcher.R')
    methylkit_cmd = [''] * nContext
    for c in range(0, nContext):
        methylkit_cmd[c] = ("Rscript \"%s\" --n_threads %d "
                        "--name_group_A \"%s\" --files_A \"%s\" --ids_A \"%s\" "
                        "--name_group_B \"%s\" --files_B \"%s\" --ids_B \"%s\" "
                        "--context %s --diff %d --qv %.3f --win %d --norm %s --pool %s --out_dir \"%s\" " \
                        % (methylkit_script, methkit_threads, \
                        opt.name_A, list_str_methkit_files_by_ctx_A[c], str_id_A, \
                        opt.name_B, list_str_methkit_files_by_ctx_B[c], str_id_B, \
                        list_context[c], list_diff[c], opt.qvalue, opt.win, opt.norm, opt.pool, context_dir[c]))
        if is_win_report is True:
            methylkit_cmd[c]=("%s --win_report \"%s\" " % (methylkit_cmd[c], list_win_report[c]))
        if opt.plot:
            methylkit_cmd[c]=("%s --plot \"%s\" --plot_diff %s --plot_cov %s --plot_meth %s --plot_cor %s --plot_pca %s" \
                        % (methylkit_cmd[c], list_pdfs[c], opt.plot_diff, opt.plot_cov, opt.plot_meth, opt.plot_cor, opt.plot_pca))
        cmd_report.info('COMMAND:\n%s\n\n' % methylkit_cmd[c])


    log_report.info("4. RUN THE DIFFERENTIAL ANALYSES\n\n")

    """
    LAUNCH EACH METHYLKIT COMMAND
    """

    # Conversion is heavy, so we proceed to parallel processing
    diff_analysis_tasks = [['', '', '']] * nContext
    for c in range(0, nContext):
        diff_analysis_tasks[c] = [methylkit_cmd[c], list_context[c], log_report]
    jobs = []
    for t in diff_analysis_tasks:
        if parallel is True:
            p = multiprocessing.Process(target=launch_methylkit, args=(t))
            jobs.append(p)
            p.start()
        else:
            launch_methylkit(*t)
    # Wait until all the jobs are finished
    for j in jobs: j.join()

    """
    MERGE ALL THE PLOTS INTO ONE PDF
    """
    if opt.plot:
        merge_pdf(list_pdfs, opt.plot)

    """
    REMOVE TEMPORARY WORKING DIR
    """
    cleanup_before_exit(working_dir)

if __name__=="__main__": __main__()
