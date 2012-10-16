#!/usr/bin/python

###########################################################################################
###
### SPAdes Repeat Resolver based on PNAS paper
###
### INPUT: De Bruijn Graph + Paired Info + Edge Sequences
### OUTPUT: FASTA, some Graph
###
### TODO: what type of contigs' prolongation should we use? do we need scaffolding?
### TODO: check outputted graph
### TODO: check misassemblies on SC_LANE1
###
###########################################################################################

import sys
import os
import fastaparser
import saveparser
from test_util import TestUtils
from graph import Graph
from rectangle_set import RectangleSet
import experimental
import logging
import check_diags
from optparse import OptionParser

def makelogger(logfilename):
    # create logger with 'rectangles'
    logger = logging.getLogger('rectangles')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfilename, mode='w')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)


def resolve(input_path, output_path, test_utils, genome):
    grp_filename = os.path.join(input_path, 'late_pair_info_counted.grp')
    sqn_filename = os.path.join(input_path, 'late_pair_info_counted.sqn')
    cvr_filename = os.path.join(input_path, 'late_pair_info_counted.cvr')
    prd_filename = os.path.join(input_path,
        'late_pair_info_counted.prd' if experimental.filter != experimental.Filter.spades else 'distance_filling_cl.prd')
    pst_filename = os.path.join(input_path,
        'distance_estimation.pst') if experimental.filter == experimental.Filter.pathsets else None
    inf_filename = os.path.join(input_path, 'late_pair_info_counted_est_params.info')
    log_filename = os.path.join(output_path, 'rectangles.log')
    config = saveparser.config(inf_filename)
    #d = utils.choose_d(config)
    d = config.median - config.RL

    makelogger(log_filename)
    logger = logging.getLogger('rectangles')

    logger.info("Rectangle Resolving %s..." % input_path)
    logger.info("d = %d..." % d)

    #################################
    # PARSE INITIAL BE BRUIJN GRAPH #
    #################################

    ingraph = Graph()
    ingraph.load(grp_filename, sqn_filename, cvr_filename)
    ingraph.check()

    maxN50 = 0
    maxgraph = None
    maxbgraph = None
    maxthreshold = 0

    rs = RectangleSet(ingraph, d, test_utils, prd_filename, config)
    if experimental.filter == experimental.Filter.pathsets:
        rs.pathsets(pst_filename)
    else:
        rs.filter(prd_filename, config)
    logger.info("  RectangleSet built.")
    if experimental.filter == experimental.Filter.spades:
        thresholds = [0.0] # everything supported by paired info
    elif experimental.filter == experimental.Filter.pathsets:
        thresholds = [-1] # everything from pathsets file
    else:
        thresholds = rs.percentiles()
    logger.info("  Checking thresholds %s..." % thresholds)
    for threshold in set(thresholds):
        logger.info("  Checking threshold %f..." % threshold)
        bgraph = rs.bgraph(threshold)
        if not bgraph.diagonals:
            continue
        bgraph.build_missing_rectangles(ingraph.K, rs)
        bgraph.condense()
        outgraph = bgraph.project(output_path)
        thisN50 = outgraph.stats(d)
        if thisN50 > maxN50:
            maxN50 = thisN50
            maxgraph = outgraph
            maxbgraph = bgraph
            maxthreshold = threshold

    #maxgraph.fasta(open(os.path.join(output_path, 'rectangles_before_graph_simplifications.fasta'), 'w'))
    #maxgraph.save(os.path.join(output_path, 'rectangles'))
    #maxbgraph.save(output_path, ingraph.K)
    maxbgraph.check_tips(ingraph.K)

    outgraph = maxbgraph.project(output_path)
    #outgraph.fasta(open(os.path.join(output_path, "after_tips.fasta"), "w"))
    #outgraph = maxbgraph.project(output_path)
    #outgraph.fasta(open(os.path.join(output_path,"after_tips_expand.fasta"),"w"))
    maxbgraph.delete_loops(ingraph.K, 500, 10)
    maxbgraph.condense()
    outgraph = maxbgraph.project(output_path)
    outgraph.fasta(open(os.path.join(output_path, "rectangles.fasta"), "w"))

    if genome:
        check_diags.check(genome, maxbgraph, maxgraph.K, open(os.path.join(output_path, "check_log.txt"), "w"),
            test_utils)

    maxgraph.stats(d)
    logger.info("Best Threshold = %d" % maxthreshold)
    logger.info("Best N50 = %d" % maxN50)


if __name__ == '__main__':
    ##########
    # PARAMS #
    ##########
    parser = OptionParser()
    parser.add_option("-g", "--genome", dest="genome", help="File with genome")
    parser.add_option("-s", "--saves", dest="saves_dir", help="Name of directory with saves")
    parser.add_option("-o", "--out", dest="out_dir", help="output folder", default="out")
    parser.add_option("-d", "--debug-logger", dest="debug_logger", help="File for debug logger", default="debug_log.txt")
    parser.add_option("-k", "--k", type=int, dest="k", help="k")
    parser.add_option("-D", "--D", type=int, dest="d", help="d")
    (options, args) = parser.parse_args()
    if not os.path.exists(options.out_dir):
        os.mkdir(options.out_dir)
    
    if options.genome and not options.saves_dir:
        if not options.k or not options.d:
            print "specify k and d"
            exit(1)
        k = options.k
        ingraph = Graph()
        _, genome = fastaparser.read_fasta(options.genome).next()
        ingraph.make_graph(genome, int(k))
        ingraph.save(os.path.join(options.out_dir, "graph"))
        rs = RectangleSet(ingraph, int(options.d), None, None, None)
        rs.filter_without_prd()
        f_left = open(os.path.join(options.out_dir, "paired_genom_contigs_1.fasta"), "w")
        f_right = open(os.path.join(options.out_dir, "paired_genom_contigs_2.fasta"), "w")
        contigs_id = 0
        for key, rect in rs.rectangles.items():
            for key, diag in rect.diagonals.items():
                e1 = rect.e1.seq
                e2 = rect.e2.seq
                f_left.write(">" + str(contigs_id) + "/1\n")
                f_left.write(e1[diag.offseta:diag.offsetc])
                f_left.write("\n")
                f_right.write(">" + str(contigs_id) + "/2\n")
                f_right.write(e2[diag.offsetb:diag.offsetd])
                f_right.write("\n")
                contigs_id += 1
        bgraph = rs.bgraph_from_genome()
        bgraph.condense()
        outgraph = bgraph.project()
        outgraph.fasta(open(os.path.join(options.out_dir, 'rectangles.fasta'), 'w'))
        exit(1)
    if len(args) != 0 or not options.saves_dir:
        parser.print_help()
        sys.exit(0)
    

    reference_information_file = os.path.join(options.saves_dir, "late_pair_info_counted_etalon_distance.txt")
    test_util = TestUtils(reference_information_file, os.path.join(options.out_dir, options.debug_logger))

    resolve(options.saves_dir, options.out_dir, test_util, options.genome)
