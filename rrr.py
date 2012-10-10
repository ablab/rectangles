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
import utils
import saveparser
from test_util import TestUtils
from graph import Graph
from bigraph import BGraph
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

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    grp_filename = os.path.join(input_path, 'late_pair_info_counted.grp')
    sqn_filename = os.path.join(input_path, 'late_pair_info_counted.sqn')
    cvr_filename = os.path.join(input_path, 'late_pair_info_counted.cvr')
    prd_filename = os.path.join(input_path, 'late_pair_info_counted.prd' if experimental.filter != experimental.Filter.spades else 'distance_filling_cl.prd')
    pst_filename = os.path.join(input_path, 'distance_estimation.pst') if experimental.filter == experimental.Filter.pathsets else None
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
        #print "forward", bgraph.expand()
        #bgraph.condense()
        #print "backward", bgraph.reverse_expand()
        outgraph = bgraph.project(output_path)
        thisN50 = outgraph.stats(d)
        #outgraph.check() # TODO: check error
        if thisN50 > maxN50:
            maxN50 = thisN50
            maxgraph = outgraph
            maxbgraph = bgraph
            maxthreshold = threshold

    maxgraph.fasta(open(os.path.join(output_path, 'rectangles.fasta'), 'w'))
    maxgraph.save(os.path.join(output_path, 'rectangles'))
    maxbgraph.save(output_path, ingraph.K)
    maxbgraph.check_tips(ingraph.K)
    
    outgraph = maxbgraph.project(output_path)
    outgraph.fasta(open(os.path.join(output_path,"after_tips.fasta"),"w"))
    #maxbgraph.expand()
    #maxbgraph.condense()
    maxbgraph.find_everything_about_edges([6452, 6454, 6462, 5153, 4342, 3900], ingraph.K)
    maxbgraph.find_everything_about_edges([4697, 4027, 4562, 4502], ingraph.K)
    #outgraph = maxbgraph.project(output_path)
    #outgraph.fasta(open(os.path.join(output_path,"after_tips_expand.fasta"),"w"))
    maxbgraph.delete_loops(ingraph.K, 500, 10)
    maxbgraph.condense()
    outgraph = maxbgraph.project(output_path)
    outgraph.fasta(open(os.path.join(output_path,"after_tips_expand_delete_loops.fasta"),"w"))
    
    if genome:  
      check_diags.check(genome, maxbgraph, maxgraph.K, open(os.path.join(output_path, "check_log.txt"), "w"), test_utils) 
#    maxbgraph.diag_output(output_path, maxgraph.K)

    maxgraph.stats(d)
    logger.info("Best Threshold = %d" % maxthreshold)
    logger.info("Best N50 = %d" % maxN50)


if __name__ == '__main__':

    ##########
    # PARAMS #
    ##########
    parser = OptionParser()
    parser.add_option("-g", "--genome", dest ="genome", help = "File with genome") 
    parser.add_option("-s", "--saves", dest="saves_dir", help="Name of directory with saves")
    parser.add_option("-o", "--out", dest="out_dir", help = "output folder", default = "out")
    #parser.add_option("-e", "--etalon-distance", dest="etalon_distance", help = "File with etalon distance", default = None)
    parser.add_option("-d", "--debug-logger", dest = "debug_logger", help = "File for debug logger", default = "debug_log.txt")
    parser.add_option("-k", "--k", type = int, dest = "k", help = "k")
    parser.add_option("-D", "--D", type = int, dest="d", help = "d")
    (options, args) = parser.parse_args()
    if options.genome and not options.saves_dir:  
      if not options.k or not options.d:
        print "specify k and d"
        exit(1)
      k = options.k
      ingraph = Graph()
      genome_file = open(options.genome)
      genome = ""
      line_id = genome_file.readline()
      line = genome_file.readline()
      while line:
        genome += line.strip()
        line = genome_file.readline()
      ingraph.make_graph(genome, int(k))
      ingraph.save(os.path.join(options.out_dir,"graph"))
      rs = RectangleSet(ingraph, int(options.d), None, None, None)
      rs.filter_without_prd()
      f_left = open(os.path.join(options.out_dir, "paired_genom_contigs_1.fasta"),"w")
      f_right = open(os.path.join(options.out_dir, "paired_genom_contigs_2.fasta"),"w")
      contigs_id = 0
      for key, rect in rs.rectangles.items():
        for key, diag in rect.diagonals.items():
          e1 = rect.e1.seq
          e2 = rect.e2.seq
          f_left.write(">" + str(contigs_id) + "/1\n")
          f_left.write(e1[diag.offseta:diag.offsetc])
          f_left.write("\n")
          f_right.write(">"+str(contigs_id) + "/2\n")
          f_right.write(e2[diag.offsetb:diag.offsetd])
          f_right.write("\n")
          contigs_id += 1
      bgraph = rs.bgraph_from_genome()

      bgraph.condense()
      #print "forward", bgraph.expand()
      #bgraph.condense()
      #print "backward", bgraph.reverse_expand()
      outgraph = bgraph.project()
      #outgraph.check() # TODO: check error

      outgraph.fasta(open(os.path.join(options.out_dir, 'rectangles.fasta'), 'w'))


      exit(1)

    if len(args) != 0:
      parser.print_help()
      sys.exit(0)

    input_dir = options.saves_dir
    outpath = options.out_dir
    
    reference_information_file = os.path.join(input_dir,"late_pair_info_counted_etalon_distance.txt")
    test_util = TestUtils(reference_information_file, os.path.join(outpath, options.debug_logger))
    
    resolve(input_dir, outpath, test_util, options.genome)
    if test_util.has_ref_info:
      test_util.logger.info("unaligned " + str(test_util.unaligned) +  " true_diags " + str( test_util.true_diags) +  " not_true_diags " + str( test_util.not_true_diags) +  " join correct " + str( test_util.join_correct) +  " join incorrect"  + str( test_util.join_incorrect) + " join unaligned " + str( test_util.join_unaligned))
    test_util.print_similar_diags(10)