import logging
import sys 
import math
import os
import utils
import saveparser
from diagonal import Diagonal
from utils import conjugate
import experimental

def makelogger(logfilename):
  logger = logging.getLogger('debug')
  logger.setLevel(logging.DEBUG)
  fh = logging.FileHandler(logfilename, mode='w')
  fh.setLevel(logging.DEBUG)
  ch = logging.StreamHandler()
  ch.setLevel(logging.ERROR)
  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
  fh.setFormatter(formatter)
  ch.setFormatter(formatter)
  logger.addHandler(fh)
  logger.addHandler(ch)



def parse_ref_info(fin):
  ref_info = dict()
  for line in fin:
    block = line.split()
    edge_id = int(block[0].strip())
    ref_begin = int(block[1].strip())
    ref_end = int(block[2].strip())
    contig_begin = int(block[3].strip())
    contig_end = int(block[4].strip())
    if edge_id < 0:
      edge_id = -edge_id
      if edge_id not in ref_info:
        ref_info[edge_id] = []
      ref_info[edge_id].append(-(ref_begin - contig_begin))
    else:
      if contig_end > contig_begin:
        if edge_id not in ref_info:
          ref_info[edge_id] = []
        ref_info[edge_id].append(ref_begin - contig_begin)
  #print "ref_info", len(ref_info)
  return ref_info

class TestUtils(object):
  def __init__(self, edge_aligned_file, log_file_name):
    makelogger(log_file_name)
    self.logger = logging.getLogger('debug')
    self.has_ref_info = True  
    if edge_aligned_file:
      try:
        with open(edge_aligned_file) as f: 
           self.ref_info = parse_ref_info(open(edge_aligned_file))
      except IOError as e:
           self.ref_info = dict()
           self.has_ref_info = False
    else:
      self.ref_info = dict()
      self.has_ref_info = False
    self.unaligned = 0
    self.not_true_diags = 0
    self.true_diags = 0
    self.join_correct = 0
    self.join_incorrect = 0
    self.join_unaligned = 0
    self.similar_diags = dict()
  
  def __get_ref_info(self, eid):
    if eid in self.ref_info:
      return self.ref_info[eid]
    return []

  def add_to_diags(self, diag):
    rect = diag.rectangle
    e1 = rect.e1
    e2 = rect.e2
    D = diag.D
    if (e1, e2) not in self.similar_diags:
      self.similar_diags[(e1, e2)] = []
    self.similar_diags[(e1, e2)].append(D)

  def print_similar_diags(self, coef):
    if not self.has_ref_info:
      return
    res = []
    for i in range(coef):
      res.append([])
    for (e1, e2), list_of_D in self.similar_diags.items():
      list_of_D.sort()
      for i in range(len(list_of_D)):
        for j in range(1, coef + 1):
          """while begin >= 0 and list_of_D[i] - list_of_D[begin] == j:
            res[j-1].append((e1,e2, list_of_D[i], list_of_D[begin]))
            begin = begin -1"""
          end = i+1
          while end < len(list_of_D) and list_of_D[end] - list_of_D[i] < j:
            end += 1
          while end < len(list_of_D) and list_of_D[end] - list_of_D[i] == j:
            res[j-1].append((e1,e2, list_of_D[i], list_of_D[end]))
            end += 1
    self.logger.info(res)
    for i in range(coef):
      self.logger.info(str(i + 1) + " " + str(len(res[i])))


  def is_true_diagonal(self, diag):
    rect = diag.rectangle
    e1 = rect.e1
    e2 = rect.e2
    D = diag.D
    if e1.eid == e2.eid:
      return True
    e1_ref_infos = self.__get_ref_info(e1.eid)
    e2_ref_infos = self.__get_ref_info(e2.eid)
    if len(e1_ref_infos) == 0 or len(e2_ref_infos) == 0:
      self.unaligned += 1
      #return False
      return True
    true_diag = False
    for e1_ref_info in e1_ref_infos:
      for e2_ref_info in e2_ref_infos:
        dist = abs(e2_ref_info - e1_ref_info)
        if dist == D or dist == D + 1 or dist == D - 1 or dist == D + 2 or dist == D - 2:
          true_diag = True
          break
      if true_diag:
        break
    if true_diag:
      self.true_diags += 1
    else:
      self.not_true_diags += 1
    ###print "true_diag", true_diag
    """if not true_diag:
      print "edges", e1.eid, e2.eid
      print "e1",  e1_ref_infos
      print "e2", e2_ref_infos
      print "D", D

    print "diag support", diag.support(), "len", diag.offsetc- diag.offseta"""
    return true_diag

  def should_join(self, diag1, diag2):
     rect1 = diag1.rectangle
     rect2 = diag2.rectangle
     e1_1 = rect1.e1
     e1_2 = rect1.e2
     e2_1 = rect2.e1
     e2_2 = rect2.e2
     offseta_2 = diag2.offseta
     offsetb_2 = diag2.offsetb
     offsetc_1 = diag1.offsetc
     offsetd_1 = diag1.offsetd
     if len(self.__get_ref_info(e1_1.eid)) == 0 or len(self.__get_ref_info(e1_2.eid)) == 0 or len(self.__get_ref_info(e2_1.eid) ) == 0 or len(self.__get_ref_info(e2_2.eid)) == 0:
      self.join_unaligned += 1
      return True

     for e11_ref_info in self.__get_ref_info(e1_1.eid):
      for e12_ref_info in self.__get_ref_info(e1_2.eid):
        for e21_ref_info in self.__get_ref_info(e2_1.eid):
          for e22_ref_info in self.__get_ref_info(e2_2.eid):
            if e11_ref_info*e21_ref_info < 0 or e12_ref_info * e22_ref_info < 0:
             # print "join different strand"
              continue
            #print "join",  e11_ref_info - e21_ref_info, abs(e11_ref_info) + offsetc_1 -(  abs(e21_ref_info) + offseta_2), e12_ref_info - e22_ref_info, abs(e12_ref_info) + offsetd_1 -(abs(e22_ref_info) + offsetb_2)
            if abs(e11_ref_info) + offsetc_1 == abs(e21_ref_info) + offseta_2 and abs(e12_ref_info) + offsetd_1 == abs(e22_ref_info) + offsetb_2:
              return True
     #print offsetc_1, offseta_2, offsetd_1, offsetb_2
     #print self.__get_ref_info(e1_1.eid)
     #print self.__get_ref_info(e1_2.eid)
     #print self.__get_ref_info(e2_1.eid)
     #print self.__get_ref_info(e2_2.eid)
     return False
