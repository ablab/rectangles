import os
import saveparser
import experimental
import graph
import logging
import cStringIO
import itertools
from test_util import TestUtils
from rectangle import Rectangle  
from utils import conjugate
import utils

class BVertex(object):
  bvid = 0

  def __init__(self, key):
    self.key = key
    self.inn = []
    self.out = []
    self.bvid = BVertex.bvid
    BVertex.bvid += 1

class BEdge(object):
  beid = 0

  def __init__(self, bv1, bv2, diag):
    self.bv1, self.bv2 = bv1, bv2
    bv1.out.append(self)
    bv2.inn.append(self)
    self.diagonals = [diag]
    self.beid = BEdge.beid
    BEdge.beid += 1

  def __hash__(self):
    return self.beid

  def get_seq(self, K, d):
    (seq1, seq2) = self.get_paired_seq(K, d)
    seq = utils.seq_join(seq1, seq2).strip('N')
    return seq
  
  def get_paired_seq(self, K, d):
    seq1 = cStringIO.StringIO()
    seq2 = cStringIO.StringIO()
    seq2.write('N' * d)
    for this in self.diagonals:
      seq1.write(this.rectangle.e1.seq[this.offseta : this.offsetc])
      seq2.write(this.rectangle.e2.seq[this.offsetb : this.offsetd])
    last = self.diagonals[-1]
    seq1.write(last.rectangle.e1.seq[last.offsetc : last.offsetc + K])
    seq2.write(last.rectangle.e2.seq[last.offsetd : last.offsetd + K])
    seq1.write('N' * (len(seq2.getvalue())-len(seq1.getvalue())))
    return (seq1.getvalue(), seq2.getvalue())

  def get_cvr(self):
    cvr = 0.0
    sumlen = 0
    for this in self.diagonals:
      thiscvr = (this.rectangle.e1.cvr + this.rectangle.e2.cvr) * 0.5
      l =  this.offsetc - this.offseta
      cvr += thiscvr * l
      sumlen += l
    cvr /= sumlen
    return cvr

  def __repr__(self):
    return str((self.beid, self.diagonals))


class BGraph(object):
  def __init__(self, graph, d, test_utils):
    self.logger = logging.getLogger('rectangles')
    self.graph = graph
    self.d = d
    self.bvs = {} # key -> BVertex
    self.diagonals = set()
    self.bes = set() # BEdges
    self.test_utils = test_utils
  
  def __remove_bedge__(self, bedge):
    bv1 = bedge.bv1
    bv2 = bedge.bv2
    if bedge in bv1.out:
      bv1.out.remove(bedge)
    if (bedge in bv2.inn):
      bv2.inn.remove(bedge)
    self.__try_delete_bv(bv1)
    self.__try_delete_bv(bv2)
    if bedge in self.bes:
      self.bes.remove(bedge)
  def __try_delete_bv(self, bv):
    if len(bv.out) == 0 and len(bv.inn) == 0 and bv.key in self.bvs:
      del self.bvs[bv.key]
  
  def check_tips(self, K ):
    bv1s = set()
    bv2s = set()
    tips = set()
    for bv in self.bvs.itervalues():
        if len(bv.inn) == 1 and len(bv.out) == 0 and len(bv.inn[0].get_seq(K, self.d)) < 3 * self.d and bv.inn[0].bv1.bvid != bv.bvid:
          edge =  bv.inn[0]
          if len(edge.diagonals) == 1:
            rect = edge.diagonals[0].rectangle
            if (rect.e1.eid == rect.e2.eid):
              continue
          bv1s.add(bv)
          supp = 0
          for diag in edge.diagonals:
            supp += diag.support()
          if supp < 1000:
            tips.add(bv.inn[0])
    self.delete_tips(K,tips)

  def delete_tips(self, K, bes):
    for be in bes:
      self.__remove_bedge__(be)
      if (be.beid != be.conj.beid):
        self.__remove_bedge__(be.conj)
    self.condense()
  
  #L - big contig's len > L
  #threshold - max deep in bfs
  def delete_loops(self, K, L, threshold):
    edges_to_delete = set()
    connected_paths = set()
    for be in self.bes:
      if len(be.get_seq(K, self.d)) > L:
        self.__find_loops(be, K, threshold, L, edges_to_delete, connected_paths)
    for edge in edges_to_delete:
      if edge.beid not in connected_paths and edge.conj.beid not in connected_paths:
        self.__remove_bedge__(edge)
        self.__remove_bedge__(edge.conj)
    
  def __find_loops(self, be, K, threshold, L, edges_to_delete, connected_path):
    bv1 = be.bv2
    lst = [(bv1,0)]
    long_ends = set()
    visited_bvs = set()
    while len(lst) != 0:
      (bv, deep) = lst.pop(0)
      if deep > threshold:
        continue
      for be1 in bv.out:
        if len(be1.get_seq(K, self.d)) > L:
          long_ends.add(be1)
        else:
          visited_bvs.add(bv)
          lst.append((be1.bv2, deep + 1))
    if len(long_ends) != 1:
      return
    isoleted_component = True
    long_end = None
    for e in long_ends:
      long_end = e
    end_long_edge_id = long_end.bv2.bvid
    begin_long_edge_id = be.bv1.bvid
    for bv in visited_bvs:
      for bv2 in [e.bv2 for e in bv.out]:
        if bv2.bvid != end_long_edge_id and bv2 not in visited_bvs:
          isoleted_component = False
          return
      for bv_begin in [e.bv1 for e in bv.inn]:
        if bv_begin.bvid != begin_long_edge_id and bv_begin not in visited_bvs:
          isoleted_component = False
          return
    for bv in visited_bvs:
      if not self.__is_connected(bv, long_end, threshold):
        return
    for bv in visited_bvs:
      for e in bv.out:
        if e.bv2.bvid != end_long_edge_id:
          edges_to_delete.add(e)
      for e in bv.inn:
        bv_begin = e.bv1
        if bv_begin.bvid != begin_long_edge_id:
          edges_to_delete.add(e)
    path = self.__get_path(be.bv2, long_end, threshold) 
    for e in path:
      connected_path.add(e.beid)
    
  def __is_connected(self, bv1, be1, threshold):
    lst = [(bv1,0)]
    while len(lst) != 0:
      (bv, deep) = lst.pop(0)
      if deep > threshold:
        continue
      for be in bv.out:
        if be.beid == be1.beid:
          return True
        lst.append((be.bv2, deep + 1))
    return True
  def __get_path(self, bv1, be1, threshold):
    lst = [(bv1,0, [])]
    while len(lst) != 0:
      (bv, deep, path) = lst.pop(0)
      if deep > threshold:
        continue
      for be in bv.out:
        if be.beid == be1.beid:
          return path
        new_path = list(path)
        new_path.append(be)
        lst.append((be.bv2, deep + 1, new_path))
    return True

  def save(self, outpath, K):
    eid = 0
    left = open(os.path.join(outpath,"rectangle_paired_info_1.fasta"), "w")
    right = open(os.path.join(outpath, "rectangle_paired_info_2.fasta"), "w")
    for be in self.bes:
      (seq1, seq2) = be.get_paired_seq(K, self.d)
      seq2 = seq2.strip('N')
      seq1 = seq1.strip('N')
      left.write(">" + str(eid) + "/1\n" + seq1 + "\n")
      right.write(">" + str(eid) + "/2\n" + seq2 + "\n")
      eid += 1
    left.close()
    right.close()

  def __get_bvertex(self, key):
    if key in self.bvs:
      bv = self.bvs.pop(key)
      bv.key = bv.key.join_with(key) # transitive closure
    else:
      bv = BVertex(key)
    self.bvs[bv.key] = bv
    return bv

  def __is_bvertex(self, key):
    return key in self.bvs

  def __add_bedge(self, diag):
    bv1 = self.__get_bvertex(diag.key1)
    bv2 = self.__get_bvertex(diag.key2)
    be = BEdge(bv1, bv2, diag)
    self.bes.add(be)
    return be

  def add_diagonal(self, diag):
    if diag in self.diagonals:
      return
    be = self.__add_bedge(diag)
    conj = self.__add_bedge(diag.conj) if diag.conj != diag else be
    self.diagonals.add(diag)
    self.diagonals.add(diag.conj)
    conjugate(be, conj)
    conjugate(be.bv1, conj.bv2)
    conjugate(be.bv2, conj.bv1)
    
  def add_diagonal_and_conj(self, diag):
    for old_diag in self.diagonals:
      if diag.rectangle.e1 == old_diag.rectangle.e1 and diag.rectangle.e2 == old_diag.rectangle.e2:
        if diag.D == old_diag.D:
          return
    rect = diag.rectangle
    rect_conj = Rectangle(rect.e2.conj, rect.e1.conj)
    conjugate(rect, rect_conj)
        
    D = diag.D - diag.rectangle.e1.len + diag.rectangle.e2.len       
    pathset = diag.pathset.conj() if experimental.filter == experimental.Filter.pathsets else None
    rect_conj.add_diagonal(self.d, D, pathset)
    diag_conj = rect.conj.diagonals[D, pathset]       
    conjugate(diag, diag_conj)
        
    self.add_diagonal(diag)
 
  def __join_biedges(self, be1, be2):
        ## u ---be1---> v ---be2---> w
        ## z <--be4---- y <--be3---- x
        ## transforms to:
        ## u --------beA--------> w
        ## z <-------beB--------- x
        be3 = be2.conj
        be4 = be1.conj
        u, v, w = be1.bv1, be1.bv2, be2.bv2
        x, y, z = be3.bv1, be3.bv2, be4.bv2
        assert be1.bv2 == be2.bv1
        assert 1 == len(v.inn) == len(v.out) == len(y.inn) == len(y.out), (len(v.inn), len(v.out), len(y.inn), len(y.out))
        assert be1 != be3, "=> (v == y) => (in-degree(v) > 1)"
        assert be2 != be4, "=> (v == y) => (out-degree(v) > 1)"

        if be1 == be4 and be2 == be3:
            assert z == v == x
            assert u == y == w
            assert False
            return # TODO: think how to condense better, rare case

        if be2 == be3: # loop on the right: be1->be2=be3->be4
            assert v == x
            assert y == w
            beA = BEdge(u, z, None)
            beA.diagonals = be1.diagonals + be2.diagonals + be4.diagonals
            first_connect =  self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
            second_connect =  self.test_utils.should_join(be2.diagonals[-1],be4.diagonals[0]) 
            if first_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect += 1
            if second_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect +=1
            #print "connect diagonals be2==be3",first_connect,second_connect            
            conjugate(beA, beA)
            self.bes.add(beA)
            u.out.remove(be1)
            w.inn.remove(be2)
            z.inn.remove(be4)
            self.bes.remove(be1)
            self.bes.remove(be2)
            self.bes.remove(be4)
        elif be1 == be4: # loop on the left: be3->be1=be4->be2
            assert u == y
            assert z == v
            beA = BEdge(x, w, None)
            beA.diagonals = be3.diagonals + be1.diagonals + be2.diagonals
            first_connect =  self.test_utils.should_join(be3.diagonals[-1], be1.diagonals[0])
            second_connect =  self.test_utils.should_join(be1.diagonals[-1],be2.diagonals[0]) 
            if first_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect += 1
            if second_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect +=1
            #print "connect diagonals be1==be4",first_connect,second_connect            
            conjugate(beA, beA)
            self.bes.add(beA)
            u.out.remove(be1)
            w.inn.remove(be2)
            x.out.remove(be3)
            self.bes.remove(be1)
            self.bes.remove(be2)
            self.bes.remove(be3)
        else: # most usual case
            assert len({be1, be2, be3, be4}) == 4, (be1, be2, be3, be4) # all different
            if u == w:
                assert z == x
                assert len({u, v, w, x, y, z}) == 4, (u, v, w, x, y, z) # same ends, ok
            elif u == x:
                assert z == w
                assert len({u, v, w, x, y, z}) == 4, (u, v, w, x, y, z) # conjugated ends, ok
            else:
                assert len({u, v, w, x, y, z}) == 6, (u, v, w, x, y, z) # all different
            # TODO: check (x == u and w == z)
            beA = BEdge(u, w, None)
            beA.diagonals = be1.diagonals + be2.diagonals
            first_connect =  self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
            #print "connect diagonals usual be1, be2",first_connect 
            #if not first_connect:
              #print "shouldn't join", be1.diagonals[-1].support(), be2.diagonals[0].support()
            #else:
              #print "should",  be1.diagonals[-1].support(), be2.diagonals[0].support()
            second_connect =  self.test_utils.should_join(be3.diagonals[-1],be4.diagonals[0]) 
            if first_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect += 1
            if second_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect +=1
            beB = BEdge(x, z, None)
            beB.diagonals = be3.diagonals + be4.diagonals
            #print "connect diagonals usual be3, be4", second_connect
            #if not second_connect:
              #print "shouldn't join", be3.diagonals[-1].support(), be4.diagonals[0].support()
            #else:
              #print "should", be3.diagonals[-1].support(), be4.diagonals[0].support()
            
            conjugate(beA, beB)
            self.bes.add(beA)
            self.bes.add(beB)
            u.out.remove(be1)
            w.inn.remove(be2)
            x.out.remove(be3)
            z.inn.remove(be4)
            self.bes.remove(be1)
            self.bes.remove(be2)
            self.bes.remove(be3)
            self.bes.remove(be4)

        v.inn, v.out = [], []
        y.inn, y.out = [], []
        self.bvs.pop(v.key)
        self.bvs.pop(y.key)
 
  def condense(self):
        l = len(self.bvs)
        for bv in self.bvs.values(): # copy because can be: "Set changed size during iteration"
            if len(bv.inn) == 1 and len(bv.out) == 1 and (bv.inn[0] != bv.out[0]):
                self.__join_biedges(bv.inn[0], bv.out[0])
                self.__check()
        self.logger.info("Condensed %d bi-vertices (out of %d). %d bi-vertices left." % (l - len(self.bvs), l,
                                                                                         len(self.bvs)))

  def project(self, outpath):
        g = graph.Graph()
        for be in self.bes:
            # v ---------be--------> w
            # y <-----be.conj------- x
            v, w = be.bv1, be.bv2
            x, y = be.conj.bv1, be.conj.bv2
            #assert be != be.conj
            #assert v != w
            #assert v != x
            seq = be.get_seq(self.graph.K, self.d)
            cvr = be.get_cvr()
            g.add_vertex(v.bvid, y.bvid)
            g.add_vertex(y.bvid, v.bvid)
            g.add_vertex(w.bvid, x.bvid)
            g.add_vertex(x.bvid, w.bvid)
            g.add_edge(be.beid, v.bvid, w.bvid, len(seq) - self.graph.K, be.conj.beid)
            g.add_seq(be.beid, seq)
            g.add_cvr(be.beid, cvr)
        g.update_K()
        maxv = BVertex.bvid
        maxe = BEdge.beid
        taken = set()
        for diag in self.diagonals:
            taken.add(diag.rectangle.e1)
            taken.add(diag.rectangle.e2)
        for e in self.graph.es.itervalues():
            if e not in taken:
                # v ---e1---> w
                # x <--e2---- y
                assert e.conj not in taken
                e1 = e
                e2 = e.conj
                v = e1.v1.vid + maxv
                w = e1.v2.vid + maxv
                y = e2.v1.vid + maxv
                x = e2.v2.vid + maxv
                g.add_vertex(v, x)
                g.add_vertex(x, v)
                g.add_vertex(w, y)
                g.add_vertex(y, w)
                seq = e1.seq
                g.add_edge(e1.eid + maxe, v, w, len(seq) - self.graph.K, e2.eid + maxe)
                g.add_seq(e1.eid + maxe, seq)
        return g


  def __check(self):
        for edge in self.bes:
            for this, next in itertools.izip(edge.diagonals, edge.diagonals[1:]):
                assert this.key2 == next.key1, (this, "->", next)

  def build_missing_rectangles(self, K, rectangles):
    return
    threshold = self.d
    self.test_utils.logger.info( "treshold " + str( threshold))
    bv1s = set()
    bv2s = set()
    for bv in self.bvs.itervalues():
        if len(bv.inn) == 1 and len(bv.out) == 0 and len(bv.inn[0].get_seq(K, self.d)) > 3 * self.d:
          bv1s.add(bv)
        if len(bv.inn) == 0 and len(bv.out) == 1 and len(bv.out[0].get_seq(K, self.d)) > 3 * self.d:
          bv2s.add(bv)
    assert len(bv1s) == len(bv2s) # because of rev-compl
    self.test_utils.logger.info("bv1s.len "+ str( len(bv1s)))
         
    all_paired_paths = []
    for bv1 in bv1s:
      be1 = bv1.inn[0]
      diag1 = be1.diagonals[-1]
      for bv2 in bv2s:
        be2 = bv2.out[0]
        if (be1.beid == be2.beid):
          continue
        diag2 = be2.diagonals[0]
        paths1 = graph.find_paths(diag1.rectangle.e1.v1, diag2.rectangle.e1.v1, diag1.rectangle.e1, threshold + diag1.rectangle.e1.len)
        paths2 = graph.find_paths(diag1.rectangle.e2.v1, diag2.rectangle.e2.v1, diag1.rectangle.e2, threshold + diag1.rectangle.e2.len)
        paired_paths = find_pair_paths(paths1, paths2, diag1, diag2)
        if len(paired_paths) != 0:
          all_paired_paths.append((paired_paths, diag1, diag2))
    can_find_one_path_more = True
    added_paths = []
    while can_find_one_path_more:
      the_best_path = None
      the_best_support = 0
      can_find_one_path_more = False

      for paired_paths in all_paired_paths:
        (best_support, best_len, best_rectangles, best_diags, best_path) = self.choose_best_path(paired_paths[0], rectangles, paired_paths[1], paired_paths[2], self.d, added_paths)
        if best_support/best_len > 0.001 and best_support/best_len > the_best_support:
          the_best_support = best_support/best_len
          the_best_path = (best_support, best_len, best_rectangles, best_diags, best_path)
      if the_best_path:
        added_paths.append(the_best_path[-1])
        (best_support, best_len, best_rectangles, best_diags, best_path) = the_best_path
        can_find_one_path_more = True
        prev_diag = best_diags[0]
        


  def choose_best_path(self, paired_paths, rectangeles_set, diag1, diag2, d, added_paths):
      best_support = 0
      best_len = 10000
      best_rectangles = []
      best_diags = []
      best_path = paired_paths[0]
      best_not_supported = 0
        
      for paired_path in paired_paths:
        (path1, path2, path_len) = paired_path
        if paired_path in added_paths:
          continue
        first_shift = diag1.offseta
        second_shift = diag1.offsetb
        path1.append(diag2.rectangle.e1)
        path2.append(diag2.rectangle.e2)
        rectangles = []
        diags = []
        not_supported = []
        path_support = 0
        pos_first_path = 0
        pos_second_path = 0
        first_len = first_shift
        make_less_N50 = False
        while not make_less_N50 and first_len < path_len + diag2.offseta:
          ed1 = path1[pos_first_path]
          ed2 = path2[pos_second_path]
          rectangle = Rectangle(ed1,ed2)
          rectangle.add_diagonal(d, d + first_shift - second_shift)
          rect_diag = rectangle.get_closest_diagonal(d + first_shift - second_shift) 
          if (not (rect_diag.key1 == diag1.key1 and rect_diag.key2 == diag1.key2) and not(rect_diag.key1 == diag2.key1 and rect_diag.key2 == diag2.key2)):
            can_use = [diag1.key1, diag1.key2, diag2.key1, diag2.key2]
            if (rect_diag.key1 in self.bvs and rect_diag.key1 not in can_use) or  (rect_diag.key2 in self.bvs and rect_diag.key2 not in can_use):
              make_less_N50 = True
              continue
          diags.append(rect_diag)
          rectangles.append(rectangle)
          rectangeles_set.use_prd_diag(rect_diag)
          if rect_diag.prd_support < 0.00001 and (ed2.len > 10 and ed1.len > 10):
            make_less_N50 = True
            continue
          path_support += rect_diag.prd_support 
          if ed2.len - second_shift < ed1.len - first_shift:
            pos_second_path += 1
            first_shift += ed2.len - second_shift
            first_len += ed2.len - second_shift
            if rect_diag.prd_support < 0.000000001:
              not_supported.append(ed2.len - second_shift)
            second_shift = 0
          elif ed1.len - first_shift < ed2.len - second_shift:
            pos_first_path += 1
            first_len += ed1.len - first_shift
            second_shift += ed1.len - first_shift
            if rect_diag.prd_support < 0.000000001:
              not_supported.append(ed1.len - first_shift)
            first_shift = 0
          else:
            first_len += ed1.len - first_shift
            pos_second_path += 1
            pos_first_path += 1
            first_shift = 0
            second_shift = 0
        if not make_less_N50 and path_len > 1 and  path_support / path_len < 1000 and  path_support / path_len > best_support:
          best_support = path_support 
          best_len = path_len
          best_rectangles = rectangles
          best_diags = diags
          best_path = (path1, path2, path_len)
          best_not_supported = not_supported
      return (best_support,best_len, best_rectangles, best_diags, best_path)


def find_pair_paths(paths1, paths2, diag1, diag2):
  paired_paths = []
  for (path1, len1) in paths1:
    for (path2, len2) in paths2:
      if path1[0] != diag1.rectangle.e1 or path2[0] != diag1.rectangle.e2:
        continue
      if len1 - diag1.offseta + diag2.offseta == len2 - diag1.offsetb + diag2.offsetb: 
        paired_paths.append((path1, path2, len1))
  return paired_paths
  
