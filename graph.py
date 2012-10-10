import logging
import N50
import sys
import utils
import saveparser
from utils import conjugate

##################
# DEBRUIJN GRAPH #
##################
alphabet = "ACGT"

def findOutEdges(vertexBody, kmerMap):
    nextKmer = [(vertexBody + base) for base in alphabet]
    return [kmer for kmer in nextKmer if kmer in kmerMap]
 
def findInEdges(vertexBody, kmerMap):
    nextKmer = [(base + vertexBody) for base in alphabet]
    return [kmer for kmer in nextKmer if kmer in kmerMap]

def extendForward(vertexBody, kmerMap):    
    inEdge = findInEdges(vertexBody, kmerMap)
    outEdge = findOutEdges(vertexBody, kmerMap)
    if len(inEdge) == 1 and len(outEdge) == 1:
        return outEdge[0]
    return None

def extendBackward(vertexBody, kmerMap):    
    inEdge = findInEdges(vertexBody, kmerMap)
    outEdge = findOutEdges(vertexBody, kmerMap)
    if len(inEdge) == 1 and len(outEdge) == 1:
        return inEdge[0]
    return None

 
class Graph(object):
    class Vertex(object):
        def __init__(self, vid, conj):
            self.vid = vid
            self.inn = []
            self.out = []
            conjugate(self, conj)

        def __hash__(self):
            return self.vid

        def seq(self, K):
            if self.out:
                return self.out[0].seq[ :K]
            else:
                return self.inn[0].seq[-K:]

    def __repr__(self):
            return "V%d" % (self.vid)

    class Edge(object):
        def __init__(self, eid, v1, v2, edge_len, conj):
            self.eid = eid
            self.v1 = v1
            self.v2 = v2
            self.len = edge_len
            conjugate(self, conj)
            self.seq = None
            self.cvr = 0

        def __hash__(self):
            return self.eid

        def __repr__(self):
            return "E%d(%d)" % (self.eid, self.len)

    def __init__(self):
        self.vs = {} # vid -> Vertex
        self.es = {} # eid -> Edge
        self.max_eid = 0
        self.etalon_dist = dict()
        self.logger = logging.getLogger('rectangles')

    def add_vertex(self, vid, conj_id):
        assert vid != conj_id, "Vertex can't be self-conjugated"
        conj = self.vs.get(conj_id, None)
        v = Graph.Vertex(vid, conj)
        self.vs[vid] = v
        return v

    def add_edge(self, eid, v1id, v2id, edge_len, conj_id):
        #assert eid != conj_id, "Self-conjugate edges are not supported yet"
        if eid in self.es:
          return self.es[eid]
        if eid > self.max_eid or conj_id > self.max_eid:
            self.max_eid = max(eid, conj_id)
        conj = self.es.get(conj_id, None)
        v1 = self.vs[v1id]
        v2 = self.vs[v2id]
        e = Graph.Edge(eid, v1, v2, edge_len, conj)
        v1.out.append(e)
        v2.inn.append(e)
        self.es[eid] = e
        if eid == conj_id:
          conjugate(e, e)
        return e

    def add_seq(self, eid, seq):
        self.es[eid].seq = seq

    def add_cvr(self, eid, cvr):
        self.es[eid].cvr = cvr

    def get_edge(self, eid):
        return self.es[eid]

    def update_K(self):
        assert len(self.es) > 0, "Empty graph"
        any_edge = (e for e in self.es.itervalues()).next()
        K = len(any_edge.seq) - any_edge.len
        self.K = K

    def check(self):
        for v in self.vs.itervalues():
            assert v.conj, "Some vertex have no conjugate"
        for e in self.es.itervalues():
            assert e.conj, "Some edge have no conjugate"
        for e in self.es.itervalues():
            assert self.K == len(e.seq) - e.len, "Inconsistent K"
        for e in self.es.itervalues():
            assert e.seq == utils.rc(e.conj.seq), (e.seq, utils.rc(e.conj.seq))


    def fasta(self, stream=sys.stdout, is_all=False):
        contig_id = 0 
        for edge in self.es.itervalues():
            if is_all or edge.conj.eid <= edge.eid: # non-conjugate
                #if 'N' in edge.seq: continue # hack - remove small edges with Ns :)
                for id_contig, contig in enumerate(edge.seq.split('N')):
                    if not contig: continue
                    l = len(contig)
                    print >>stream, '>contig_%d_%d_%d_%d_l=%06d' % (contig_id, edge.eid, edge.conj.eid, id_contig, l)
                    contig_id += 1
                    for l in xrange(0, l, 60):
                        print >>stream, contig[l:l + 60]

    def stats(self, d):
        ls = []
        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
                ls.append(len(edge.seq))
        ls.sort()
        self.logger.info('Edges   = %s' % ls)
        self.logger.info('#Edges  = %d' % len(ls))
        self.logger.info('Edg N50 = %dbp' % N50.N50(ls))
        ls = []
        for edge in self.es.itervalues():
            #if 'N' in edge.seq: continue
            if edge.conj.eid <= edge.eid: # non-conjugate
                lens = filter(None, map(len, edge.seq.split('N')))
                #lens = [len(edge.seq)]
                ls += lens
        ls.sort()
        self.logger.info('Contigs = %s' % ls)
        self.logger.info('#Contigs= %d' % len(ls))
        self.logger.info('K       = %d' % self.K)
        self.logger.info('d       = %d' % d)
        self.logger.info('Total   = %dbp' % sum(ls))
        self.logger.info('N50     = %d' % N50.N50(ls))
        self.logger.info('We split small edges (with Ns) to multiple contigs, so #edges < #contigs')
        self.logger.info('N50>1000 = %d' % N50.N50([x for x in ls if x >= 1000]))
        return N50.N50(ls)

    def save(self, filename):
        # Graph save
        grp = open(filename + '.grp', 'w')
        print >>grp, len(self.vs), len(self.es)
        for vertex in self.vs.itervalues(): # Vertex 2 ~ 1 .
            print >>grp, 'Vertex %d ~ %d .' % (vertex.vid, vertex.conj.vid)#, vertex.inn, vertex.out, vertex.seq(15)
        print >>grp  # empty line
        for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
            print >>grp, 'Edge %d : %d -> %d, l = %d ~ %d .' % (
                edge.eid, edge.v1.vid, edge.v2.vid, len(edge.seq), edge.conj.eid)#, self.etalon_dist[edge.eid], edge.seq
        grp.close()
  
        # Sequences save
        sqn = open(filename + '.sqn', 'w')
        contigs = open(filename + ".contigs", "w")
        print >>sqn, len(self.es)
        visited = set()
        for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
            if edge.eid not in visited:
              visited.add(edge.eid)
              visited.add(edge.conj.eid)
              print >> contigs, ">" + str(edge.eid) + "\n" + edge.seq.strip()
            print >>sqn, '%d %s .' % (edge.eid, edge.seq)
        sqn.close()
        contigs.close()
        # Coverage save
        cvr = open(filename + '.cvr', 'w')
        print >>cvr, len(self.es)
        for edge in self.es.itervalues():  # 15 1.234 .
            print >>cvr, '%d %f .' % (edge.eid, edge.cvr)
        cvr.close()
        # Stupid pos file
        pos = open(filename + '.pos', 'w')
        print >>pos, len(self.es)
        for edge in self.es.itervalues():  # 15 0
            print >>pos, '%d %d' % (edge.eid, 0)
        pos.close()

    def load(self, grp_filename, sqn_filename, cvr_filename):
        for vid, conj in saveparser.grp_vertices(grp_filename):
            self.add_vertex(vid, conj)
        for eid, v1id, v2id, l, conj in saveparser.grp_edges(grp_filename):
            self.add_edge(eid, v1id, v2id, l, conj)
        for eid, seq in saveparser.sqn(sqn_filename):
            self.add_seq(eid, seq)
        for eid, cvr in saveparser.cvr(cvr_filename):
            self.add_cvr(eid, cvr)
        self.update_K()

    
   
    def make_graph(self, genome,k):
        #print genome,"\n", utils.rc(genome)
        self.K = k
        kmers = dict()
        for i in range(len(genome) - k + 1):
          kmer = genome[i: i +k]
          if kmer not in kmers:
            kmers[kmer] = [i+1]
            kmers[utils.rc(kmer)] = [-(len(genome) - i -k + 1)]
          else:
            kmers[kmer].append(i+1)
            kmers[utils.rc(kmer)].append(-(len(genome)- i - k+1))

        visit = set()
        vid = 0
        eid = 0
        edges = set()
        verts = dict()
        for key in kmers:
            if key in visit:
                continue
            body = [key[-1]]
            endVertex = key[1:]                    
            while True:
                nextKmer = extendForward(endVertex, kmers)
                if nextKmer == None:
                    break
                body.append(nextKmer[-1])
                endVertex = nextKmer[1:]
                visit.add(nextKmer)
                visit.add(utils.rc(nextKmer))
                
            beginVertex = key[:-1]
            while True:
                nextKmer = extendBackward(beginVertex, kmers)
                if nextKmer == None:
                    break
                body.insert(0, nextKmer[-1])
                beginVertex = nextKmer[0:-1]
                visit.add(nextKmer)
                visit.add(utils.rc(nextKmer))
                
            body = beginVertex + ''.join(body)
            if beginVertex not in verts:
              beginRef = self.add_vertex(vid, vid+1)
              r_endRef = self.add_vertex(vid+1, vid)
              verts[beginVertex] = beginRef.vid
              verts[utils.rc(beginVertex)] = r_endRef.vid
              vid +=2
            if endVertex not in verts:
							endRef = self.add_vertex(vid, vid+1)
							r_beginRef = self.add_vertex(vid+1, vid)
							verts[endVertex] = endRef.vid
							verts[utils.rc(endVertex)] = r_beginRef.vid
							vid +=2
            bv = verts[beginVertex]
            ev = verts[endVertex]
            rbv = verts[utils.rc(endVertex)]
            rev = verts[utils.rc(beginVertex)]
            if (bv, ev) not in edges:
              if (bv,ev) == (rbv, rev) and body == utils.rc(body):
                self.add_edge(eid, bv, ev, len(body) -k +1 , eid)
                edges.add((bv,ev))
                self.add_seq(eid, body)
                #print kmers[body[:k]], kmers[utils.rc(body)[:k]]
                self.etalon_dist[eid] = kmers[body[:k]] + kmers[utils.rc(body)[:k]]
                eid += 1
              else:
                self.add_edge(eid, bv, ev, len(body) - k + 1, eid +1)
                self.add_edge(eid +1, rbv, rev, len(body) -k +1, eid)
                edges.add((bv,ev)) 
                edges.add((rbv, rev))
                self.add_seq(eid, body)
                self.add_seq(eid +1, utils.rc(body))
                #print "edge\n", body,"\n", utils.rc(body),"\n", beginVertex,"\n", endVertex, "\n",utils.rc(endVertex), "\n", utils.rc(beginVertex)

                #print kmers[body[:k]],  kmers[utils.rc(body)[:k]]
                self.etalon_dist[eid] = kmers[body[:k]]
                self.etalon_dist[eid+1] = kmers[utils.rc(body)[:k]]
                eid += 2
				
			
 
   
    def dfs_with_etalon_dist(self, e, d):
      limit1 = d - e.len
      limit2 = d
      if e.len > d:
        yield e, 0
      ls = [set() for _ in xrange(limit2)]
      ls[0].add(e.v2)
      all_dist  = dict()
      all_dist[(0,e.v2.vid)] = (e, self.etalon_dist[e.eid])
      for pos in xrange(limit2):
        for v in ls[pos]:
          (prev_e, dists) = all_dist[(pos, v.vid)]
          for e2 in v.out:
            new_dists = []
            for dist in dists:
              if dist>= 0 and dist + prev_e.len  + 1 in self.etalon_dist[e2.eid]:
                new_dists.append(dist + prev_e.len + 1)
              elif dist<= 0 and -(-dist + prev_e.len  + 1) in self.etalon_dist[e2.eid]:
                new_dists.append(dist - prev_e.len -1 )
            if len(new_dists) == 0:
              continue

            pos2 = pos + e2.len
            if pos2 < limit2:
              ls[pos2].add(e2.v2)
              if (pos2, e2.v2.vid) in all_dist:
                new_dists += all_dist[(pos2, e2.v2.vid)][1]
              all_dist[(pos2, e2.v2.vid)] = (e2, new_dists)
            if pos + e2.len > limit1:
              yield e2, pos + e.len



    def dfs(self, e, d):
      limit1 = d - e.len
      limit2 = d
      if e.len > d:
          yield e, 0
      ls = [set() for _ in xrange(limit2)]
      ls[0].add(e.v2)
      for pos in xrange(limit2):
          for v in ls[pos]:
              #v = e.v2
              for e2 in v.out:
                pos2 = pos + e2.len
                if pos2 < limit2:
                  ls[pos2].add(e2.v2)
                if pos + e2.len > limit1:
                  yield e2, pos + e.len


def number_of_pathes(v1, v2, limit):
    limit += 1 # limit inclusive
    ls = [{} for _ in xrange(limit)]
    ls[0][v1] = 1
    for pos in xrange(limit):
        for v, cnt in ls[pos].iteritems:
            if v == v2: return cnt
            for e in v.out:
                pos2 = pos + e.len
                if pos2 < limit:
                    v2 = e.v2
                    ls[pos2][v2] = ls[pos2].get(v2, 0) + cnt
    return 0

def find_paths(v1, v2, e1, threshold):
  return __find_all_paths(v1, v2, e1, [], 0, [], threshold)

def __find_all_paths(v1, v2, e1, path, path_len, paths, threshold):
    for e in v1.out:
        if (len(path) == 0) and e != e1:
          continue
        if e.len > threshold:
            return paths
        new_path = list(path)
        new_path.append(e)
        new_path_len = path_len + e.len
        if e.v2 == v2:
            paths.append((new_path, new_path_len))
            return paths
        else:
            if (len(path) > 10):
              return paths
            __find_all_paths(e.v2, v2, e1, new_path, new_path_len, paths, threshold - e.len)
    return paths
 
      
