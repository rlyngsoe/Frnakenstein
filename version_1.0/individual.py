# The different classes used to construct individuals in the genetic
# algorithm population.
from random import uniform, randint, shuffle, choice
from math import log
from types import IntType, StringType
from copy import deepcopy
from structure import int2bp, int2str, str2int, bp2int, posdist, posdif, Structure
from mfe import fold, boltzmann, energy, design, pfold, pboltzmann, pkfold ,design
from collections import defaultdict
import operator
from time import time

### INDIVIDUAL; general class containing the core of a member of the
### GA population.
class Individual:

  def __init__(self, seq = None, target=None, dist=None, id = None, softtype = "RNAfold",ismutant=True, motif = None):
    # Initialise attributes capturing sequence and structure information
    self.seq = seq
    self.struc = {}
    self.prob = defaultdict(tuple)
    self.fit = None
    self.obj = None
    self.optscore = defaultdict(tuple)
    self.structurescore = {}
    self.ensemble = {}
    self.target=target
    self.dist=dist
    self.id = id
    self.softtype=softtype
    self.ismutant=ismutant
    self.motif = motif
  
  # Length of individual is length of its sequence
  def __len__(self):
    if self.seq == None:
      return 0
    else:
      return len(self.seq)
  # Generate string representation of individual
  def __str__(self):
    return self.getseq()
  # Retrieve current sequence of individual
  def getseq(self):
    return int2str(self.seq)
  def gettarget(self):
      return 
  # Retrieve MFE structure of current sequence at temperature T - if
  # necessary, it is computed
  def getstruct(self, T = None):    
    if not self.struc.has_key(T):
      if self.softtype == 'PPfold':
        s, opt = pfold(self.seq, T)
        self.optscore[T]=(opt,)
        self.struc[T] = Structure(s,T)
      elif self.softtype == 'pknotsRG':
        s, self.optscore[T] = pkfold(self.seq, T)
        self.struc[T] = Structure(s,T) 
      elif self.softtype == 'Both':
        s, opt = fold(self.seq, T)
        s2, opt2 = pfold(self.seq, T)
        if (opt < opt2):
          self.struc[T]=Structure(s,T)
        else:
          self.struc[T]=Structure(s2,T)
        self.optscore[T]=(opt,opt2) 
      else:
        s, opt = fold(self.seq, T)
        self.optscore[T]=(opt,)
        self.struc[T] = Structure(s,T)
      
    return self.struc[T]
  # Retrieve identifier of individual
  def getid(self):
    return self.id
  # Set identifier of individual
  def setid(self, id):
    self.id = id
  # Retrieve base pairing probabilities for current sequence at
  # temperature T - if necessary, it is computed - using specified
  # software
  def getprobability(self, T = None):
    if not self.prob.has_key(T):
      if self.softtype == 'PPfold':
        prob, s, opt = pboltzmann(self.seq, T)
        self.prob[T]=(prob,)
        self.optscore[T]=(opt,)
        self.struc[T] = Structure(s, T)
      elif self.softtype == 'Both':
        s, opt = fold(self.seq, T)
        s2, opt2 = pfold(self.seq, T)
        if (opt < opt2 ):
          prob, s_b, opt, ens = boltzmann(self.seq, T)
          prob2, s_p, opt2, ens = boltzmann(self.seq, T, s_b)
          self.struc[T] = Structure(s_b, T)
        else:
          prob2, s_p, opt2 = pboltzmann(self.seq, T)
          prob, s_b, opt, ens = boltzmann(self.seq, T, s_p)
          self.struc[T] = Structure(s_p, T)
          
        self.ensemble[T]=ens          
        self.prob[T]=(prob,prob2)
        self.optscore[T]=(opt,opt2)
      else:
        prob, s, opt, self.ensemble[T] = boltzmann(self.seq, T)
        self.prob[T]=(prob,)
        self.optscore[T]=(opt,)
        self.struc[T] = Structure(s, T)
        
    return self.prob[T]                                
  # Retrieve minimum free energy taken over all structures for current
  # sequence at temperature T - if necessary, it is computed.
  def getoptscore(self, T = None):    
    if not self.optscore.has_key(T):
    
      if self.softtype == 'PPfold':
              stuc, opt = pfold(self.seq, T)
      	      self.struc[T]=(stuc,)
      	      self.optscore[T] = (opt,)
      	      
      elif self.softtype == 'pknotsRG':
              self.struc[T], self.optscore[T] = pkfold(self.seq, T)
        
      elif self.softtype == 'Both':
              stuc, opt= fold(self.seq, T)
              stuc2, opt2 = pfold(self.seq, T)
              if (opt < opt2):
                self.struc[T] = Structure(struc, T)
              else:
                self.struc[T] = Structure(struc2, T)
                
              self.optscore[T]=(opt,opt2)
      else:
	      stuc, opt = fold(self.seq, T)
      	      self.struc[T]=(stuc,)
      	      self.optscore[T] = (opt,)

    return self.optscore[T]
  
  # Retrieve free energy for structure t at temperature T
  def getstructurescore(self, t, T = None):
    
    if not self.structurescore.has_key((t, T)):
        self.structurescore[(t, T)] = energy(self.seq, t, T)
        
    return self.structurescore[(t, T)]
  
  # Retrieve ensemble free energy at temperature T
  def getensemblescore(self, T = None):
     
    if not self.ensemble.has_key(T):
      if self.softtype == 'PPfold':
        self.prob[T], s, self.optscore[T] = pboltzmann(self.seq, T)
        self.struc[T] = Structure(s, T)
      elif self.softtype == "pknotsRG":
        self.ensemble[T] = pkfold(self.seq)[1]
      else:
        self.prob[T], s, self.optscore[T], self.ensemble[T] = boltzmann(self.seq, T)
        self.struc[T] = Structure(s, T)
        
    return self.ensemble[T]
  
  # Retrieve objective value of individual
  def getobjective(self):
    if self.obj == None:
      self.obj = self.objective()
    return self.obj
  # Retrieve fitness of individual
  def getfitness(self):
    if self.fit == None:
      self.fit = self.fitness()
    return self.fit
  
  # Retrieve set list of fitness
  def getparetofitness(self):
    return self.paretofitness()

  # Retrieve secondary structure folding software
  def getsofttype(self):
      if self.softtype==None:
          return "RNAfold"
      else:
          return self.softtype

  # Calculate "GC" base pair for given sequence      
  def getgccontent(self):
      numc=self.getseq().count('C')
      numg=self.getseq().count('G')
      return float(numc+numg)/len(self.seq)

  # Retrieve motif
  def getmotif(self):
    if not self.motif:
      return None
    else:
      return self.motif

  # Draw nucleotides for all positions in dependency structure
  # containing idx from dist and update sequence seq accordingly. If
  # position idx is defined (not None), it is assumed to be a position
  # forced to mutate, and draws are done according to a mutation,
  # including effects of current nucleotide. If not defined, draws are
  # done according to an initialisation where no current nucleotide is
  # assumed. In the latter case, all relevant positions in seq are
  # assumed to be None.
  def __drawondependencystructure__(self, dist, seq, target, idx, visited = None):
    self.motif = self.getmotif()
    # List dependencies, 
    temp = list(target.getdependencies())
    for item in temp:
      temp2 = list(item)
      if idx in temp2:
        temp2.remove(idx)
        deps = temp2
    if len(deps) == 0:
      # If motif is desired and nucleotide is specified at idx, set equal
      # to specified nucleotide
      if self.motif:
        if len(self.motif) != len(seq):
          raise ValueError, "Motif length is not equal to target length"
        for letter in self.motif:
          if letter not in "ACGU?":
            raise ValueError, "Disallowed characters in motif (must be A,C,G,U, or ?)"
        if self.motif[idx] != "?":
          seq[idx] = str2int(self.motif[idx])[0]
        else:
          seq[idx] = dist.drawunpaired(seq[idx])
      # Single unpaired position in structure and not in motif
      else:
        seq[idx] = dist.drawunpaired(seq[idx])
    else:
      if len(deps) == 1:
        # Dependency structure consists of single base pair
        if self.motif:
          if self.motif[idx] != "?":
            seq[idx] = str2int(self.motif[idx])[0]
            if self.motif[deps[0]] != "?":
              seq[deps[0]] = str2int(self.motif[deps[0]])[0]
            else:  
              # If motif specified and non-blank here, draw paired base in pairing
              # index with paired base probability as specified
              seq[deps[0]] = dist.drawpaired(neighbours = str2int(self.motif[idx]))
          else:
            if self.motif[deps[0]] != "?":
              seq[idx] = dist.drawpaired(neighbours = str2int(self.motif[deps[0]]))
            else:
              seq[idx], seq[deps[0]] = int2bp(dist.drawbasepair(seq[idx]))
        else:
          seq[idx], seq[deps[0]] = int2bp(dist.drawbasepair(seq[idx]))
      else:
        # Complex dependency structure
        # Remember whether we are mutating or drawing from scratch and consider
        # motif
        # Create list of tuples detailing the index and nucleotide (in integer form)
        # of each specified motif entry 
        pairings = [[0,0,0,1],[0,0,1,0],[0,1,0,1],[1,0,1,0]]
        possible_nucs = [[1,1,1,1] for i in xrange(len(seq))]
        motif_nucs = []
        # Assume start structure is compatible and not fixed (we check both later)
        compatible = True
        fitted = False

        # Define dot product
        def multiply(self,seq):
          for i in xrange(4):
            self[i] *= seq[i]
          return self

        # Function restricting nucleotide to two possibilities...
        def or_multiply(self,seq1,seq2):
          options1 = multiply(self,seq1)
          options2 = multiply(self,seq2)
          for i in xrange(4):
            self[i] = min(options1[i]+options2[i], 1)
          return self

        # ... and three possibilities
        def or_or_multiply(self,seq1,seq2,seq3):
          options1 = multiply(self,seq1)
          options2 = multiply(self,seq2)
          options3 = multiply(self,seq3)
          for i in xrange(4):
            self[i] = min(options1[i]+options2[i]+options3[i], 1)
          return self

        # Simple function summing column of a matrix
        def sumcolumn(i,mat):
          total = 0
          for k in xrange(len(mat[i])):
            total += mat[i][k]
          return total

        # Sum matrix
        def summatrix(mat):
          matsum = 0
          for q in range(len(mat)):
            matsum += sumcolumn(q,mat)
          return matsum

        # Simple function returning list of non-zero indices in a matrix column
        def nonzero(i,mat):
          nonzero_list = []
          for r in xrange(len(mat[i])):
            if mat[i][r] != 0:
              nonzero_list.append(r)
          return nonzero_list

        # Function returning a set of pairing elements for any index in the structure
        def index_pairs(i,targets):
          pairs = set([])
          for tar in targets:
            for j in Structure.basepairs(tar):
              # Check if left or right bracket
              if i == j[0]:
                pairs.add(j[1])
              elif i == j[1]:
                pairs.add(j[0])
          return [i,pairs]

        # Function running index_pairs over all elements
        def pairs(targets):
          pairs = []
          for i in xrange(len(seq)):
            if index_pairs(i,targets) not in pairs:
              pairs.append(index_pairs(i,targets))
          return pairs

        # Create list of tuples of specified nucleotide and index in structure
        # structure
        for i in xrange(len(seq)):
          if self.motif[i] != "?":
            motif_nucs.append([str2int(self.motif[i])[0],i])

        # Fix all specified nucleotides from motif
        for i in motif_nucs:
          for j in xrange(4):
            possible_nucs[i[1]][j] = 0
          possible_nucs[i[1]][i[0]] = 1

        # We wish to check possibilities whilst structure is compatible and not complete
        pairs_list = pairs(target)
        while compatible and not fitted:
          # Check current possibilities
          old_sum = summatrix(possible_nucs)
          new_sum = -1
          while old_sum != new_sum:
            old_sum = new_sum
            for i in pairs_list:
              choices = sumcolumn(i[0],possible_nucs)
              k = nonzero(i[0],possible_nucs)
              for j in i[1]:
                # No choice left for this index
                if choices == 0:
                  compatible == False
                  raise NameError, "Motif invalid with target structures. Non-canonical base pairs are disallowed"
                # Index fixed
                elif choices == 1:
                  possible_nucs[j] = multiply(possible_nucs[j],pairings[k[0]])
                # Two possibilities
                elif choices == 2:
                  possible_nucs[j] = or_multiply(possible_nucs[j],pairings[k[0]],pairings[k[1]])
                # Three possibilities
                elif choices == 3:
                  possible_nucs[j] = or_or_multiply(possible_nucs[j],pairings[k[0]],pairings[k[1]],pairings[k[3]])
                if sumcolumn(j,possible_nucs) == 0:
                  compatible == False
                  raise NameError, "Motif invalid with target structures. Non-canonical base pairs are disallowed"
            new_sum = summatrix(possible_nucs)

          # Check that we have no zero columns (otherwise structure is incompatible) and
          # if every column total is one (in which case we are done)
          if compatible:
            correct_count = 0
            for i in deps:
              if sumcolumn(i,possible_nucs) == 0:
                compatible = False
                raise NameError, "Motif invalid with target structures. Non-canonical base pairs are disallowed"
              elif sumcolumn(i,possible_nucs) == 1:
                correct_count += 1
                seq[i] = (nonzero(i,possible_nucs)[0])
            if correct_count == len(deps):
              fitted = True

          # If not done, randomly choose column from dependency structure and randomly fix element in it to possible
          # nucleotide
          if compatible and not fitted:
            unfixed_columns = []
            for i in deps:
              if sumcolumn(i,possible_nucs) > 1:
                unfixed_columns.append(i)
            col = choice(unfixed_columns)
            nucs = []
            for i in xrange(4):
              if possible_nucs[col][i] == 1:
                nucs.append(i)
            possible_nucs[col] = [0,0,0,0]
            possible_nucs[col][choice(nucs)] = 1

  # Function for mutating in the dependency structure containing
  # position i; if dist is specified, new bases/base pairs are chosen
  # according to this.
  def mutate(self, i, target = None, dist = None, softtype= None):
    if target == None:
      target = self.target
    if dist == None:
      dist = self.dist
    if softtype == None:
      softtype = self.softtype
      
    s=deepcopy(self)
    s.ismutant=True
    seq = deepcopy(self.seq)

    
    
    Individual.__init__(s,target=self.target,dist=self.dist, softtype=self.softtype, motif = self.motif)
    self.__drawondependencystructure__(dist,seq, target, i)
    if self.id != None:
      s.id = (self.id,)

    s.seq=seq
    return s


  # Function for creating cross-over with other individual. If
  # xoverpoint is a pair of integers, the region from and including
  # the first index and to but not including the second index is taken
  # from other, with the remaining material from self. Otherwise,
  # xoverpoint is assumed to be a set of sets of positions where
  # material is taken from other with the remaining taken from self.
  def crossover(self, other, xoverpoint):
    s = deepcopy(self)
    s.ismutant=True
    seq=deepcopy(self.seq)
    Individual.__init__(s,target=self.target,dist=self.dist, softtype=self.softtype, motif = self.motif)
    if len(xoverpoint) == 2 and type(xoverpoint[0]) == IntType:
      # Normal cross-over, switch to other and back again
      for i in xrange(xoverpoint[0], xoverpoint[1]):
        seq[i] = other.seq[i]
    else:
      # Cross-over is defined in terms of a set of sets of positions
      # that should be copied from other
      for c in xoverpoint:
        for i in c:
          seq[i] = other.seq[i]
    # All done creating cross-over sequence, create and return new
    # object of the same type with this new sequence
    
    if self.id != None or other.id != None:
      s.id = (self.id, other.id)
    s.seq=seq
    return s

# Function for selecting an index at random with probability
# proportional to the contents of array.
def _weightedselect(array):
    #when using both software programs, use only the results from RNA fold
    #to compute the weightedselection
    if type(array) is tuple:
      array=map(lambda x: x / 2, (map(operator.add,array[0],array[1]))) 
    s = sum(array)
    if s == 0:
      # All entries are 0
      return randint(0, len(array) - 1)
    x = uniform(0, sum(array))
    for i in xrange(len(array) - 1):
      if x < array[i]:
        return i
      x -= array[i]
    return len(array) - 1


### Initialisation; initialise sequence for an individual

# Class initialising sequence to be a random (but not uniformly
# picked) sequence compatible with dependency structure of targets.
class RandomInitialise:
  # Function for (re-)initialising the sequence of an individual from
  # scratch, compatible with target. New bases and base pairs are
  # drawn from the distribution specified by dist.
  def initialise(self, target = None, dist = None, softtype = None):
    if target == None:
      target = self.target
    if dist == None:
      dist = self.dist
    if softtype == None:
      softtype = self.softtype
    # Start with an undefined sequence, then draw nucleotides for each
    # dependency structure in turn.
    if self.seq != None:
      self.seq = len(self) * [None]
    else:
      self.seq = target.size() * [None]
    visited = len(self) * [None]
    order = range(len(self))
    shuffle(order)
    for i in order:
      if visited[i] is None:
        self.__drawondependencystructure__(dist, self.seq, target, i, visited)
    # Reset values so we don't use old ones
    tmp = self.id
    self.__init__(self.seq, target=target,dist=dist, softtype=softtype, motif = self.motif)
    self.id = tmp
    from subprocess import Popen, PIPE
    import os
    os.environ['PATH']=os.environ['PATH']+':/usr/local/bin'
    p = Popen("RNAeval", stdin = PIPE)
    print >> p.stdin, int2str(self.seq)
    p.stdin.close()

# Class initialising sequence by using the initial function from the
# fold module (assumed to design a sequence for a single target
# structure). If target consists of multiple structures, a random one
# is used.
class DesignInitialise:
  def initialise(self, target = None, dist = None, softtype = None, motif = None):
    if target == None:
      target = self.target
    if dist==None:
      dist=self.dist
    if softtype == None:
      softtype = self.softtype
    # Remember id so we can revert to this after resetting values
    tmp = self.id
    # Design sequence for random target structure
    self.__init__(str2int(design(target.structures[randint(0, len(target) - 1)])), target=target,dist=dist, softtype=softtype )
    self.id = tmp

# Class initialising sequence by reading a single line from file f
def FileInitialise(f):
  # Check to see whether f is provided as file or file name
  if type(f) == StringType:
    try:
      f = open(f)
    except IOError:
      raise IOError, "Could not open file " + f + "for reading initial sequences"
  class _FileInitialise(RandomInitialise):
    def initialise(self, target = None, dist = None, softtype = None ):
      if target == None:
          target = self.target
      if dist == None:
          dist = self.dist
      if softtype == None:
          softtype = self.softtype
      # Remember id so we can revert to this after resetting values
      tmp = self.id
      try:
        self.__init__(str2int("".join(f.readline().split())), self.getpopulation())
        # Check to see if sequence length matches target length
        if len(self) != self.getpopulation().gettarget().size():
          raise IOError
      except IOError:
        # If no valid sequence could be read, use random initialisation
        # instead
        RandomInitialise.initialise(self)
      self.id = tmp
    # Close file sequences are read from
    def __fileinitialiseclosefile__(self):
      f.close()
  return _FileInitialise

### MUTATIONAL WEIGHTS; interpretation is as a measure of how wrong a
### position is, on a scale from 0 to 1

# Class providing method for scoring sites based on the fraction of
# target structures it has an incorrect predicted structure for.
class WrongPrediction:
  def __getwrongprediction__(self, target = None):
    
    try:
      if self.__wrongprediction__ == None or self.ismutant:
        raise AttributeError
    except AttributeError:
      self.ismutant=False
      if target == None:
        target = self.target
      self.__wrongprediction__ = [0] * len(self.seq)
      # For each target, add to scores of positions with wrong position
      for t in target:
        #print "using cal_once from Wrong Pred"
        for i in posdif(t, self.getstruct(t.temperature)):
         # print t,'  vs  ',self.getstruct(t.temperature)
          
          self.__wrongprediction__[i] += 1.0 / len(target)
          #print self.__wrongprediction__
    return self.__wrongprediction__

# Class providing method for scoring sites based on the probability of
# it being wrong, averaged over all target structures.
class ProbabilityPrediction:
  def __getprobabilityprediction__(self, target = None):
    try:
      if self.__probabilityprediction__ == None or self.ismutant:
        raise AttributeError
    except AttributeError:
      self.ismutant=False
      if target == None:
        target = self.target
      self.__probabilityprediction__ = [0] * len(self.seq)
      self.__probabilitypredictionb__ = [0] * len(self.seq)
      # For each target, add probability of incorrect structure to each position
      for t in target:
        #Make the call to the software package of choice only once to save time
        #print "using cal_once from probability"
        cal_once=self.getprobability(t.temperature)
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            self.__probabilityprediction__[i] += (1.0 - cal_once[0][p[0]][p[1]]) / len(target)
            if len(cal_once) > 1: #if both software types are being used
              self.__probabilitypredictionb__[i] += (1.0 - cal_once[1][p[0]][p[1]]) / len(target)
          # Update weights for unpaired positions in t
          for i in t.unpaired():
            self.__probabilityprediction__[i] += (1.0 - cal_once[0][i][None]) / len(target)
            if len(cal_once) > 1:
              self.__probabilitypredictionb__[i] += (1.0 - cal_once[1][i][None]) / len(target)
    #return the pair of probabilities if both soft types are being used 
    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
      return self.__probabilityprediction__, self.__probabilitypredictionb__
    
    return self.__probabilityprediction__

# Class providing method for scoring sites based on the probability of
# it being wrong, returned as tuple for all the target structures.
class ProbabilityTuplePrediction:
  def __getprobabilitytupleprediction__(self, target = None):
    try:
      if self.__probabilitytupleprediction__ == None or self.ismutant:
        raise AttributeError
    except AttributeError:
      self.ismutant=False
      if target == None:
        target = self.target
      self.__probabilitytupleprediction__ = [[] for i in xrange(len(self.seq))]
      self.__probabilitytuplepredictionb__ = [[] for i in xrange(len(self.seq))]
      # For each target, add probability of incorrect structure to each position
      for t in target:
        cal_once=self.getprobability(t.temperature)
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:            
            self.__probabilitytupleprediction__[i].append(1.0 - cal_once[0][p[0]][p[1]])
            if len(cal_once) > 1:
              self.__probabilitytuplepredictionb__[i].append(1.0 - cal_once[1][p[0]][p[1]])
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          self.__probabilitytupleprediction__[i].append(1.0 - cal_once[0][i][None])
          if len(cal_once) > 1:
              self.__probabilitytuplepredictionb__[i].append(1.0 - cal_once[1][i][None])

    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
      return self.__probabilitytupleprediction__,self.__probabilitytuplepredictionb__
    #print "leaving probability tuple prediction"        
    return self.__probabilitytupleprediction__

# Class providing method for scoring sites based on 1 minus product of
# probabilities of having the correct structure, taken over all target
# structures.
class ProductProbabilityPrediction:
  def __getproductprobabilityprediction__(self, target = None):
    try:
      if self.__productprobabilityprediction__ == None or self.ismutant:
        raise AttributeError
      #Second else statement because mutate method makes a deep copy.  We 
      #want it to recalculate the prob prediction in either case
    except AttributeError:
      self.ismutant=False
      if target == None:
        target = self.target
      self.__productprobabilityprediction__ = [1] * len(self.seq)
      self.__productprobabilitypredictionb__ = [1] * len(self.seq)
      # For each target, multiply by probability of correct structure
      # at each position.
      for t in target:
        cal_once = self.getprobability(t.temperature)
        
        # Update values for base paired positions in t
        for p in t.basepairs():
          for i in p:
            self.__productprobabilityprediction__[i] *= cal_once[0][p[0]][p[1]]
            if len(cal_once) > 1:
              self.__productprobabilitypredictionb__[i] *= cal_once[1][p[0]][p[1]]
        # Update values for unpaired positions in t
        for i in t.unpaired():
          self.__productprobabilityprediction__[i] *= cal_once[0][i][None]
          if len(cal_once) > 1:
              self.__productprobabilitypredictionb__[i] *= cal_once[1][i][None]
              
      self.__productprobabilityprediction__ = map(lambda x: 1 - x, self.__productprobabilityprediction__)
      self.__productprobabilitypredictionb__ = map(lambda x: 1 - x, self.__productprobabilitypredictionb__)

    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
      return self.__productprobabilityprediction__, self.__productprobabilitypredictionb__
    
    return self.__productprobabilityprediction__

# Class providing method for scoring sites based on the negative log
# probability of it being correct, across the target
# structures. Values are truncated at 100 and normalised to the
# interval from 0 to 1
class LogProbabilityPrediction:
  def __getlogprobabilityprediction__(self, target = None):
    try:
      if self.__logprobabilityprediction__ == None or self.ismutant:
        raise AttributeError
    except AttributeError:
      self.ismutant=False
      if target == None:
        target = self.target
      self.__logprobabilityprediction__ = [0] * len(self.seq)
      self.__logprobabilitypredictionb__ = [0] * len(self.seq)
      # For each target, add probability of incorrect structure to each position
      for t in target:
        cal_once = self.getprobability(t.temperature)
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            try:
              self.__logprobabilityprediction__[i] += min(100, -log(cal_once[0][p[0]][p[1]])) * .01 / len(target)
              if len(cal_once) > 1:
                self.__logprobabilitypredictionb__[i] += min(100, -log(cal_once[1][p[0]][p[1]])) * .01 / len(target)
            except ValueError:
              self.__logprobabilityprediction__[i] += 1.0 / len(target)
              if len(cal_once) > 1:
                self.__logprobabilitypredictionb__[i] += 1.0 / len(target)
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          try:
            self.__logprobabilityprediction__[i] += min(100, -log(cal_once[0][i][None])) * .01 / len(target)
            if len(cal_once) > 1:
                self.__logprobabilitypredictionb__[i] += min(100, -log(cal_once[1][i][None])) * .01 / len(target)
          except ValueError:
            self.__logprobabilityprediction__[i] += 1.0 / len(target)
            if len(cal_once) > 1:
                self.__logprobabilitypredictionb__[i] += 1.0 / len(target)

    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
      return self.__logprobabilityprediction__,self.__logprobabilitypredictionb__
    return self.__logprobabilityprediction__

# Class providing method for scoring sites based on the probability of
# it being wrong not exceeding a threshold, averaged over the target
# structures.                                                       
class ThresholdProbabilityPrediction:
  def __getthresholdprobabilityprediction__(self, threshold, target = None):
    try:
      try:
        if threshold not in self.__thresholdprobabilityprediction__ or self.ismutant:
          # Scores have not previously been computed for this threshold
          self.__thresholdprobabilityprediction__[threshold] = [0] * len(self.seq)
          self.__thresholdprobabilitypredictionb__[threshold] = [0] * len(self.seq)
          raise KeyError
      except AttributeError:
        # No threshold scores available
        self.__thresholdprobabilityprediction__ = {threshold: [0] * len(self.seq)}
        self.__thresholdprobabilitypredictionb__[threshold] = [0] * len(self.seq)
        raise KeyError
    except KeyError:
      self.ismutant=False
      # Compute scores
      if target == None:
        target = self.target
      # For each target, add probability of incorrect structure to each position
      for t in target:
        cal_once = self.getprobability(t.temperature)
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            self.__thresholdprobabilityprediction__[threshold][i] += float(cal_once[0][p[0]][p[1]] < threshold) / len(target)
            if len(cal_once) > 1:
              self.__thresholdprobabilitypredictionb__[threshold][i] += float(cal_once[1][p[0]][p[1]] < threshold) / len(target)
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          self.__thresholdprobabilityprediction__[threshold][i] += float(cal_once[0][i][None] < threshold) / len(target)
          if len(cal_once) > 1:
            self.__thresholdprobabilitypredictionb__[threshold][i] += float(cal_once[1][i][None] < threshold) / len(target)
    #print "leaving threshould probability prediction"

    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
      return self.__thresholdprobabilityprediction__[threshold], self.__thresholdprobabilitypredictionb__[threshold]
    return self.__thresholdprobabilityprediction__[threshold]

# Class for finding probabilities of target structure and most
# probable alternative structure for all positions.
class RelativeProbabilityPrediction:
  def __getrelativeprobabilityprediction__(self, target = None):
    try:
      if self.__relativeprobabilityprediction__ == None is self.ismutant:
        raise AttributeError
    except AttributeError:
      self.ismutant=False
      # Compute positional scores
      if target == None:
        target = self.target
      self.__relativeprobabilityprediction__ = [[] for i in xrange(len(self.seq))]
      self.__relativeprobabilitypredictionb__ = [[] for i in xrange(len(self.seq))]
      for t in target:
        # Get probabilities
        B = self.getprobability(t.temperature)
        # Find probabilities for base paired positions in t
        for p in t.basepairs():
          for i in xrange(2):
            try:
              m = max(map(lambda x: x[1], filter(lambda x: x[0] != p[1 - i], B[0][p[i]].items())))
              if len(B) > 1:
                m2 = max(map(lambda x: x[1], filter(lambda x: x[0] != p[1 - i], B[1][p[i]].items())))
            except ValueError:
              # No other probabilities for this position
              m = 0
              m2 = 0
            self.__relativeprobabilityprediction__[p[i]].append((B[0][p[0]][p[1]], m))
            self.__relativeprobabilitypredictionb__[p[i]].append((B[1][p[0]][p[1]], m2))
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          try:
            m = max(map(lambda x: x[1], filter(lambda x: x[0] != None, B[0][i].items())))
            if len(B) >1:
              m2 = max(map(lambda x: x[1], filter(lambda x: x[0] != None, B[1][i].items())))
          except ValueError:
            # No other probabilities fpr this position
            m = 0
            m2 = 0
          self.__relativeprobabilityprediction__[i].append((B[0][i][None], m))
          if len(B) > 1:
            self.__relativeprobabilitypredictionb__[i].append((B[1][i][None], m2))
    #print "leaving realtive proability predction"
    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
      return self.__relativeprobabilityprediction__, self.__relativeprobabilitypredictionb__
    
    return self.__relativeprobabilityprediction__

# Class providing method for scoring sites based on difference between
# probability of most probable wrong configuration and probability of
# correct base pair minus one third of cube of this difference
# (integral of 1 - x^2 from -1 to difference, resulting in largest
# effect of changes around 0), rescaled to range from 0 to 1.
class TransformedRelativeProbabilityPrediction(RelativeProbabilityPrediction):
  def __gettransformedrelativeprobabilityprediction__(self, target = None):
    try:
      if self.__transformedrelativeprobabilityprediction__ == None or self.ismutant:
        raise AttributeError
    except AttributeError:
      self.ismutant=False
      if target == None:
        target = self.target
      # For each target structure, add score based on probability of
      # correct structure vs maximum probability of alternative
      # structure for each position.
      def _score(d):
        return .5 + .75 * d - .25 * pow(d, 3)
      # Compute positional scores
      #look at this again, does not work with BOTH option!!
      comp_once =  self.__getrelativeprobabilityprediction__(target)
      self.__transformedrelativeprobabilityprediction__ = map(sum, map(lambda y: map(lambda x: _score(x[1] - x[0]) / float(len(target)), y), comp_once[0]))
    if target==None:
        target=self.target
    if len(self.prob[target[0].temperature])>1:
        self.__transformedrelativeprobabilityprediction2__ = map(sum, map(lambda y: map(lambda x: _score(x[1] - x[0]) / float(len(target)), y), comp_once[1]))
        return self.__transformedrelativeprobabilityprediction__,self.__transformedrelativeprobabilityprediction2__
      
    return self.__transformedrelativeprobabilityprediction__


### CHOOSE POSITION; classes for selecting position for mutation

# Class choosing a position uniformly at random
class RandomPosition:
  # Return random position
  def getposition(self, *args):
    return randint(0, len(self.seq) - 1)
  def pweight(self, i = None, target = None):
    if i == None:
      return 0.5 * len(self.seq)
    else:
      return 0.5
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, 0.5)

# Class choosing a position based on how many target structures it is wrong in
# go use the array casting trick and fix all of the POSITION STUFF
class WrongPosition(WrongPrediction):
    
  def getposition(self, target = None):
    return _weightedselect(self.__getwrongprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      if type(self.__getwrongprediction__(target)) is tuple:
        hold=self.__getwrongprediction__(target)
        return sum(map(operator.add,hold[0],hold[1]))
      return sum(self.__getwrongprediction__(target))
    else:
      if type(self.__getwrongprediction__(target)) is tuple:
        hold=self.__getwrongprediction__(target)
        return sum(map(operator.add,hold[0],hold[1]))
      return self.__getwrongprediction__(target)[0][i]
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, WrongPosition.pweight(self, i, target))

# Class choosing a position based on probability of it not having the
# target structure(s)
class ProbabilityPosition(ProbabilityPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__getprobabilityprediction__(target))
  
  def pweight(self, i = None, target = None):
    if i == None:
      if type(self.__getprobabilityprediction__(target)) is tuple:
        hold=self.__getprobabilityprediction__(target)
        return sum(map(lambda x: x / 2, (map(operator.add,hold[0],hold[1])))) 
      return sum(self.__getprobabilityprediction__(target))
    
    else:
      if type(self.__getprobabilityprediction__(target)) is tuple:
        hold=self.__getprobabilityprediction__(target)
        return (map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))))[i]
  
      return self.__getprobabilityprediction__(target)[0][i]
    
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, ProbabilityPosition.pweight(self, i, target))

# Class choosing a position based on maximum probabilitiy of it
# not having the target structures
class MinProbabilityPosition(ProbabilityTuplePrediction):
  def getposition(self, target = None):
    return _weightedselect(map(max, self.__getprobabilitytupleprediction__(target)))
  def pweight(self, i = None, target = None):
    if i == None:
      if type(self.__getprobabilitytupleprediction_(target)) is tuple:
        hold=self.__getprobabilitytupleprediction__(target)
        return sum(map(max, map(lambda x: x / 2, (map(operator.add,hold[0],hold[1])))))
      return sum(map(max, self.__getprobabilitytupleprediction__(target)))
    else:
      if type(self.__getprobabilitytupleprediction__(target)) is tuple:
        hold=self.__getprobabilitytupleprediction__(target)
        return (max(lambda x: x / 2, (map(operator.add,hold[0],hold[1])))[i])
      return max(self.__getprobabilitytupleprediction__(target)[0][i])
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, MinProbabilityPosition.pweight(self, i, target))

# Class choosing a position based on negative log probability of it
# having the correct target structure
class LogProbabilityPosition(LogProbabilityPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__getlogprobabilityprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      if type(self.__getlogprobabilityprediction__(target)) is tuple:
        hold = self.__getlogprobabilityprediction__(target)
        return sum(map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))))
      return sum(self.__getlogprobabilityprediction__(target))
    else:
      if type(self.__getlogprobabilityprediction__(target)) is tuple:
        hold = self.__getlogprobabilityprediction__(target)
        return map(lambda x: x / 2, (map(operator.add,hold[0],hold[1])))[i]
      return self.__getlogprobabilityprediction__(target)[0][i]
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, LogProbabilityPosition.pweight(self, i, target))

# Class choosing a position where the probability of being wrong
# exceeds a given threshold. With multiple target structures, the
# probability of choosing a position is proportional to the fraction
# of target structures for which the porobability of being wrong
# exceeds the threshold.
def ThresholdProbabilityPosition(threshold):
  class _ThresholdProbabilityPosition(ThresholdProbabilityPrediction):
    def getposition(self, target = None):
      return _weightedselect(self.__getthresholdprobabilityprediction__(threshold, target))
    def pweight(self, i = None, target = None):
      if i == None:
        if type(self.__getthresholdprobabilityprediction__(threshold, target)) is tuple:
          hold = self.__getthresholdprobabilityprediction__(threshold, target)
          return sum(map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))))
        return sum(self.__getthresholdprobabilityprediction__(threshold, target))
      else:
        if type(self.__getthresholdprobabilityprediction__(threshold, target)) is tuple:
          hold = self.__getthresholdprobabilityprediction__(threshold, target)
          return (map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))))[i]
        return self.__getthresholdprobabilityprediction__(threshold, target)[0][i]
    def pweights(self, target = None):
      for i in xrange(len(self.seq)):
        yield (i, _ThresholdProbabilityPosition.pweight(self, i, target))
  return _ThresholdProbabilityPosition

# Class choosing position based on score calculated from probability
# of correct structure vs maximum probability of an incorrect
# structure
class RelativeProbabilityPosition(TransformedRelativeProbabilityPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__gettransformedrelativeprobabilityprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      if type(self.__gettransformedrelativeprobabilityprediction__(target)) is tuple:
          hold = self.__gettransformedrelativeprobabilityprediction__(target)
          return sum(map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))))
      return sum(self.__gettransformedrelativeprobabilityprediction__(target))
    else:
      if type(self.__gettransformedrelativeprobabilityprediction__(target)) is tuple:
          hold = self.__gettransformedrelativeprobabilityprediction__(target)
          return (map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))))[i]
      return self.__gettransformedrelativeprobabilityprediction__(target)[0][i]
    
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, RelativeProbabilityPosition.pweight(self, i, target))

# Class using a weighted mixture of other schemes to choose a position
# for mutation, where schemes should be a list of pairs of a weight
# and a class for choosing positions
def CombinedPosition(schemes):
  if schemes == []:
    raise ValueError, "At least one scheme should be specified for CombinedPosition"
  # Split weights and schemes to allow use of _weightedselect for
  # randomly selecting a scheme.
  __combinedpositionweights__ = map(lambda x: x[0], schemes)
  __combinedpositionclasses__ = map(lambda x: x[1], schemes)
  # Set up base class for combined scheme
  class _CombinedPosition:
    # As positional weights are normalised by sum of positional
    # weights for each scheme, only compute them once.
    __pweight__ = None
    def getposition(self, target = None):
      # Randomly choose a scheme to apply
      c = __combinedpositionclasses__[_weightedselect(__combinedpositionweights__)]
      # Choose a position using this scheme
      return c.getposition(self, target)
    def pweight(self, i = None, target = None):
      if self.__pweight__ == None:
        # Compute positional weights
        self.__pweight__ = len(self.seq) * [0]
        w = sum(__combinedpositionweights__)
        for j in xrange(len(__combinedpositionclasses__)):
          c = __combinedpositionclasses__[j]
          # Normalise contribution from each class by sum of
          # positional weights for this class and sum of weights for
          # classes.
          z = float(c.pweight(self, None, target) * w)
          if z != 0:
            z = __combinedpositionweights__[j] / z
            for j in xrange(len(self.seq)):
              self.__pweight__[j] += c.pweight(self, j, target) * z
      if max(self.__pweight__) <= 0:
        # No position with positive weight, reset all to weight 1
        self.__pweight__ = len(self.seq) * [1]
      # We now have weights
      if i == None:
        return sum(self.__pweight__)
      else:
        return self.__pweight__[i]
    def pweights(self, target = None):
      for i in xrange(len(self.seq)):
        yield (i, self.pweight(i, target))
  # Add the schemes specified as superclasses, so we can invoke their methods
  for s in schemes:
    class _CombinedPosition(_CombinedPosition, s[1]):
      pass
  # Finished building the combined class
  return _CombinedPosition

### CUT WEIGHTS; classes for selecting cuts for cross overs

# Generic class for selecting cross over cuts. This should be
# specialised to define __weightf.
class CutWeight:
  # Function for computing weight of cuts
  def __computecutweight(self, target = None):
    if target == None:
      target = self.target
    self.__cweight = {}
    # Compute prefix array of positional weights for efficient
    # computation of cut weights - note that prefix[i] is sum up to
    # index i - 1.
    prefix = [0] * (len(self.seq) + 1)
    for i in xrange(len(self.seq)):
      prefix[i + 1] = prefix[i] + self._weightf(target)[i]
    # For every pair of compatible positions, compute their cross-over score
    for c in target.getcuts():
      # A cut is represented by a set of positions that are pairwise
      # compatible as cuts.
      for i in c:
        for j in c:
          # Ensure only distinct pairs of compatible positions are
          # considered, and only considered once.
          if i < j:
            # Weight of inside structure
            self.__cweight[(i, j)] = (prefix[j] - prefix[i]) / float(len(c) - 1)
            # Weight of outside structure
            self.__cweight[(j, i)] = (prefix[-1] + prefix[i] - prefix[j]) / float(len(c) - 1)
  # Iterator over all cut weights
  def cweights(self, target = None):
    try:
      if self.__cweight == None:
        raise AttributeError
    except AttributeError:
      self.__computecutweight(target)
    for c in self.__cweight:
      yield (c, self.__cweight[c])
  # Return weight of specific cut, or sum of all weights
  def cweight(self, i = None,all=False, target = None):
    try:
      if self.__cweight == None:
        raise AttributeError
    except AttributeError:
      self.__computecutweight(target)
    if i != None:
      return self.__cweight[i]
    else:
      if all:
          return (self.__cweight)
      else:
          return sum(self.__cweight.values())

# Class weighting cuts equally
class RandomCut(CutWeight):
  # Function for computing position weights
  def _weightf(self, target):
    class __Weightf:
      n = len(self.seq)
      def __getitem__(self, idx):
        if type(idx) == IntType and idx >= -self.n and idx < self.n:
          return 0.5
        else:
          raise IndexError, "Cannot use " + repr(idx) + " as index"
    return __Weightf()

# Class weighting cuts by number of positions that are right in
# predicted structure(s).
class CorrectCut(CutWeight, WrongPrediction):
  def _weightf(self, target):
    # A prediction is right if it is not wrong
    if type(self.__getwrongprediction__(target)) is tuple:
      hold=self.__getwrongprediction__(target)
      ave=map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))) 
      return map(lambda x: 1 - x, ave)
    return map(lambda x: 1 - x, self.__getwrongprediction__(target))

# Class weighting cuts by probabilities of having right structure(s)
class ProbabilityCut(CutWeight, ProbabilityPrediction):
  def _weightf(self, target):
    if type(self.__getprobabilityprediction__(target)) is tuple:
      hold=self.__getprobabilityprediction__(target)
      ave=map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))) 
      return map(lambda x: 1 - x, ave)

      # A prediction is right if it is not wrong
    return map(lambda x: 1 - x, self.__getprobabilityprediction__(target))

# Class weighting cuts by minimum probability of having right
# structure, taken over all targets
class MinProbabilityCut(CutWeight, ProbabilityTuplePrediction):
  def _weightf(self, target):
    if type(self.__getprobabilitytupleprediction__(target)) is tuple:
      hold=self.__getprobabilitytupleprediction__(target)
      ave=map(lambda x: x / 2, (map(operator.add,hold[0],hold[1]))) 
      return map(lambda x: 1 - x, ave)
    
    return map(lambda x: 1 - max(x), self.__getprobabilitytupleprediction__(target))

# Class weighting cuts by sum of products of having correct structure,
# where the product is taken over target structures and the sum over
# positions.
class ProductProbabilityCut(CutWeight, ProductProbabilityPrediction):
  def _weightf(self, target):
    return self.__getproductprobabilityprediction__(target)

# Class weighting cuts by score based on probability of correct
# structure vs maximum probability of an incorrect structure
class RelativeProbabilityCut(CutWeight, TransformedRelativeProbabilityPrediction):
  def _weightf(self, target):
    # A prediction is right if it is not wrong
    return map(lambda x: 1 - x, self.__gettransformedrelativeprobabilityprediction__(target))

# Class weighting cuts by number of positions where probability of
# correct structures exceeds the specified threshold.
def ThresholdProbabilityCut(threshold):
  class _ThresholdProbabilityCut(CutWeight, ThresholdProbabilityPrediction):
    def _weightf(self, target):
      # A prediction is right if it is not wrong
      return map(lambda x: 1 - x, self.__getthresholdprobabilityprediction__(threshold, target))
  return _ThresholdProbabilityCut

# Class using a weighted mixture of other schemes to weight cross-over
# cuts, where schemes should be a list of pairs of a weight and a
# class for choosing positions.
def CombinedCut(schemes):
  if schemes == []:
    raise ValueError, "At least one scheme should be specified for CombinedCut"
  __combinedcutschemes__ = deepcopy(schemes)
  class _CombinedCut(CutWeight):
    # Gather positional scores from all specified schemes first time
    # it is needed.
    __combinedcutweight__ = None
    def _weightf(self, target):
      if self.__combinedcutweight__ == None:
        # Weights yet uncomputed
        for s in __combinedcutschemes__:
          tmp = s[1]._weightf(self, target)
          z = float(sum(tmp))
          if z > 0:
            # Scheme provides positive scores
            if self.__combinedcutweight__ == None:
              # First scheme with positive scores
              self.__combinedcutweight__ = map(lambda x: x * s[0] / z, tmp)
            else:
              # Already have some scores compiled
              for i in xrange(len(self.seq)):
                self.__combinedcutweight__[i] += tmp[i] * s[0] / z
        if self.__combinedcutweight__ == None:
          # No positive scores
          self.__combinedcutweight__ = len(self.seq) * [1.0]
      # Scores either pre-existing or just computed
      return self.__combinedcutweight__
  # Add classes in schemes as superclasses, so their _weightf methods
  # can be invoked.
  for s in schemes:
    class _CombinedCut(_CombinedCut, s[1]):
      pass
  return _CombinedCut

### Simple class for maintaining an array of weights and their sum, and
### for drawing a random element from it.
class Weights:
  def __init__(self, w):
    self.weights = w
    self.sum = sum(w)
  def __getitem__(self, i):
    return self.weights[i]
  def __setitem__(self, i, x):
    self.sum += x - self.weights[i]
    self.weights[i] = x
  # Make a weighted draw of an index, ignoring any prohibited ones
  def draw(self, prohibited = None):
    # Draw random real number in the range of sum of weights of
    # accessible indeces.
    if prohibited == None:
      # If no indeces are prohibited, we don't need to check
      x = uniform(0, self.sum)
      for i in xrange(len(self.weights)):
        if x < self.weights[i]:
          return i
        x -= self.weights[i]
      # Overflow
      if len(self.weights) > 0:
        j = len(self.weights) - 1
      else:
        j = None
    else:
      # Remember last allowed index, in case of overflow
      j = None
      # Draw a random item, but avoid prohibited indeces
      if type(prohibited) == IntType:
        # Only one prohibited index
        x = uniform(0, self.sum - self.weights[prohibited])
        for i in xrange(len(self.weights)):
          if i != prohibited:
            if x < self.weights[i]:
              return i
            x -= self.weights[i]
            j = i # Remember last index, in case of overflow
      else:
        # List of prohibited indeces
        s = self.sum
        for i in prohibited:
          s -= self.weights[i]
        x = uniform(0, s)
        for i in xrange(len(self.weights)):
          if i not in prohibited:
            if x < self.weights[i]:
              return i
            x -= self.weights[i]
            j = i # Remember last index, in case of overflow
    # Overflow - return highest allowed index
    return j

### Default distribution to draw unpaired bases and base pairs from
class BaseDistribution:
  def __init__(self, unpaired = None, basepair = None, paired = None):
    # Unpaired is distribution an unpaired position is drawn from
    if unpaired == None:
      self.unpaired = Weights(4 * [1])
    else:
      if isinstance(unpaired, Weights):
        self.unpaired = unpaired
      else:
        self.unpaired = Weights(unpaired)
    # Basepair is the distribution a base pair is drawn from - first
    # pair in a dependency structure is drawn from this distribtuion.
    if basepair == None:
      self.basepair = Weights(6 * [1])
    else:
      if isinstance(basepair, Weights):
        self.basepair = basepair
      else:
        self.basepair = Weights(basepair)
    # Paired is the distribution positions paired to other positions
    # that have already been fixed are drawn from = subsequent
    # positions in dependency structures are drawn from this
    # distribution.
    if paired == None:
      self.paired = Weights(4 * [1])
    else:
      if isinstance(paired, Weights):
        self.paired = paired
      else:
        self.paired = Weights(paired)
  # Draw a base for unpaired position - if current is not None, this
  # will be prohibited.
  def drawunpaired(self, current = None):
    return self.unpaired.draw(current)
  # Draw a basepair - if current is an integer, any pair with that as
  # first base is avoided, if current is a pair, that specific pair is
  # avoided.
  def drawbasepair(self, current = None):
    if type(current) == IntType:
      if current > 1:
        return self.basepair.draw([current, current + 2])
      else:
        return self.basepair.draw(current)
    elif current != None:
      return self.basepair.draw(bp2int(current))
    else:
      return self.basepair.draw()
  # Draw a paired base, given bases fixed at paired positions are
  # given by the set neighbours. If current is an integer, that
  # particular base is avoided, if current is a pair it is assumed to
  # be an integer representing the current base and the probability
  # the position is correct; if current is one of two possibilities, it
  # is reweighted by this probability.
  def drawpaired(self, neighbours = None, current = None):
    if neighbours == None or len(neighbours) == 0:
      # No fixed neighbours, draw from paired distribution
      if type(current) == IntType or current == None:
        return self.paired.draw(current)
      else:
        # For the first position changed in a dependency structure, we
        # definitely want a new nucleotide, as this is the mutating
        # site.
        return self.paired.draw(current[0])
    elif 0 in neighbours:
      # A neighbouring position has an A as fixed nucleotide
      return 3
    elif 1 in neighbours:
      # A neighbouring position has a C as fixed nucleotide
      return 2
    else:
      # Position has fixed neighbours, but none that uniquely defines
      # this one.
      if 2 in neighbours:
        # Valid nucleotides are C and U
        offset = 1
      else:
        # Valid nucleotides are A and G
        offset = 0
      if type(current) == IntType:
        # Current is an integer, make sure we don't choose this
        w = [int(current != offset + 2 * i) * self.paired[offset + 2 * i] for i in xrange(2)]
      elif current != None:
        # Current is a pair, so modify by weights
        w = [(current[1] + int(current[0] == offset + 2 * i) * (1 - 2 * current[1])) * self.paired[offset + 2 * i] for i in xrange(2)]
      else:
        # No restrictions
        w = [self.paired[offset + 2 * i] for i in xrange(2)]
      x = uniform(0, sum(w))
      if x < w[0]:
        return offset
      else:
        return offset + 2

### MUTATE DEPENDENT; class for adding information from current
### status of position to distribution from which new nucleotide is
### drawn when mutation is caused by path of dependencies to position
### with forced mutation.

# Simple class not adding information
class MutateDependent:
  def __mutatedependent__(self, i, dist, neighbours):
    return dist.drawpaired(neighbours)

# Class using positional weights to bias draw
class WeightedMutateDependent:
  def __mutatedependent__(self, i, dist, neighbours):
    return dist.drawpaired(neighbours, (self.seq[i], self.pweight(i)))

# Class retaining current nucleotide, if possible
class RetainMutateDependent:
  def __mutatedependent__(self, i, dist, neighbours):
    for n in neighbours:
      if self.seq[i] + n not in [3, 5]:
        # Cannot retain current nucleotide
        return dist.drawpaired(neighbours)
    # Can retain current nucleotide
    return self.seq[i]

### Auxiliary classes for computing distance to target
class PredictionErrors(WrongPrediction):
  def __predictionerrors__(self):
    
    if type(self.__getwrongprediction__()) is tuple:
        
      return sum(map(lambda x: x / 2, (map(operator.add,self.__getwrongprediction__()[0],self.__getwrongprediction__()[1]))) )
    return sum(self.__getwrongprediction__())

class ProbabilityErrors(ProbabilityPrediction):
  def __probabilityerrors__(self):
    if type(self.__getprobabilityprediction__()) is tuple: 
      return sum(map(lambda x: x / 2, (map(operator.add,self.__getprobabilityprediction__()[0],self.__getprobabilityprediction__()[1]))) )
    return sum(self.__getprobabilityprediction__())

class MinProbabilityErrors(ProbabilityTuplePrediction):
  def __minprobabilityerrors__(self):
    if type(self.__getprobabilitytupleprediction__()) is tuple:
      array=map(lambda x: x / 2, (map(operator.add,self.__getprobabilitytupleprediction__()[0],self.__getprobabilitytupleprediction__()[1]))) 
      return (sum(map(max, array)))
    return sum(map(max, self.__getprobabilitytupleprediction__()))

class ProductProbabilityErrors(ProductProbabilityPrediction):
  def __productprobabilityerrors__(self):
    if type(self.__getproductprobabilityprediction__()) is tuple:
      array=map(lambda x: x / 2, (map(operator.add,self.__getproductprobabilityprediction__()[0],self.__getproductprobabilityprediction__()[1]))) 
      return sum(array)
    return sum(self.__getproductprobabilityprediction__())

class RelativeProbabilityErrors(TransformedRelativeProbabilityPrediction):
  def __relativeprobabilityerrors__(self):
    if type(self.__gettransformedrelativeprobabilityprediction__()) is tuple:
      array=map(lambda x: x / 2, (map(operator.add,gettransformedrelativeprobabilityprediction__[0],self.__gettransformedrelativeprobabilityprediction__()[1]))) 
      return sum(array)
    return sum(self.__gettransformedrelativeprobabilityprediction__())

class ThresholdProbabilityErrors(ThresholdProbabilityPrediction):
  def __thresholdprobabilityerrors__(self, threshold):
    if type(self.__getthresholdprobabilityprediction__(threshold)) is tuple:
      array=map(lambda x: x / 2, (map(operator.add,self.__getthresholdprobabilityprediction__(threshold)[0],self.__getthresholdprobabilityprediction__(threshold)[1]))) 
      return sum(array)       
    return sum(self.__getthresholdprobabilityprediction__(threshold))

def EnsembleProbability(w):
  R = 1.9858775e-3
  class _EnsembleProbability:
    def __ensembleprobability__(self):
      a, b = 0, 0
      target = self.target
      for t in target:
        e = self.getstructurescore(t, t.temperature)
        z = self.getensemblescore(t.temperature)
        if t.temperature == None:
          T = 310.15
        else:
          T = t.temperature + 273.15
        a += (e - z) / (R * T)
        b += pow((e - z) / (R * T), 2)
      return a / len(target) + w * (b / len(target) - pow(a / len(target), 2))
  return _EnsembleProbability


### OBJECTIVE; classes for computing the objective of an individual
# Simple class measuring objective as the number of errors. A perfect
# solution is one with no errors.
class WrongObjective(PredictionErrors):
  def objective(self):
    return self.__predictionerrors__()
  def perfectsolution(self):
    return self.__predictionerrors__() == 0

# Class measuring objective as expected number of errors, according to
# probability distribution. A perfect solution is one where each
# position is incorrect with probability at most precision.
def ProbabilityObjective(precision = 0):
  class _ProbabilityObjective(ProbabilityErrors):
    def objective(self):
      return self.__probabilityerrors__()
    def perfectsolution(self):
      if precision == 0:
        return self.__probabilityerrors__() == 0
      elif precision >= 1:
        # No demands on perfect solution
        return True
      else:
        for i in self.__getprobabilityprediction__():
          if i > precision:
            return False
        return True
  return _ProbabilityObjective

# Class measuring objective as maximum expected incorrectness,
# according to probability distribution and maximised over all targets
# for each position. A perfect solution is one where each position has
# a maximum expected error of at most precision.
def MinProbabilityObjective(precision = 0):
  class _MinProbabilityObjective(MinProbabilityErrors):
    def objective(self):
      return self.__minprobabilityerrors__()
    def perfectsolution(self):
      if precision == 0:
        return self.__minprobabilityerrors__() == 0
      elif precision >= 1:
        # No demands on perfect solution
        return True
      else:
        return reduce(lambda x, y: x and y, map(lambda x: max(x) <= precision, self.__getprobabilitytupleprediction__()), True)
  return _MinProbabilityObjective

# Class measuring objective as sum of negative log probability of each
# position having the correct structure. A perfect solution is one
# where the log probability of correct structure does not exceed
# precision for any position.
def ProductProbabilityObjective(precision = 0):
  class _ProductProbabilityObjective(ProductProbabilityErrors):
    __productprobabilityobjectiveprecision__ = precision
    def objective(self):
      return self.__productprobabilityerrors__()
    def perfectsolution(self):
      if precision == 0:
        return self.__productprobabilityerrors__() == 0
      elif precision >= 1:
        # No demands on perfect solution
        return True
      else:
        for i in self.__getproductprobabilityprediction__():
          if i > self.__productprobabilityobjectiveprecision__:
            return False
        return True
  return _ProductProbabilityObjective

# Class measuring objective by difference between probability of
# correct structure and most probable incorrect structure. A perfect
# solution is one where the difference between these probabilities is
# at least precision at every position.
def RelativeProbabilityObjective(precision = 1):
  class _RelativeProbabilityObjective(RelativeProbabilityErrors):
    def objective(self):
      return self.__relativeprobabilityerrors__()
    def perfectsolution(self):
      if precision == 1:
        return self.__relativeprobabilityerrors__() == 0
      elif precision <= -1:
        # No demands on perfect solution
        return True
      else:
        for i in xrange(len(self.seq)):
          if sum(map(lambda x: x[0] - x[1], self.__getrelativeprobabilityprediction__()[i])) < precision * len(self.__getrelativeprobabilityprediction__()[i]):
            # Position i has smaller average difference than required
            return False
        # Did not encounter any positions with too small average difference
        return True
  return _RelativeProbabilityObjective

# Class measuring objective as fraction of positions not having a
# probability of being correct of at least threshold. A perfect
# solution is a solution where no positions fail this criteria.
def ThresholdProbabilityObjective(threshold):
  class _ThresholdProbabilityObjective(ThresholdProbabilityErrors):
    def objective(self):
      return self.__thresholdprobabilityerrors__(threshold)
    def perfectsolution(self):
      return self.__thresholdprobabilityerrors__(threshold) == 0
  return _ThresholdProbabilityObjective

# Class measuring objective by log probabilities of target structures
# in ensemble of all structures. The variance in these probabilities
# is added with weight w, and a solution is considered perfect if the
# score does not exceed threshold.
def EnsembleProbabilityObjective(w, threshold):
  class _EnsembleProbabilityObjective(EnsembleProbability(w)):
    def objective(self):
      return self.__ensembleprobability__()
    def perfectsolution(self):
      return self.__ensembleprobability__() < threshold
  return _EnsembleProbabilityObjective

# Generic class for combining multiple objectives, combine specifies
# how perfect solutions should be determined.
def CombinedObjective(schemes, combine):
  if schemes == []:
    raise ValueError, "At least one scheme should be specified for CombinedObjective"
  # Set up base class for combined scheme
  class _CombinedObjective:
    def __combinedobjectivecombine__(self, x, y):
      return combine(x, y)
    # Objective is normalised sum of objectives of schemes
    def objective(self):
      z = float(sum(map(lambda x: x[0], schemes)))
      obj = 0
      for c in schemes:
        try:
          obj += c[0] * c[1].objective(self) / z
        except ZeroDivisionError:
          # Normalising constant is 0
          if c[0] * c[1].objective(self) < 0:
            return -1
          elif c[0] * c[1].objective(self) > 0:
            return 1
      return obj
    def perfectsolution(self):
      return reduce(self.__combinedobjectivecombine__, map(lambda x: (x[0], x[1].perfectsolution(self)), schemes))
  # Add schemes to base class
  for c in schemes:
    class _CombinedObjective(_CombinedObjective, c[1]):
      pass
  return _CombinedObjective

# Class combining multiple objectives, where a solution is perfect if
# it is perfect under all schemes.
def ConjunctiveObjective(schemes):
  return CombinedObjective(schemes, lambda x, y: x[1] and y[1])

# Class combining multiple objectives, where a solution is perfect if
# it is perfect under at least one scheme.
def DisjunctiveObjective(schemes):
  return CombinedObjective(schemes, lambda x, y: x[1] or y[1])


### FITNESS; classes for computing the fit of an individual
# Simple class measuring fraction of positions that are incorrect
class WrongFitness(PredictionErrors):
  def fitness(self):
    return self.__predictionerrors__() / float(len(self.seq))

# Class measuring fit as fraction of expected number of incorrect
# positions, according to probability distribution.
class ProbabilityFitness(ProbabilityErrors):
  def fitness(self):
    return self.__probabilityerrors__() / float(len(self.seq))

# Class measuring fit as maximum expected incorrectness over all
# targets, i.e. penalising for positions with a low minimum expected
# correctness.
class MinProbabilityFitness(MinProbabilityErrors):
  def fitness(self):
    return self.__minprobabilityerrors__() / float(len(self.seq))

# Class measuring fit as sum of negative log probabilities
class ProductProbabilityFitness(ProductProbabilityErrors):
  def fitness(self):
    return self.__productprobabilityerrors__() / float(len(self.seq))

# Class measuring fit as a score based on probability of correct
# structure compared with maximum probability of incorrect structure.
class RelativeProbabilityFitness(RelativeProbabilityErrors):
  def fitness(self):
    return self.__relativeprobabilityerrors__() / float(len(self))

# Class measuring fit as fraction of positions not having a
# probability of being correct of at least threshold.
def ThresholdProbabilityFitness(threshold):
  class _ThresholdProbabilityFitness(ThresholdProbabilityErrors):
    def fitness(self):
      return self.__thresholdprobabilityerrors__(threshold) / float(len(self.seq))
  return _ThresholdProbabilityFitness

# Class measuring fit by log probabilities of target structures in
# ensemble of all structures divided by target length. The variance in
# these probabilities is added with weight w.
def EnsembleProbabilityFitness(w):
    class _EnsembleProbabilityFitness(EnsembleProbability(w)):
        def fitness(self):
            return self.__ensembleprobability__() / float(len(self.target))
            
    return _EnsembleProbabilityFitness

# Class combining multiple fitnesses
def CombinedFitness(schemes):
  if schemes == []:
    raise ValueError, "At least one scheme should be specified for CombinedFitness"
  # Set up base class for combined scheme
  class _CombinedFitness:
    # Fitness is normalised sum of fitnesses of schemes
    def fitness(self):
      z = float(sum(map(lambda x: x[0], schemes)))
      fitness = 0
      fitnessb = 0
      for c in schemes:
        fits = c[1].fitness(self)
        #if Not using both softwares at the same time do the standard thing
        if type(fits) is tuple:
          fits_single = fits[0]
        else:
          fits_single = fits
    
        try:
          fitness += c[0] * fits_single / z
        except ZeroDivisionError:
          print "A BAD THING HAPPEND!!!"
            # Normalising constant is 0
          if c[0] * fits_single < 0:
            return -1
          elif c[0] * fits_single > 0:
            return 1
      
        if type(fits) is tuple:
          fits_both = fits[1]
          try:
            fitnessb += c[0] * fits_both / z
          except ZeroDivisionError:
            print "A BAD THING HAPPEND!!!"
            # Normalising constant is 0
            if c[0] * fits_both < 0:
              return -1
            elif c[0] * fits_both > 0:
              return 1

      if(self.getsofttype() == 'Both'):
        return fitness, fitnessb 
      return fitness #a list of all of the fitnesses

    # Use Pareto sorting on list of fitness scores
    def paretofitness(self):
      fitness=[ ]
      for c in schemes:
          fitness.append(c[0] * c[1].fitness(self))
      return fitness


    # Dominated fitness score
    def dominated(self, other):
        dominated = False
        flag = [False, False]
        for i in range(0, len(self.getparetofitness())):
            if self.getparetofitness()[i] < other.getparetofitness()[i]:
                flag[0] = True
            elif self.getparetofitness()[i] > other.getparetofitness()[i]:
                flag[1] = True
        if flag[0] == True and flag[1] == False:
            dominated = True 
        return dominated

      
  # Add schemes as superclasses to base class
  for c in schemes:
    class _CombinedFitness(_CombinedFitness, c[1]):
      pass
  return _CombinedFitness

# Functions for standardising target structure without or without pseudoknots
def transform(self):
    t = ""
    for i in self:
      if i in "[{(<":
        t += '('
      elif i in "]})>":
        t += ')'
      elif i in ".-_:":
        t += '.'
      elif not i.isspace():
        raise ValueError, "Unexpected character " + i + "in target specification"
    return t

def pktransform(self):
    t = ""
    for i in self:
      if i not in "(){}[]<>.-_:":
        raise ValueError, "Unexpected character " + i + "in target specification"
      if i in ".-_:":
        t += '.'
      else:
        t += str(i)
    return t
