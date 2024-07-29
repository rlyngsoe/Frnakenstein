# The different classes used to construct a population in the genetic
# algorithm.
from random import uniform, randint, shuffle, choice, sample
from math import sqrt
from types import IntType
from copy import deepcopy
from bisect import bisect
from structure import int2bp, bp2int
from individual import BaseDistribution
from sys import stderr

### POPULATION; general class containing the core of a GA population
class Population:
  # Initialise population with n individuals. The target structure(s)
  # for the population is given by target, the class structure for
  # individuals by indclass, and basedist is used to draw nucleotides
  # when initialising and mutating sequences.
  def __init__(self, target, indclass, basedist = BaseDistribution(), n = 0):
    self.members = set([]) # Current set of individuals
    self.new = set([]) # Set of individuals created waiting to be added
    self.target = target # Target structure(s) for population
    self.__individualclass__ = indclass # Class for creating individuals
    self.__basedistribution__ = basedist # Distribution to draw bases from
    for i in xrange(n):
      self.members.add(indclass(population = self))
  # Return target structure(s) for this population
  def gettarget(self):
    return self.target
  # Return distribution for drawing new nucleotides
  def getbasedistribution(self):
    return self.__basedistribution__
  # Add n more individuals to population
  def addrandom(self, n, idfunc = None):
    if idfunc == None:
      idfunc = (lambda x: None)
    for i in xrange(n):
      new = self.__individualclass__(population = self, id = idfunc(i))
      new.initialise()
      self.new.add(new)
  # Add waiting new individuals to current population
  def addnew(self):
    self.members.update(self.new)
    self.new = set([])
  # Size of population
  def __len__(self):
    return len(self.members)
  # Iterator for population
  def __iter__(self):
    for i in self.members:
      yield i

### MUTATE; classes for selecting a set of mutation events creating
### new individuals.

### DEBUG
from structure import int2str

# Simple class that uses each sequence in current population equally many times
class Mutate:
  # Add n mutated sequences to population (one for each individual in
  # current population if n is undefined), using each current
  # individual as starting point the same number of times.
  def mutate(self, n = None):
    if n == None:
      n = len(self)
    for i in self.members:
      for j in xrange(n / len(self)):
        self.new.add(i.mutate(i.getposition()))
    if n % len(self) != 0:
      # Some current individuals have to be used as staring point one
      # more time - choose these at random.
      c = sample(range(len(self)), n % len(self))
      c.sort()
      j = 0
      k = 0
      for i in self.members:
        if c[j] == k:
          self.new.add(i.mutate(i.getposition()))
          j += 1
          if j == len(c):
            break
        k += 1

# Class choosing starting points independently uniformly at random
class RandomMutate:
  ### Add n mutated sequences (the same number as in the current
  ### population if n is undefined)
  def mutate(self, n = None):
    if n == None:
      n = len(self)
    m = list(self.members)
    for j in xrange(n):
      i = choice(m)
      self.new.add(i.mutate(i.getposition()))

# Class choosing starting points according to fitness
class FitnessWeightedMutate:
  ### Add n mutated sequences (the same number as in the current
  ### population if n is undefined)
  def mutate(self, n = None):
    if n == None:
      n = len(self)
    m = list(self.members)
    # Generate prefix sum array of fitnesses
    w = [1 - m[0].fitness()] + (len(self) - 1) * [0]
    for i in xrange(1, len(self)):
      w[i] = w[i - 1] + 1 - m[i].fitness()
    # Now create mutants
    for j in xrange(n):
      i = m[bisect(w, uniform(0, w[-1]), hi = len(self) - 1)]
      self.new.add(i.mutate(i.getposition()))

### RECOMBINE; classes for selecting a set of crossover events creating
### new individuals.

# Simple class just choosing pairs at random
class SelectPair:
  def __selectpairs__(self, n):
    m = list(self.members)
    return [sample(m, 2) for i in xrange(n)]

# Class choosing n pairs based on fitness - it is assumed that fitness
# is a number between 0 and 1, with 0 being a perfectly fit
# individual.
class WeightedSelectPair:
  def __selectpairs__(self, n):
    pairs = []
    m = list(self.members)
    # Generate prefix sum array of fitnesses
    w = [1 - m[0].fitness()] + (len(self) - 1) * [0]
    for i in xrange(1, len(self)):
      w[i] = w[i - 1] + 1 - m[i].fitness()
    for i in xrange(n):
      # Choose first individual
      a = bisect(w, uniform(0, w[-1]), hi = len(self) - 1)
      # Choose second individual
      x = uniform(0, w[-1] - m[a].fitness())
      if x >= w[a - 1]:
        # Choose beyond a
        x += m[a].fitness()
      b = bisect(w, x, hi = len(self) - 1)
      pairs.append((m[a], m[b]))
    return pairs

# Class choosing pairs based on combined fitness
class CombinedSelectPair:
  def __selectpairs__(self, n):
    # Compute weight for combining all pairs
    m = list(self.members)
    w = len(m) * (len(m) - 1) / 2 * [0.0]
    for a in xrange(1, len(m)):
      for b in xrange(a):
        # For each pair, the weight is the combined fitness summed
        # over all crossover points.
        idx = a * (a - 1) / 2 + b
        for c in self.gettarget().getcuts():
          for i in c:
            for j in c:
              if i != j:
                w[idx] += pow(float(m[a].cweight((i, j)) + m[b].cweight((j, i))), 2) / len(c)
    # Choose n pairs according to the computed weights
    pairs = []
    # Change w to prefix sum array
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    for i in xrange(n):
      j = bisect(w, uniform(0, w[-1]), hi = len(w) - 1)
      a = int(sqrt(2 * j)) + 1 # Either a or a + 1
      a = int(sqrt(2 * j + a)) # 2x is between a^2 - a and a^2 + a - 2
      b = j - a * (a - 1) / 2  # Once we know a, b is easy
      pairs.append((m[b], m[a]))
    return pairs

# Class for choosing cross over points uniformly at random
class SelectCut:
  # Select a cut for recombining a and b such that positions that can
  # be part of a cut are chosen uniformly at random.
  def __selectcut__(self, a, b):
    cuts = filter(lambda x: len(x) > 1, self.gettarget().getcuts())
    w = map(lambda x: len(x), cuts)
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    # Select dependency class at random weighted by size
    i = bisect(w, randint(0, w[-1] - 1), hi = len(w) - 1)
    # Choose random pair of positions in dependency class
    return sample(cuts[i], 2)

# Class for choosing cross over points weighted by current fitness
class WeightedSelectCut:
  def __selectcut__(self, a, b):
    # Compute prefix sums of fitnesses
    w = []
    cuts = []
    for c in self.gettarget().getcuts():
      for i in c:
        for j in c:
          if i != j:
            w.append(a.cweight((j, i)) + b.cweight((i, j)))
            cuts.append((i, j))
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    # Select a cut proportional to its weight
    i = bisect(w, uniform(0, w[-1]), hi = len(w) - 1)
    return cuts[i]

class Recombine:
  # Add n sequences obtained by crossover (the same number as in the
  # current population if n is undefined).
  def recombine(self, n = None):
    # Choose pairs to create recombinants from
    if n == None:
      n = len(self)
    pairs = self.__selectpairs__(n)
    # Choose recombination point for each pair, and create recombinant
    for p in pairs:
      c = self.__selectcut__(p[0], p[1])
      self.new.add(p[0].crossover(p[1], c))

### Reduce: Classes for reducing the population size

# Method for eliminating duplicates in m, as long as we still have
# at least n remaining elements.
def _removeduplicates(m, n):
  s = set([])
  k = len(m) - 1
  i = 0
  while i <= k:
    if tuple(m[i].seq) in s:
      # Seen this sequence before, eliminate from m
      m[i], m[k] = m[k], m[i]
      k -= 1
      if k < n:
        # Can't afford to eliminate any more
        break
    else:
      # New sequece, remember and proceed with next one
      s.add(tuple(m[i].seq))
      i += 1
  del m[k + 1:]
  return m

# Simple class just retaining the fittest individuals
class Reduce:
  def reduce(self, n):
    if len(self) > n:
      m = _removeduplicates(list(self.members), n)
      shuffle(m)
      m.sort(lambda x, y: int(x.getfitness() >= y.getfitness()) - int(x.getfitness() < y.getfitness()))
      self.members = set(map(lambda x: x, m[:n]))

# Class selecting based on a weighted combination of fitness and
# diversity. This class is returned by a function that takes the
# weight as parameter, where the score of an individual in a
# particular stage of the culling is computed as the fitness plus the
# weight times the average fraction of positions where a sequence
# differs from the already selected sequences.
def DiversityReduce(weight):
  class __DiversityReduce:
    __diversityweight__ = weight
    def reduce(self, n):
      if len(self) > n:
        def sim(a, b):
          c = 0
          for i in xrange(len(a.seq)):
            if a.seq[i] == b.seq[i]:
              c += 1
          return float(c) / len(a.seq)
        m = _removeduplicates(list(self.members), n)
        shuffle(m)
        w = [i.getfitness() for i in m]
        # Find most fit individual
        keep = 0
        for i in xrange(1, len(m)):
          if w[i] < w[keep]:
            keep = i
        ### DEBUG
        if w[keep] > min(w):
          print >> stderr, "WARNING: Didn't detect minimum"
        if m[keep].getfitness() > min(w):
          print >> stderr, "WARNING: Didn't store minimum"

        m[keep], w[keep], keep = m[-1], w[-1], [m[keep]]
        ### DEBUG
        if keep[0].getfitness() > min(w):
          print >> stderr, "WARNING: Didn't remember minimum"

        # Select remaining n - 1 individuals
        d = [sim(i, keep[0]) for i in m[:-1]]
        for i in xrange(1, n):
          # Find best score
          best = w[0] + d[0] * self.__diversityweight__ / float(i)
          besti = 0
          for j in xrange(1, len(m) - i):
            s = w[j] + d[j] * self.__diversityweight__ / float(i)
            if s < best:
              best = s
              besti = j
          # Move best scoring individual to keep and remove it from remaining
          keep.append(m[besti])
          if i != n - 1:
            m[besti] = m[-1 - i]
            w[besti] = w[-1 - i]
            d[besti] = d[-1 - i]
            # Update diversity scores
            for j in xrange(len(m) - 1 - i):
              d[j] += sim(m[j], keep[-1])
        # All done selecting the keepers
        self.members = set(keep)
  return __DiversityReduce
