# The different classes used to construct individuals in the genetic
# algorithm population.
from random import uniform, randint, shuffle
from math import log
from types import IntType
from copy import deepcopy
from structure import int2bp, int2str, bp2int, posdist, posdif, Structure
from fold import fold, boltzmann

### POSITIONAL WEIGHTS; interpretation is as a measure of how wrong a
### position is, on a scale from 0 to 1

# Function for selecting an index at random with probability
# proportional to the contents of array.
def _weightedselect(array):
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

# Class providing method for scoring sites based on the fraction of
# target structures it has an incorrect predicted structure for.
class WrongPrediction:
  def __getwrongprediction__(self, target = None):
    try:
      if self.__wrongprediction__ == None:
        raise AttributeError
    except AttributeError:
      if target == None:
        target = self.getpopulation().gettarget()
      self.__wrongprediction__ = [0] * len(self.seq)
      # For each target, add to scores of positions with wrong position
      for t in target:
        for i in posdif(t, self.getstruc(t.temperature)):
          self.__wrongprediction__[i] += 1.0 / len(target)
    return self.__wrongprediction__

# Class providing method for scoring sites based on the Boltzmann
# probability of it being wrong, averaged over all target structures.
class BoltzmannPrediction:
  def __getboltzmannprediction__(self, target = None):
    try:
      if self.__boltzmannprediction__ == None:
        raise AttributeError
    except AttributeError:
      if target == None:
        target = self.getpopulation().gettarget()
      self.__boltzmannprediction__ = [0] * len(self.seq)
      # For each target, add probability of incorrect structure to each position
      for t in target:
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            self.__boltzmannprediction__[i] += (1.0 - self.getboltzmann(t.temperature)[p[0]][p[1]]) / len(target)
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          self.__boltzmannprediction__[i] += (1.0 - self.getboltzmann(t.temperature)[i][None]) / len(target)
    return self.__boltzmannprediction__

# Class providing method for scoring sites based on the Boltzmann
# probability of it being wrong, returned as tuple for all the target
# structures.
class BoltzmannTuplePrediction:
  def __getboltzmanntupleprediction__(self, target = None):
    try:
      if self.__boltzmanntupleprediction__ == None:
        raise AttributeError
    except AttributeError:
      if target == None:
        target = self.getpopulation().gettarget()
      self.__boltzmanntupleprediction__ = [[] for i in xrange(len(self.seq))]
      # For each target, add probability of incorrect structure to each position
      for t in target:
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            self.__boltzmanntupleprediction__[i].append(1.0 - self.getboltzmann(t.temperature)[p[0]][p[1]])
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          self.__boltzmanntupleprediction__[i].append(1.0 - self.getboltzmann(t.temperature)[i][None])
    return self.__boltzmanntupleprediction__

# Class providing method for scoring sites based on the negative log
# Boltzmann probability of it being correct, across the target
# structures. Values are truncated at 100 and normalised to the
# interval from 0 to 1
class LogBoltzmannPrediction:
  def __getlogboltzmannprediction__(self, target = None):
    try:
      if self.__logboltzmannprediction__ == None:
        raise AttributeError
    except AttributeError:
      if target == None:
        target = self.getpopulation().gettarget()
      self.__logboltzmannprediction__ = [0] * len(self.seq)
      # For each target, add probability of incorrect structure to each position
      for t in target:
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            try:
              self.__logboltzmannprediction__[i] += min(100, -log(self.getboltzmann(t.temperature)[p[0]][p[1]])) * .01 / len(target)
            except ValueError:
              self.__logboltzmannprediction__[i] += 1.0 / len(target)
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          try:
            self.__logboltzmannprediction__[i] += min(100, -log(self.getboltzmann(t.temperature)[i][None])) * .01 / len(target)
          except ValueError:
            self.__logboltzmannprediction__[i] += 1.0 / len(target)
    return self.__logboltzmannprediction__

# Class providing method for scoring sites based on the Boltzmann
# probability of it being wrong not exceeding a threshold, averaged
# over the target structures.
class ThresholdBoltzmannPrediction:
  def __getthresholdboltzmannprediction__(self, threshold, target = None):
    try:
      try:
        if threshold not in self.__thresholdboltzmannprediction__:
          # Scores have not previously been computed for this threshold
          self.__thresholdboltzmannprediction__[threshold] = [0] * len(self.seq)
          raise KeyError
      except AttributeError:
        # No threshold scores available
        self.__thresholdboltzmannprediction__ = {threshold: [0] * len(self.seq)}
        raise KeyError
    except KeyError:
      # Compute scores
      if target == None:
        target = self.getpopulation().gettarget()
      # For each target, add probability of incorrect structure to each position
      for t in target:
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in p:
            self.__thresholdboltzmannprediction__[threshold][i] += float(self.getboltzmann(t.temperature)[p[0]][p[1]] < threshold) / len(target)
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          self.__thresholdboltzmannprediction__[threshold][i] += float(self.getboltzmann(t.temperature)[i][None]) / len(target)
    return self.__boltzmannprediction__[threshold]

# Class providing method for scoring sites based on difference between
# probability of most probable wrong configuration and probability of
# correct base pair minus one third of cube of this difference
# (integral of 1 - x^2 from -1 to difference, resulting in largest
# effect of changes around 0), rescaled to range from 0 to 1.
class RelativeBoltzmannPrediction:
  def __getrelativeboltzmannprediction__(self, target = None):
    try:
      if self.__relativeboltzmannprediction__ == None:
        raise AttributeError
    except AttributeError:
      # Compute positional scores
      if target == None:
        target = self.getpopulation().gettarget()
      self.__relativeboltzmannprediction__ = [0] * len(self.seq)
      # For each target, add score based on probability of correct
      # structure vs maximum probability of alternative structure for
      # each position.
      def _score(d):
        return .5 + .75 * d - .25 * pow(d, 3)
      for t in target:
        # Get Boltzmann probabilities
        B = self.getboltzmann(t.temperature)
        # Update weights for base paired positions in t
        for p in t.basepairs():
          for i in xrange(2):
            try:
              m = max(map(lambda x: x[1], filter(lambda x: x[0] != p[1 - i], B[p[i]].items())))
            except ValueError:
              # No other probabilities fpr this position
              m = 0
            self.__relativeboltzmannprediction__[p[i]] += _score(m - B[p[0]][p[1]]) / len(target)
        # Update weights for unpaired positions in t
        for i in t.unpaired():
          try:
            m = max(map(lambda x: x[1], filter(lambda x: x[0] != None, B[i].items())))
          except ValueError:
            # No other probabilities fpr this position
            m = 0
          self.__relativeboltzmannprediction__[i] += _score(m - B[i][None]) / len(target)
    return self.__relativeboltzmannprediction__

# Class providing method for scoring sites based on a weighted mixture
# of incorrect predicted structure, Boltzmann expectation of incorrect
# structure, and a score computed from the difference between the
# Boltzmann probability of the correct structure and the maximum
# Boltzmann probability of an incorrect structure.
class CombinedPrediction(WrongPrediction, BoltzmannPrediction, RelativeBoltzmannPrediction):
  __combinedprediction = None
  def __getcombinedprediction__(self, target = None, wrongweight = 1, boltzmannweight = 1, relativeboltzmannweight = 1, normalise = False):
    if self.__combinedprediction == None:
      self.__combinedprediction = {}
    if (wrongweight, boltzmannweight, relativeboltzmannweight, normalise) not in self.__combinedprediction:
      # Compute positional weights
      if boltzmannweight != 0:
        # Initialise with contributions from boltzmann probability of
        # incorrect structure.
        w = boltzmannweight
        if normalise:
          # Adjust for expected number of wrong predictions
          x = sum(self.__getboltzmannprediction__(target))
          if x != 0:
            w = boltzmannweight / float(x)
        self.__combinedprediction[(wrongweight, boltzmannweight, relativeboltzmannweight, normalise)] = map(lambda x: w * x, self.__getboltzmannprediction__(target))
      else:
        # Initialise as all 0
        self.__combinedprediction[(wrongweight, boltzmannweight, relativeboltzmannweight, normalise)] = len(self) * [0]
      if relativeboltzmannweight != 0:
        # Add contribution from score based on difference between
        # probability of correct structure and probability of most
        # probable incorrect structure.
        w = relativeboltzmannweight
        if normalise:
          # Adjust for sum of relative scores
          x = sum(self.__getrelativeboltzmannprediction__(target))
          if x != 0:
            w = relativeboltzmannweight / float(x)
        for j in xrange(len(self)):
          self.__combinedprediction[(wrongweight, boltzmannweight, relativeboltzmannweight, normalise)][j] += w * self.__getrelativeboltzmannprediction__(target)[j]
      if wrongweight != 0:
        # Add contribution from incorrect position in predicted structure
        w = wrongweight
        if normalise:
          # Adjust for number of wrong predictions
          x = sum(self.__getwrongprediction__(target))
          if x != 0:
            w = wrongweight / float(x)
        for j in xrange(len(self)):
          self.__combinedprediction[(wrongweight, boltzmannweight, relativeboltzmannweight, normalise)][j] += w * self.__getwrongprediction__(target)[j]
      # Normalise by sum of weights - assume sum of weights is
      # non-zero, this class really shouldn't be used with all-zero or
      # negative weights.
      x = wrongweight + boltzmannweight + relativeboltzmannweight
      for i in xrange(len(self)):
        self.__combinedprediction[(wrongweight, boltzmannweight, relativeboltzmannweight, normalise)][i] /= float(x)
    return self.__combinedprediction[(wrongweight, boltzmannweight, relativeboltzmannweight, normalise)]

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
class WrongPosition(WrongPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__getwrongprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      return sum(self.__getwrongprediction__(target))
    else:
      return self.__getwrongprediction__(target)[i]
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, self.pweight(i, target))

# Class choosing a position based on probability of it not having the
# target structure(s) in the Boltzmann ensemble
class BoltzmannPosition(BoltzmannPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__getboltzmannprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      return sum(self.__getboltzmannprediction__(target))
    else:
      return self.__getboltzmannprediction__(target)[i]
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, self.pweight(i, target))

# Class choosing a position based on maximum probabilitiy of it
# not having the target structures in the Boltzmann ensemble
class MinBoltzmannPosition(BoltzmannTuplePrediction):
  def getposition(self, target = None):
    return _weightedselect(map(max, self.__getboltzmanntupleprediction__(target)))
  def pweight(self, i = None, target = None):
    if i == None:
      return sum(map(max, self.__getboltzmanntupleprediction__(target)))
    else:
      return max(self.__getboltzmanntupleprediction__(target))
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, self.pweight(i, target))

# Class choosing a position based on negative log probability of it
# having the correct target structure in the Boltzmann ensemble
class LogBoltzmannPosition(LogBoltzmannPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__getlogboltzmannprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      return sum(self.__getlogboltzmannprediction__(target))
    else:
      return self.__getlogboltzmannprediction__(target)[i]
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, self.pweight(i, target))

# Class choosing position based on score calculated from probability
# of correct structure vs maximum probability of an incorrect
# structure under the Boltzmann distribution.
class RelativeBoltzmannPosition(RelativeBoltzmannPrediction):
  def getposition(self, target = None):
    return _weightedselect(self.__getrelativeboltzmannprediction__(target))
  def pweight(self, i = None, target = None):
    if i == None:
      return sum(self.__getrelativeboltzmannprediction__(target))
    else:
      return self.__getrelativeboltzmannprediction__(target)[i]
  def pweights(self, target = None):
    for i in xrange(len(self.seq)):
      yield (i, self.pweight(i, target))

# Class using a weighted mixture of WrongPosition, BoltzmannPosition,
# and RelativeBoltzmannPosition to define correctness at a
# position. When choosing a position, rather than using the positional
# weight, the contributions are rescaled so the overall contribution
# of each of the three schemes corresponds to their weight, i.e. the
# effect of an overall small sum of positional scores from a scheme is
# eliminated.
def CombinedPosition(wrongweight = 1, boltzmannweight = 1, relativeboltzmannweight = 1):
  class _CombinedPosition(CombinedPrediction):
    __combinedpositionwrongweight = wrongweight
    __combinedpositionboltzmannweight = boltzmannweight
    __combinedpositionrelativeboltzmannweight = relativeboltzmannweight
    def getposition(self, target = None):
      return _weightedselect(self.__getcombinedprediction__(target, self.__combinedpositionwrongweight, self.__combinedpositionboltzmannweight, self.__combinedpositionrelativeboltzmannweight, True))
    def pweight(self, i = None, target = None):
      if i == None:
        return sum(self.__getcombinedprediction__(target, self.__combinedpositionwrongweight, self.__combinedpositionboltzmannweight, self.__combinedpositionrelativeboltzmannweight, False))
      else:
        return self.__getcombinedprediction__(target, self.__combinedpositionwrongweight, self.__combinedpositionboltzmannweight, self.__combinedpositionrelativeboltzmannweight, False)[i]
    def pweights(self, target = None):
      for i in xrange(len(self)):
        yield(i, self.pweight(i, target))
  return _CombinedPosition

### CUT WEIGHTS; classes for selecting cuts for cross overs

# Generic class for selecting cross over cuts. This should be
# specialised to define __weightf.
class CutWeight:
  # Function for computing weight of cuts
  def __computecutweight(self, target = None):
    if target == None:
      target = self.getpopulation().gettarget()
    self.__cweight = {}
    # Compute prefix array of positional weights for efficient
    # computation of cut weights - note that prefix[i] is sum up to
    # index i - 1.
    prefix = [0] * (len(self.seq) + 1)
    for i in xrange(len(self.seq)):
      prefix[i + 1] = prefix[i] + 1 - self._weightf(target)[i]
    # The weight of a cut is sum of one minus positional weight for
    # each position
    for c in target.getcuts():
      # A cut is represented by a set of positions that are pairwise
      # compatible as cuts.
      for i in c:
        for j in c:
          if i < j:
            # Penalise heavily biased cuts
            w = min(1, (j - i) * 10.0 / len(self))
            # Weight of inside structure
            self.__cweight[(i, j)] = w * (prefix[j] - prefix[i])
            # Weight of outside structure
            self.__cweight[(j, i)] = w * (prefix[-1] + prefix[i] - prefix[j])
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
  def cweight(self, i = None, target = None):
    try:
      if self.__cweight == None:
        raise AttributeError
    except AttributeError:
      self.__computecutweight(target)
    if i != None:
      return self.__cweight[i]
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
    return self.__getwrongprediction__(target)

# Class weighting cuts by Boltzmann probabilities of having right structure(s)
class BoltzmannCut(CutWeight, BoltzmannPrediction):
  def _weightf(self, target):
    return self.__getboltzmannprediction__(target)

# Class weigting cuts by minimum Boltzmann probability of having right
# structure, taken over all targets
class MinBoltzmannCut(CutWeight, BoltzmannTuplePrediction):
  def _weightf(self, target):
    return map(lambda x: 1 - max(x), self.__getboltzmanntupleprediction__(target))

# Class weighting cuts by 
class RelativeBoltzmannCut(CutWeight, RelativeBoltzmannPrediction):
  def _weight(self, target):
    return self.__getrelativeboltzmannprediction__(target)

# Class taking weighted average of number of correct positions in
# prediction, Boltzmann expected number of correct predictions, and
# score based on difference between Boltzmann probability of correct
# structure and maximum Boltzmann probability of incorrect structure.
def CombinedCut(wrongweight = 1, boltzmannweight = 1, relativeboltzmannweight = 1):
  class _CombinedCut(CutWeight, CombinedPrediction):
    __combinedcutwrongweight = wrongweight
    __combinedcutboltzmannweight = boltzmannweight
    __combinedcutrelativeboltzmannweight = relativeboltzmannweight
    def _weightf(self, target):
      return self.__getcombinedprediction__(target, self.__combinedcutwrongweight, self.__combinedcutboltzmannweight, self.__combinedcutrelativeboltzmannweight, False)
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
        w = [(current[1] + int(current[0] != offset + 2 * i) * (1 - 2 * current[1])) * self.paired[offset + 2 * i] for i in xrange(2)]
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

### INDIVIDUAL; general class containing the core of a member of the
### GA population.
class Individual:
  def __init__(self, seq = None, population = None, id = None):
    # Initialise attributes capturing sequence and structure information
    self.seq = seq
    self.struc = {}
    self.boltz = {}
    self.fit = None
    self.obj = None
    self.mfe = {}
    self.tfe = {}
    self.population = population
    self.id = id
  # Length of individual is length of its sequence
  def __len__(self):
    if self.seq == None:
      return -1
    else:
      return len(self.seq)
  # Retrieve current sequence of individual
  def getseq(self):
    return int2str(self.seq)
  # Retrieve MFE structure of current sequence at temperature T - if
  # necessary, it is computed
  def getstruc(self, T = None):
    if not self.struc.has_key(T):
      s, self.mfe[T] = fold(self.seq, T)
      self.struc[T] = Structure(s, T)
    return self.struc[T]
  # Retrieve identifier of individual
  def getid(self):
    return self.id
  # Set identifier of individual
  def setid(self, id):
    self.id = id
  # Retrieve Boltzmann probabilities for current sequence at
  # temperature T - if necessary, it is computed
  def getboltzmann(self, T = None):
    if not self.boltz.has_key(T):
      class __Getboltzmann__:
        def __init__(self, p):
          self.probabilities = p
        def __getitem__(self, index):
          if type(index) == IntType:
            return self.probabilities[index][None]
          else:
            if index[1] in self.probabilities[index[0]]:
              return self.probabilities[index[0]][index[1]]
      boltz, s, self.mfe[T] = boltzmann(self.seq, T)
      self.struc[T] = Structure(s, T)
      self.boltz[T] = boltz
    return self.boltz[T]
  # Retrieve minimum free energy taken over all structures for current
  # sequence at temperature T - if necessary, it is computed.
  def getminimumenergy(self, T = None):
    if not self.mfe.has_key(T):
      self.struc[T], self.mfe[T] = fold(self.seq, T)
    return self.mfe[T]
  # Retrieve free energy for structure t at temperature T
  def gettargetenergy(self, t = None, T = None):
    if t == None:
      for i in self.getpopulation().gettarget():
        yield self.gettargetenergy(self.seq, i.bracket(), i.temperature)
    if not self.tfe.has_key((t, T)):
      self.tfe[(t, T)] = energy(self.seq, t, T)
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
  # Retrieve population individual is part of
  def getpopulation(self):
    return self.population
  # Set population individual is part of
  def setpopulation(self, population):
    self.population = population
  # Draw nucleotides for all positions in dependency structure
  # containing idx from dist and update sequence seq accordingly. If
  # position idx is defined (not None), it is assumed to be a position
  # forced to mutate, and draws are done according to a mutation,
  # including effects of current nucleotide. If not defined, draws are
  # done according to an initialisation where no current nucleotide is
  # assumed. In the latter case, all relevant positions in seq are
  # assumed to be None.
  def __drawondependencystructure__(self, dist, seq, target, idx, visited = None):
    deps = list(target.getdependency(idx))
    # Choose new nucleotide for position i
    if len(deps) == 0:
      # Single unpaired position in structure
      seq[idx] = dist.drawunpaired(seq[idx])
    else:
      if len(deps) == 1 and len(target.getdependency(deps[0])) == 1:
        # Dependency structure consists of single base pair
        seq[idx], seq[deps[0]] = int2bp(dist.drawbasepair(seq[idx]))
      else:
        # Complex dependency structure
        # Remember whether we are mutating or drawing from scratch
        if seq[idx] == None:
          unrestricted = True
        else:
          unrestricted = False
        # First choose nucleotide for starting position
        seq[idx] = dist.drawpaired(current = seq[idx])
        if visited == None:
          visited = len(seq) * [None]
        visited[idx] = True
        # Update nucleotides for all neighbours, recursively
        # Start by setting direct neighbours up for update
        stack = [j for j in target.getdependency(idx)]
        for j in stack:
          visited[j] = False
        # As long as there are positions left to update, update them
        # and check their direct neighbours.
        while stack != []:
          # Remove a random element from stack, rather than the top one,
          # to avoid biases from the order in which positions are visited
          j = randint(0, len(stack) - 1)
          if j == len(stack) - 1:
            j = stack.pop()
          else:
            # Rather than delete an element in the middle of the list,
            # move last element to place of element we want to remove.
            stack[j], j = stack.pop(), stack[j]
          # Determine context of selected position
          neighbours = set([])
          for k in target.getdependency(j):
            if visited[k]:
              neighbours.add(seq[k])
            # We might as well also check whether direct neighbours
            # need to be added to positions that need to be updated.
            elif visited[k] is None:
              stack.append(k)
              visited[k] = False
          # Choose new nucleotide for selected position
          if unrestricted:
            seq[j] = dist.drawpaired(neighbours)
          else:
            seq[j] = self.__mutatedependent__(j, dist, neighbours)
          # Mark as visited
          visited[j] = True
    # All done updating the sequence throughout dependency structure

  # Function for mutating in the dependency structure containing
  # position i; if dist is specified, new bases/base pairs are chosen
  # according to this.
  def mutate(self, i, target = None, dist = None):
    if target == None:
      target = self.getpopulation().gettarget()
    if dist == None:
      dist = self.getpopulation().getbasedistribution()
    s = deepcopy(self.seq)
    self.__drawondependencystructure__(dist, s, target, i)
    # Create new individual from updated sequence
    new = self.getpopulation().__individualclass__(s, self.getpopulation())
    if self.id != None:
      new.id = (self.id,)
    #new = GenerateIndividual.__new__(type(self))
    #new.__init__(s)
    return new

  # Function for creating cross-over with other individual. If
  # xoverpoint is a pair of integers, the region from and including
  # the first index and to but not including the second index is taken
  # from other, with the remaining material from self. Otherwise,
  # xoverpoint is assumed to be a set of sets of positions where
  # material is taken from other with the remaining taken from self.
  def crossover(self, other, xoverpoint):
    s = deepcopy(self.seq)
    if len(xoverpoint) == 2 and type(xoverpoint[0]) == IntType:
      # Normal cross-over, switch to other and back again
      for i in xrange(xoverpoint[0], xoverpoint[1]):
        s[i] = other.seq[i]
    else:
      # Cross-over is defined in terms of a set of sets of positions
      # that should be copied from other
      for c in xoverpoint:
        for i in c:
          s[i] = other.seq[i]
    # All done creating cross-over sequence, create and return new
    # object of the same type with this new sequence
    new = self.getpopulation().__individualclass__(s, self.getpopulation())
    if self.id != None or other.id != None:
      new.id = (self.id, other.id)
    #new = GenerateIndividual.__new__(type(self))
    #new.__init__(s)
    return new

  # Function for (re-)initialising the sequence of an individual from
  # scratch, compatible with target. New bases and base pairs are
  # drawn from the distribution specified by dist.
  def initialise(self, target = None, dist = None):
    if target == None:
      target = self.getpopulation().gettarget()
    if dist == None:
      dist = self.getpopulation().getbasedistribution()
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
    self.__init__(self.seq, self.getpopulation())
    self.id = tmp

### Auxiliary classes for computing distance to target
class PredictionErrors(WrongPrediction):
  def __predictionerrors__(self):
    return sum(self.__getwrongprediction__())

class BoltzmannErrors(BoltzmannPrediction):
  def __boltzmannerrors__(self):
    return sum(self.__getboltzmannprediction__())

class RelativeBoltzmannErrors(RelativeBoltzmannPrediction):
  def __relativeboltzmannerrors__(self):
    return sum(self.__getrelativeboltzmannprediction__())

class CombinedErrors(CombinedPrediction):
  def __combinederrors__(self, wrongweight = 1, boltzmannweight = 1, relativeboltzmannweight = 1):
    return sum(self.__getcombinedprediction__(None, wrongweight, boltzmannweight, relativeboltzmannweight, False))

class TresholdBoltzmannErrors(ThresholdBoltzmannPrediction):
  def __thresholdboltzmannerrors__(self, threshold):
    return sum(self.__getthresholdboltzmannprediction__(threshold))

class LogBoltzmannErrors(LogBoltzmannPrediction):
  def __logboltzmannerrors__(self):
    return sum(self.__getlogboltzmannprediction__())

### OBJECTIVE; classes for computing the objective of an individual
# Simple class measuring the number of errors
class Objective(PredictionErrors):
  def objective(self):
    return self.__predictionerrors__()
  def perfectsolution(self):
    return self.objective() == 0

# Class measuring objective as expected number of errors, according to
# Boltzmann distribution.
def BoltzmannObjective(precission = 0):
  class __BoltzmannObjective(BoltzmannErrors):
    __boltzmannobjectiveprecission__ = 1 - precission
    def objective(self):
      return self.__boltzmannerrors__()
    def perfectsolution(self):
      if self.__boltzmannobjectiveprecission__ == 1:
        return self.objective() == 0
      else:
        for t in self.getpopulation().gettarget():
          p = self.getboltzmann(t.temperature)
          # Check paired positions
          for bp in t.basepairs():
            if bp[1] not in p[bp[0]] or p[bp[0]][bp[1]] < self.__boltzmannobjectiveprecission__:
              return False
          # Check unpaired positions
          for i in t.unpaired():
            if p[i] < self.__boltzmannobjectiveprecission__:
              return False
        return True
  return __BoltzmannObjective

# Class measuring objective by difference between Boltzmann
# probability of correct structure and most probable incorrect
# structure.
def RelativeBoltzmannObjective(precission = 0):
  class __RelativeBoltzmannObjective(RelativeBoltzmannErrors):
    __relativeboltzmannobjectiveprecission__ = precission
    def objective(self):
      return self.__relativeboltzmannerrors__()
    def perfectsolution(self):
      for t in self.getpopulation().gettarget():
        # Get Boltzmann probabilities
        B = self.getboltzmann(t.temperature)
        # Check paired positions
        for p in t.basepairs():
          for i in xrange(2):
            if B[p[0]][p[1]] - max([j[1] for j in filter(lambda x: x[0] != p[1 - i], B[p[i]].items())]) < self.__relativeboltzmannobjectiveprecission__:
              return False
        # Check unpaired positions
        for i in t.unpaired():
          if B[i][None] - max([j[1] for j in filter(lambda x: x[0] != None, B[i].items())]) < self.__relativeboltzmannobjectiveprecission__:
            return False
        # All positions have sufficiently high probability for correct
        # structure compared to incorrect structure for all targets.
        return True
  return __RelativeBoltzmannObjective

### FITNESS; classes for computing the fitness of an individual
# Simple class measuring fraction of positions that are incorrect
class Fitness(PredictionErrors):
  def fitness(self):
    return self.__predictionerrors__() / float(len(self))

# Class measuring fitness as fraction of expected number of incorrect
# positions, according to Boltzmann distribution.
class BoltzmannFitness(BoltzmannErrors):
  def fitness(self):
    return self.__boltzmannerrors__() / float(len(self))

# Class measuring fitness as sum of negative log Boltzmann probabilities
class LogBoltzmannFitness(LogBoltzmannErrors):
  def fitness(self):
    return self.__logboltzmannerrors__() / float(len(self))

# Class measuring fitness as fraction of positions not having a
# Boltzmann probability of being correct higher than threshold
def ThresholdBoltzmannFitness(threshold):
  class _ThresholdBoltzmannFitness(ThresholdBoltzmannErrors):
    __thresholdboltzmannfitnessthreshold = threshold
    def fitness(self):
      return self.__thresholdboltzmannerrors__(self.__thresholdboltzmannfitnessthreshold) / float(len(self))

# Class measuring fitness as a score based on Boltzmann probability of
# correct structure compared with maximum probability of incorrect
# structure.
class RelativeBoltzmannFitness(CombinedErrors):
  def fitness(self):
    return self.__relativeboltzmannerrors__() / float(len(self))

# Class measuring fitness by a weighted combination of fraction of
# incorrect positions, Boltzmann expected fraction of incorrect
# positions, and a length normalised sum of scores measuring the
# relative difference between Boltzmann probability of correct
# structure and maximum probability of incorrect structure.
def CombinedFitness(wrongweight = 1, boltzmannweight = 1, relativeboltzmannweight = 1):
  class _CombinedFitness(CombinedErrors):
    __combinedfitnesswrong = wrongweight
    __combinedfitnessboltzmann = boltzmannweight
    __combinedfitnessrelativeboltzmann = relativeboltzmannweight
    def fitness(self):
      return self.__combinederrors__(self.__combinedfitnesswrong, self.__combinedfitnessboltzmann, self.__combinedfitnessrelativeboltzmann) / float(len(self))
  return _CombinedFitness
