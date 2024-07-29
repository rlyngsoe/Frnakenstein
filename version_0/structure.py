from copy import deepcopy
from types import IntType, ListType, TupleType, StringType

### Mapping from integers to bases and base pairs, and from integer
### representations of base pairs to pairs of integer representations
### of bases.
BASES = "ACGUT"
def BASESindex(x):
  return min(BASES.index(x), 3)
BASEPAIRS = ["AU", "CG", "GC", "UA", "GU", "UG"]

def int2str(s):
  return "".join(map(lambda x: BASES[x], s))

def str2int(s):
  return map(lambda x: BASESindex(x.upper()), s)

def int2bp(a):
  if a < 4:
    # Canonical base pair
    return (a, 3 - a)
  else:
    # Wobble base pair
    return (a - 2, 7 - a)

def bp2int(a):
  return 2 * a[0] + a[1] - 3

# Base pairs differing between two structures
def bpdif(t, u):
  return t.basepairs().symmetric_difference(u.basepairs())

# Base pair distance between two structures
def bpdist(t, u):
  return len(bpdif(t, u))

# Positions where two stuctures differ
def posdif(t, u):
  return reduce(lambda x, y: x.union(y), t.basepairs().symmetric_difference(u.basepairs()), set([]))

# Distance measured as number of positions where structures differ
def posdist(t, u):
  return len(posdif(t, u))

# Function generating a merge of two bracket notations:
#  .() : structures agree
#  <>  : . in s, () in t
#  {}  : () in s, . in t
#  []  : () in both s and t, but paired with different positions
#  \/  : () in s, )( in t
def combinestructures(s, t):
  # Return array of pairing position, None indicating unpaired position
  def pairing(s):
    p = len(s) * [None]
    stack = []
    for i in xrange(len(s)):
      if s[i] == '(':
        stack.append(i)
      elif s[i] == ')':
        j = stack.pop()
        p[j] = i
        p[i] = j
    return p
  sp = pairing(s)
  tp = pairing(t)
  c = []
  for i in xrange(len(s)):
    if sp[i] == tp[i]:
      # Match between structures
      c.append(s[i])
    elif s[i] == '.':
      # Is unpaired, should have been paired
      c.append({'(': '<', ')': '>'}[t[i]])
    elif s[i] == t[i]:
      # Base paired the right way, but to wrong position
      c.append({'(': '[', ')': ']'}[s[i]])
    elif t[i] == '.':
      # Is base paired, should have been unpaired
      c.append({'(': '{', ')': '}'}[s[i]])
    else:
      # Base paired, but in wrong direction
      c.append({'(': '\\', ')': '/'}[s[i]])
  return "".join(c)

# Class for maintaining a structure in various formats
class Structure:
  # Function for initialising a structure
  def __init__(self, struct = None, T = None):
    if struct != None:
      if type(struct) == StringType:
        # Structure passed as bracket notation
        self.__bracket = struct
      elif type(struct) in [ListType, TupleType] and len(struct) == 2 and type(struct[0]) == IntType and "__iter__" in dir(struct[1]):
        # Structure passed as pair of sequence length and set/list of base pairs
        self.__length = struct[0]
        self.__basepairs = struct[1]
    self.temperature = T
  # Length of structure is number of positions
  def __len__(self):
    if "_Structure__length" not in dir(self):
      self.__length = len(self.bracket())
    return self.__length
  # Validate that structure information is well formed
  def validate(self):
    # Check validity of bracket representation
    try:
      if self.__bracket != None:
        c = 0
        for i in self.__bracket:
          if i == '(':
            c += 1
          elif i == ')':
            c -= 1
            if c < 0:
              # More right parentheses than left parentheses at this point
              raise ValueError, "Structure " + self.__bracket + " is not well formed"
        if c > 0:
          # More left parentheses than right parentheses in structure
          raise ValueError, "Structure " + self.__bracket + " is not well formed"
    except AttributeError:
      # Structure does not have bracket representation
      pass
    # Check validity of basepair representation
    try:
      for p in self.__basepairs:
        for i in p:
          if i < 0 or i >= self.__length:
            # Position i is outside range from 0 to length - 1
            raise ValueError, "Structure " + str(self.__basepairs) + " does not have all base paired positions within range of its length of " + str(self.__length)
    except AttributeError:
      # Structure does not have basepair representation
      pass

  # Function for returning structure in bracket format - this should
  # be the go-to function for all other structure formats.
  def bracket(self):
    try:
      if self.__bracket == None:
        raise AttributeError
    except AttributeError:
      if "_Structure__basepairs" in dir(self) and self.__basepairs != None:
        b = self.__length * ['.']
        for i in self.__basepairs:
          b[i[0]] = '('
          b[i[1]] = ')'
        self.__bracket = "".join(b)
      else:
        # No source of structure information
        raise LookupError, "Attempt to look up structure that is undefined"
    return self.__bracket
  # Function for returning structure as set of base pairs
  def basepairs(self):
    try:
      if self.__basepairs == None:
        raise AttributeError
    except AttributeError:
      try:
        self.__basepairs = set([])
        stack = []
        s = self.bracket()
        self.__length = len(s)
        for i in xrange(len(s)):
          if s[i] == '(':
            stack.append(i)
          elif s[i] == ')':
            self.__basepairs.add((stack.pop(), i))
        if stack != []:
          raise IndexError
      except IndexError:
        raise ValueError, "Invalid bracket representation:" + s
    return self.__basepairs
  # Function for returning unpaired positions in structure - this does
  # not uniquely define the structure and is not considered a valid
  # structure representation.
  def unpaired(self):
    u = set(range(len(self)))
    for b in self.basepairs():
      u.difference_update(b)
    return u

# Class for maintaining target consisting of one or more structures
class Target:
  # Initialisation of target - structures should be list of target structures
  def __init__(self, structures = None):
    if structures == None:
      self.structures = []
    else:
      self.structures = [s for s in structures]
  # Length of a target is number of structures it contains
  def __len__(self):
    return len(self.structures)
  # Size of a target is number of positions it contains
  def size(self):
    if self.structures == []:
      return -1
    else:
      return len(self.structures[0])
  # Validate that target consists of well formed structures of the same length
  def validate(self):
    for s in self.structures:
      s.validate()
      if len(s) != len(self.structures[0]):
        raise ValueError, "Target structures not of the same length"
  # Iterator for target structures
  def __iter__(self):
    if self.structures != None:
      for s in self.structures:
        yield s
  # Add one more structure to target
  def addstructure(self, structure):
    self.structures.append(structure)
    # Previously computed dependencies and cuts are now void
    self.dependencies = None
  # Compute the dependency classes imposed by the structures in this target
  def computedependencies(self):
    if len(self.structures) == 0:
      raise RuntimeError, "Cannot build dependency structure with no targets"
    n = len(self.structures[0].bracket())
    self.dependencies = [set([]) for i in xrange(n)]
    for t in self.structures:
      for bp in t.basepairs():
        self.dependencies[bp[0]].add(bp[1])
        self.dependencies[bp[1]].add(bp[0])
    # Check consistency and separate into components
    self.components = []
    colour = n * [False]
    for i in xrange(n):
      if colour[i] is False:
        # Node i is not part of previously visited component; assign
        # it first colour and start building a new component
        stack = [i]
        colour[i] = 1
        self.components.append(set([]))
        # Iterate through elements in component
        while stack != []:
          j = stack.pop()
          self.components[-1].add(j)
          # Check neighbours to this element
          for k in self.dependencies[j]:
            if colour[k] is False:
              # First time we see this element, colour it with
              # opposite colour of current element and insert it as
              # one that has to be visited
              colour[k] = 3 - colour[j]
              stack.append(k)
            else:
              # We have seen this element before, make sure that its
              # colour is consistent with current element
              if colour[k] == colour[j]:
                raise ValueError, "Targets impose incompatible restrictions on sequence"
  # Retrieve dependency component containing position i
  def getdependency(self, i):
    if "dependencies" not in dir(self) or self.dependencies == None:
      self.computedependencies()
    return self.dependencies[i]
  # Iterator for all dependency components
  def getdependencies(self):
    if "components" not in dir(self) or self.dependencies == None:
      self.computedependencies()
    return self.components
  # Retrieve all crossover cuts compatible with the dependencies. A
  # cut is given by the positions before which the sequence changes,
  # and is represented as sets from which any two elements constitute
  # a cut compatible with the dependencies.
  def getcuts(self):
    if "cuts" not in dir(self) or self.dependencies == None:
      # Build cut classes - initially every position is in the same class
      self.cuts = set([tuple(range(len(self.dependencies) + 1))])
      # Make sure we have some dependencies to compute the cuts from
      if "dependencies" not in dir(self) or self.dependencies == None:
        self.computedependencies()
      # Each dependency break classes that have positions both inside
      # and outside the base pair
      for i in xrange(len(self.dependencies)):
        for j in self.dependencies[i]:
          if i < j: # Each base pair is registered in two places
            # We can't delete and add to classes while we iterate
            # through its members, so store information about classes
            # to be deleted and to be added
            new = []
            old = []
            for A in self.cuts:
              # For each class, determine positions inside and outside
              # current dependency
              inside = []
              outside = []
              for k in A:
                if k > i and k <= j:
                  inside.append(k)
                else:
                  outside.append(k)
              if len(inside) > 0 and len(outside) > 0:
                # Class split into two classes
                old.append(A)
                new.append(tuple(inside))
                new.append(tuple(outside))
            # Update cut classes with the new splits
            for A in old:
              self.cuts.remove(A)
            for A in new:
              self.cuts.add(A)
    # Now we know we have cut classes
    return self.cuts
