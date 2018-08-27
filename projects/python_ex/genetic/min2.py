#encoding=utf-8

# volume.py - useage example
#
# the fittest individual will have a chromosome consisting of 40 '1's
#

from __future__ import division

import genetic as genetic
import math
# using maplotlib
# matplotlib is not ready for python 3.0
#import matplotlib.pyplot as plt

def logN(X, base=math.e, epsilon=1e-12):
  # logN is logarithm function with the default base of e
  integer = 0
  if X < 1 and base < 1:
    raise (ValueError, "logarithm cannot compute")
  while X < 1:
    integer -= 1
    X *= base
  while X >= base:
    integer += 1
    X /= base
  partial = 0.5               # partial = 1/2
  # list = []                   # Prepare an empty list, it seems useless
  X *= X                      # We perform a squaring
  decimal = 0.0
  while partial > epsilon:
    if X >= base:             # If X >= base then a_k is 1
      decimal += partial      # Insert partial to the front of the list
      X = X / base            # Since a_k is 1, we divide the number by the base
    partial *= 0.5            # partial = partial / 2
    X *= X                    # We perform the squaring again
  return (integer + decimal)

#此一加總函式在 volume 最大化中,並未使用
def sum(seq):
    def add(x,y): return x+y
    return reduce(add, seq, 0)

class Min1(genetic.Individual):
    #optimization = genetic.MAXIMIZE
    def evaluate(self, optimum=None):

        self._getvar(self.chromosome)

        x1 = self.var[0]
        #x2 = self.var[1]
        y = self.var[1]


        # 這裡要考量 constrainted 的條件下,設計 penality function

        z=-x1 - logN(x1/2)+y
        t=2*x1-y-logN(x1/2)

        if z <=0:
            p=0
        else:
            p=1000000

        if x1<=1.4 and x1>=0.5:
            d=0
        else:
            d=1000000
        if y<=1 and y>=0:
            e=0
        else:
            e=1000000

        fitness_value = t+p+d+e
        #print w,z

        self.score = fitness_value

    def mutate(self, gene):
        self.chromosome[gene] = not self.chromosome[gene] # bit flip

if __name__ == "__main__":
    env = genetic.Environment(Min1, size=200, maxgenerations=500)
    env.run()
    # matplotlib is not ready for python 3.0 yet
    # use matplotlib to plot the convergence figure
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(env.x, env.y)
    plt.show()
    '''