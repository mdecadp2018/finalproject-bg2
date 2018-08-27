#encoding=utf8

# genetic.py
#

import random
import operator

MAXIMIZE, MINIMIZE = 11, 22

class Individual:
    chromosome = None
    score = None
    # Here the size of var depends on var_number
    var = []
    var_number = 2
    for i in range(var_number):
        var.append(0)

    alleles = (0,1)
    # 以下為參數可負數時的編碼考量
    #前10為小數,後10為整數,第21則為正負號
    #0~9表示小數,10~19表示整數,而指標第20則表示第一數的正號或負號,若為0則表示正,若為1表示負號.
    #21~30表示第二數的小數部分,31~40則表示第二數的整數部分,第41指標則表示第二數的正號或負號
    #42~51表示第三數的小數部分,52~61則表示第二數的整數部分,第62指標則表示第三數的正號或負號
    # -1023 ~ 1023
    #length = 21*var_number,若接受負數參數,則必須同步修改 20->21
    length = 20*var_number
    seperator = ''
    optimization = MINIMIZE

    def __init__(self, chromosome=None):
        self.chromosome = chromosome or self._makechromosome()
        self.score = None  # set during evaluation

    def _getvar(self,chromosome=None):
        x = 0
        for i in range(0,self.var_number):
            for j in range(i*20,i*20+10):
                x +=self.chromosome[j]<<(j-(i*20))
            if (x>999):
                x=999
            x/=1000.
            for j in range(i*20+10,i*20+20):
                x +=self.chromosome[j]<<(j-(i*20+10))
            self.var[i] = x
        return self.var

        ''' for -1023 ~ 1023,當設計變數可以接受負值時使用,每一變數使用21個 bit strings
#for design variable -1023 ~1023
        for i in range(self.var_number):
            x = 0
            for j in range(i*21,i*21+10):
                x +=self.chromosome[j]<<(j-(i*21))
            if (x>999):
                x=999
            x/=1000.
            for j in range(i*(21)+10,i*(21)+20):
                x +=self.chromosome[j]<<(j-(i*21+10))
            if(self.chromosome[i*(21)+20] == 1):
                self.var[i] = -x
            else:
                self.var[i] = x
            x = 0
        return self.var
        '''
    def _makechromosome(self):
        "makes a chromosome from randomly selected alleles."
        return [random.choice(self.alleles) for gene in range(self.length)]

    def evaluate(self, optimum=None):
        "this method MUST be overridden to evaluate individual fitness score."
        pass

    def crossover(self, other):
        "override this method to use your preferred crossover method."
        return self._twopoint(other)

    def mutate(self, gene):
        "override this method to use your preferred mutation method."
        self._pick(gene)

    # sample mutation method
    def _pick(self, gene):
        "chooses a random allele to replace this gene's allele."
        self.chromosome[gene] = random.choice(self.alleles)

    # sample crossover method
    def _twopoint(self, other):
        "creates offspring via two-point crossover between mates."
        left, right = self._pickpivots()
        def mate(p0, p1):
            chromosome = p0.chromosome[:] # 交配時,以p0的基因為基礎(複製整個 p0 的染色體內容
            chromosome[left:right] = p1.chromosome[left:right] # 接續上一個 p0 的染色體內容,將索引 left 至 right 的內容,替換成 p1 的基因
            #child = p1.__class__(chromosome) 這是原先的程式,但是應該子代要指向 p0 的內容才對
            child = p0.__class__(chromosome)
            child._repair(p0, p1)
            return child
        return mate(self, other), mate(other, self)

    # some crossover helpers ...
    def _repair(self, parent1, parent2):
        "override this method, if necessary, to fix duplicated genes."
        pass

    def _pickpivots(self):
        left = random.randrange(1, self.length-2)
        right = random.randrange(left, self.length-1)
        return left, right

    #
    # other methods
    #

    def __repr__(self):
        "returns string representation of self"
        '''
        return '<%s chromosome="%s" score=%s var=%s>' % \
               (self.__class__.__name__,
                self.seperator.join(map(str,self.chromosome)), self.score,self._getvar(self.chromosome))
        '''
        return '<%s score=%s var=%s>' % \
               (self.__class__.__name__,self.score,self._getvar(self.chromosome))
    # since the __cmp__ special function is gone  use the __lt__ in stead
    # use the expression (a > b) - (a < b) as the equivalent for cmp(a, b)
    #def __cmp__(self, other):
    # these are for python 3
    def __cmp__(self, other):
        if self.optimization == MINIMIZE:
            #return cmp(self.score, other.score)
            return (self.score > other.score) - (self.score < other.score)
        else: # MAXIMIZE
            #return cmp(other.score, self.score)
            return (other.score > self.score) - (other.score < self.score)
            
    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __le__(self, other):
        return self.__cmp__(other) <= 0

    def __gt__(self, other):
        return self.__cmp__(other) > 0

    def __ge__(self, other):
        return self.__cmp__(other) >= 0 

    def copy(self):
        twin = self.__class__(self.chromosome[:])
        twin.score = self.score
        return twin


class Environment(object):
    x = [0]
    y = [0]
    def __init__(self, kind, population=None, size=100, maxgenerations=100,
                 crossover_rate=0.90, mutation_rate=0.07, optimum=None):
        self.kind = kind
        self.size = size
        self.optimum = optimum
        self.population = population or self._makepopulation()
        for individual in self.population:
            individual.evaluate(self.optimum)
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.maxgenerations = maxgenerations
        self.generation = 0
        self.report()

    def _makepopulation(self):
        return [self.kind() for individual in range(self.size)]

    def run(self):
        while not self._goal():
            self.step()

    def _goal(self):
        return self.generation > self.maxgenerations or \
               self.best.score == self.optimum

    def step(self):
        # this sort is not working with python 3.0, modification is needed
        self.population.sort()
        self._crossover()
        self.generation += 1
        self.report()
        self.x.append(self.generation)
        # 設定為只附加所選定範圍的值,這裡只取大於或等於 0 的 score 值
        if self.best.score <=5:
            self.y.append(self.best.score)
        else:
            self.y.append(5)

    def _crossover(self):
        next_population = [self.best.copy()]
        while len(next_population) < self.size:
            mate1 = self._select()
            if random.random() < self.crossover_rate:
                mate2 = self._select()
                offspring = mate1.crossover(mate2)
            else:
                offspring = [mate1.copy()]
            for individual in offspring:
                self._mutate(individual)
                individual.evaluate(self.optimum)
                next_population.append(individual)
        self.population = next_population[:self.size]

    def _select(self):
        "override this to use your preferred selection method"
        return self._tournament()

    def _mutate(self, individual):
        for gene in range(individual.length):
            if random.random() < self.mutation_rate:
                individual.mutate(gene)

    #
    # sample selection method
    #
    def _tournament(self, size=8, choosebest=0.90):
        competitors = [random.choice(self.population) for i in range(size)]
        competitors.sort()
        if random.random() < choosebest:
            return competitors[0]
        else:
            return random.choice(competitors[1:])

    def best():
        doc = "individual with best fitness score in population."
        def fget(self):
            return self.population[0]
        return locals()
    best = property(**best())

    def report(self):
        print ("="*70)
        print ("generation: ", self.generation)
        print ("best:       ", self.best)
