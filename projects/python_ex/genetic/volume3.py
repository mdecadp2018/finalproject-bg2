#encoding=utf8

# volume.py - useage example
#
# the fittest individual will have a chromosome consisting of 40 '1's
#

import genetic

#此一加總函式在 volume 最大化中,並未使用
def sum(seq):
    def add(x,y): return x+y
    return reduce(add, seq, 0)

class Volume(genetic.Individual):
    optimization = genetic.MAXIMIZE
    def evaluate(self, optimum=None):
        SURFACE = 80
        # self.score is the fitness value
        # x 由0.001~1023.999,利用二位元的移位運算子<<進行計算
        # 第一階段,根據變數個數,自行計算整數與小數加總
        '''
        for j in range(10):
            x +=self.chromosome[j]<<(j-0)
        if (x>999):
            x=999
        x/=1000.
        for j in range(10,20):
            x +=self.chromosome[j]<<(j-10)
        # calculate y
        for j in range(20,30):
            y +=self.chromosome[j]<<(j-20)
        if (y>999):
            y=999
        y/=1000.
        for j in range(30,40):
            y +=self.chromosome[j]<<(j-30)
            
        z=(SURFACE - x*y)/(2.*(x+y))
        fitness_value = x*y*z
        #self.score = sum(self.chromosome[10:20])
        self.score = fitness_value
        #第二階段,將各變數計算以變數個數,利用迴圈進行運算
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
        '''
        #第三階段,在實用上,將上述運算,寫進個別成員類別中,透過呼叫成員函式,取得各變數值
        self._getvar(self.chromosome)
        
        x = self.var[0]
        y = self.var[1]
        z=(SURFACE - x*y)/(2.*(x+y))
        fitness_value = x*y*z
        
        self.score = fitness_value
        
    def mutate(self, gene):
        self.chromosome[gene] = not self.chromosome[gene] # bit flip
   
if __name__ == "__main__":
    env = genetic.Environment(Volume, size=500, maxgenerations=100)
    env.run()
