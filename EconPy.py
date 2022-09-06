' Engineering Economy Formulas '

from math import e

'''
Chapter 2 Fast Forms
'''
ie = lambda i, m : 100*(pow(1 + (i/100) / m, m) - 1) # discrete compound periods
iec = lambda i : (e**(i/100) - 1)*100 # continuous compounding
sld = lambda P, S, N : (P - S) / N # straight line dep. amount
sldbv = lambda n, P, S, N : P - n*((P - S) / N) # book value after n years of SLD
sldAcc = lambda n, P, S, N : n*((P - S) / N) # accumulated SLD
dbdbv = lambda P, d, n : P*pow(1 - d/100, n) # declining balance depreciation
dbdn = lambda P, d, n : dbdbv(P, d, n - 1)*d/100 # depreciation in nth year

'''
Chapter 2 Print Forms
'''
# Effective Interest Rate
def ie_p(i, compoundPeriod, cout = "yes"):
    i = i / 100
    if compoundPeriod == "continuous":
        ie = e**i - 1
        if cout == "yes":
            print(f"ie ({i}%, {compoundPeriod}): {ie*100}")
        return ie*100
    elif compoundPeriod == "quarterly":
        m = 4
    elif compoundPeriod == "monthly":
        m = 12
    elif compoundPeriod == "weekly":
        m = 52
    elif compoundPeriod == "daily":
        m = 365
    ie = pow(1 + i / m, m) - 1
    if cout == "yes":
        print(f"ie ({i*100}%, {compoundPeriod}): {ie*100}")
    return ie*100

# Straight Line Depreciation
def SLD_p(P, S, N, cout = "yes"):
    sld = (P - S) / N
    if cout == "yes":
        print(f"SLD (${P}, ${S}, {N}): {sld}\n")
    return sld

# Straight Line Book Value
def SLDBV_p(n, P, S, N, cout = "yes"):
    sldbv = P - n*((P - S) / N)
    if cout == "yes":
        print(f"SLDBV ({n}, ${P}, ${S}, {N}): {sldbv}\n")
    return sldbv

# Straight Line Accumulated Depreciation
def SLDAcc_p(n, P, S, N, cout = "yes"):
    sldacc = n*((P - S) / N)
    if cout == "yes":
        print(f"SLDBV ({n}, ${P}, ${S}, {N}): {sldacc}\n")
    return sldacc

# Declining Balance Book Value
def DBDBV_p(P, d, n, cout = "yes"):
    d = d / 100
    dbdbv = P*pow(1 - d, n)
    if cout == "yes":
        print(f"DBDBV (${P}, {d*100}%, {n}): {dbdbv}\n")
    return dbdbv

# Declining Depreciation for Period n
def DBDdn_p(P, d, n, cout = "yes"):
    d = d / 100
    dbdn = DBDBV_p(P, d, n - 1)*d
    if cout == "yes":
        print("DBDdn (${P}, {d*100}%, {n}): {dbdn}\n")
    return dbdn

'''
Chapter 3 Fast Forms
'''
FP = lambda i, N : pow(1+(i / 100), N) # Compound Amount Factor
PF = lambda i, N : 1 / FP(i, N) # Present Worth Factor
AF = lambda i, N : (i / 100) / (FP(i, N) - 1) # Sinking Fund Factor
FA = lambda i, N : 1 / AF(i, N) # Uniform Series Compund Amount Factor
AP = lambda i, N : AF(i, N)*FP(i, N) # Capital Recovery Factor
PA = lambda i, N : 1 / AP(i, N) # Series Present Worth Factor
AG = lambda i, N : (1 / (i / 100)) - N / (FP(i, N) - 1) # Arithmetic Gradient
A = lambda P, S, i, N : (P-S)*AP(i, N)+S*i/100 # Capitalized Recovery Formula
def PAg(g, i, N): # Geometric Gradient
    if i != g:
        i = i / 100
        g = g / 100
        iO = ((1+i) / (1+g)) - 1
        print(f"iO: {iO*100}%")
        factor = ((pow(1+iO, N) - 1) / (iO*pow(1+iO, N))) * (1 / (1 + g))
        return factor
    else:
        g = g / 100
        factor = N / (1+g)
        return factor
bondP = lambda FV, CR, MR, N, m : FV*PF(MR/m, N*m)+(FV*0.01*CR/m)*PA(MR/m, N*m) # Present bond value
'''
Chapter 3 Print Forms
'''

# Find i or N for a particular function
def iterator(func, factor, error = 0.0001, i = 0, N = 0, roots=1):
    ans = 0.1
    ansArr = []
    if roots == 1:
        if i == 0 and N != 0:  
            while abs((func(ans, N) - factor) / func(ans, N)) > error:
                ans += 0.0001
        elif i != 0 and N == 0:
            while abs((func(i, ans) - factor) / func(i, ans))  > error:
                ans += 0.0001
        return ans
    else:
        error = 1
        for ans in range(1000):
            if abs(func(ans, N) - factor) < error:
                ansArr.append(ans)
        return ansArr

# Compound Amount Factor
def FP_p(i, N, P = 0, cout = "yes"):
    # Solve for F/P or F
    if i != 0 and N != 0:
        i = i / 100
        factor = pow(1+i, N)
        if cout == "yes":
            print(f"F/P ({i*100}%, {N}): {factor}\n")
        if P != 0:
            if cout == "yes":
                print(f"F ({i*100}%, {N}, ${P}): {factor * P }\n")
            return factor * P
        else:
            return factor

# Present Worth Factor
def PF_p(i, N, F = 0, cout = "yes"):
    # Solve for P/F or P
    if i != 0 and N != 0:
        i = i / 100
        factor = 1 / pow(1+i, N)
        if cout == "yes":
            print(f"P/F ({i*100}%, {N}): {factor}\n")
        if F != 0:
            if cout == "yes":
                print(f"P ({i*100}%, {N}, ${F}): {factor * F}\n")
            return factor * F
        else:
            return factor

# Sinking Fund Factor
def AF_p(i, N, F = 0, cout = "yes"):
    # Solve for A/F or A
    if i != 0 and N != 0:
        i = i / 100
        factor = i / (pow(1+i, N) - 1)
        if cout == "yes":
            print(f"A/F ({i*100}%, {N}): {factor}\n")
        if F != 0:
            if cout == "yes":
                print(f"A ({i*100}%, {N}, ${F}): {factor * F }\n")
            return factor * F
        else:
            return factor

# Uniform Series Compund Amount Factor
def FA_p(i, N, A = 0, cout = "yes"):
    # Solve for F/A or F
    if i != 0 and N != 0:
        i = i / 100
        factor = (pow(1+i, N) - 1) / i
        if cout == "yes":
            print(f"F/A ({i*100}%, {N}): {factor}\n")
        if A != 0:
            if cout == "yes":
                print(f"F ({i*100}%, {N}, ${A}): {factor * A }\n")
            return factor * A
        else:
            return factor
            
# Capital Recovery Factor
def AP_p(i, N, P = 0, cout = "yes"):
    # Solve for F/A or F
    if i != 0 and N != 0:
        i = i / 100
        factor = i*pow(1+i, N) / (pow(1+i, N) - 1)
        if cout == "yes":
            print(f"A/P ({i*100}%, {N}): {factor}\n")
        if P != 0:
            if cout == "yes":
                print(f"A ({i*100}%, {N}, ${P}): {factor * P}\n")
            return factor * P
        else:
            return factor

# Series Present Worth Factor
def PA_p(i, N, A = 0, cout = "yes"):
    # Solve for F/A or F
    if i != 0 and N != 0:
        i = i / 100
        factor = (pow(1+i, N) - 1) / (i*pow(1+i, N))
        if cout == "yes":
            print(f"P/A ({i*100}%, {N}): {factor}\n")
        if A != 0:
            if cout == "yes":
                print(f"P ({i*100}%, {N}, ${A}): {factor * A}\n")
            return factor * A
        else:
            return factor

# Arithmetic Gradient to Annuity Conversion Factor
def AG_p(i, N, G = 0):
    i = i / 100
    factor = (1 / i) - N / (pow(1 + i, N) - 1)
    print(f"A/G ({i*100}%, {N}): {factor}\n")
    if G != 0:
        print(f"A ({i*100}%, {N}, ${G}): {factor * G}\n")
        return factor * G
    else:
        return factor

# Geometric Gradient to Present Worth Conversion Factor
def PAg_p(g, i, N, A = 0):
    i = i / 100
    g = g / 100
    iO = ((1+i) / (1+g)) - 1
    factor = ((pow(1+iO, N) - 1) / (iO*pow(1+iO, N))) * (1 / (1 + g))
    print(f"PAg ({i*100}%, {N}): {factor}\n")
    return factor

'''
Project Management
'''

indexA = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "E" : 4, "F" : 5, "G" : 6,
         "H" : 7, "I" : 8, "J" : 9, "K" : 10, "L" : 11, "M" : 12, "N" : 13,
         "O" : 14, "P" : 15, "Q" : 16, "R" : 17, "S" : 18, "T" : 19, "U" : 20,
         "V" : 21, "W" : 22, "X" : 23, "Y" : 24, "Z" : 25}
index1 = {"1" : 0, "2" : 1, "3" : 2, "4" : 3, "5" : 4, "6" : 5, "7" : 6,
         "8" : 7, "9" : 8, "10" : 9, "11" : 10, "12" : 11, "13" : 12, "14" : 13,
         "15" : 14, "16" : 15, "17" : 16}
keysA = list(indexA.keys())
keys1 = list(index1.keys())

class PERT_Object:
    def __init__(self, master, predecessors, length, index):
        self.master = master
        self.predecessors = predecessors
        self.length = length
        # Need to be updated
        self.index = index
        self.limitingFactors = []
        self.elapsedTime = self.length
        self.backwardsTime = []
        self.slack = 0

class PERT_system:
    def __init__(self):
        self.objectCount = 0
        self.objects = []
        self.criticalPath = []
        self.forwardTime = 0
        self.index = {}
        self.keys = []
        self.loadCount()
        self.loadObjects()
        self.findElapsedTimes()
        self.findCompletionTime()
        self.findCriticalPath()
        self.backwardsTime()

    def loadCount(self):
        lastProcess = input("Last Process in System: ").upper()
        if lastProcess.isalpha():
            self.index = indexA
            self.keys = keysA
        else:
            self.index = index1
            self.keys = keys1
        self.objectCount = self.index[lastProcess]+1
        
    def loadObjects(self):
        for i in range(self.objectCount):
            print(f"Process {self.keys[i]}")
            predecessorStr = input("Predecessors: ").upper()
            predecessors = predecessorStr.split()
            lengthStr = input("Process Length: ")
            length = lengthStr.split()
            if len(length) > 1:
                length = (int(length[0])+4*int(length[1])+int(length[2]))/6
            else:
                length = int(length[0])
            #length = int(input("Process Length: "))
            print("=====")
            new_object = PERT_Object(self, predecessors, length, self.index)
            self.objects.append(new_object)
        
    def findElapsedTimes(self):
        for i in range(self.objectCount):
            iterCnt = len(self.objects[i].predecessors)
            if iterCnt > 0:
                maxLen = 0
                maxIndex = 0
                state = 0
                for j in range(iterCnt):
                    currentIndex = self.index[self.objects[i].predecessors[j]]
                    currentItem = self.objects[currentIndex]
                    if currentItem.elapsedTime > maxLen:
                        maxLen = currentItem.elapsedTime
                        maxIndex = currentIndex
                        state = j
                self.objects[i].elapsedTime += maxLen
                self.objects[i].limitingFactors = self.objects[maxIndex].limitingFactors + [self.objects[i].predecessors[state]]
                
    def findCompletionTime(self):
        for i in range(self.objectCount):
            if self.objects[i].elapsedTime > self.forwardTime:
                self.forwardTime = self.objects[i].elapsedTime
        print(f"Completion Time: {str(self.forwardTime)}")
            
    def findCriticalPath(self):
        for i in range(0,self.objectCount):
            if self.objects[i].elapsedTime == self.forwardTime:
                self.criticalPath = self.objects[i].limitingFactors
                self.criticalPath.append(self.keys[i])
        print(f"Critical Path: {self.criticalPath}")

    def backwardsTime(self):
        print("=====\nSlack:")
        
        cnt = self.objectCount - 1
        
        for i in range(len(self.objects[cnt].predecessors)):
            updateIndex = self.index[self.objects[cnt].predecessors[i]]
            self.objects[updateIndex].backwardsTime.append(self.objects[cnt].elapsedTime - self.objects[cnt].length)
            
        for i in range(cnt):
            if len(self.objects[cnt-i-1].predecessors) > 0:
                    if self.objects[cnt-i-1].backwardsTime != []:
                        for j in range(len(self.objects[cnt-i-1].predecessors)):
                            updateIndex = self.index[self.objects[cnt-1-i].predecessors[j]]
                            self.objects[updateIndex].backwardsTime.append(min(self.objects[cnt-1-i].backwardsTime) - self.objects[cnt-1-i].length)
                    else:
                        for j in range(len(self.objects[cnt-i].predecessors)):
                            updateIndex = self.index[self.objects[cnt-i].predecessors[j]]
                            self.objects[updateIndex].backwardsTime.append(self.objects[cnt-i].elapsedTime - self.objects[cnt-i].length)
                    
        for i in range(cnt+1):
            if self.objects[i].backwardsTime != []:
                self.objects[i].slack = min(self.objects[i].backwardsTime) - self.objects[i].elapsedTime
            print(f"{self.keys[i]}: {self.objects[i].slack}")

def avgTime(t0, tm, tp):
    return (t0 + 4*tm + tp)/6
