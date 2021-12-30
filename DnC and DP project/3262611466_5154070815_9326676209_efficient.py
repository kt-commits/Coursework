import sys
import tracemalloc
import time
def writeToFile(text):
    f=open("output.txt","w")
    f.write(text)
    f.close()
    return

def generateString(filename=sys.argv[1]):
    if filename is None or filename=='':
        filename="input.txt"
    f = open(filename,"r")
    sequence1 = f.readline().rstrip()
    index = ""
    while True:
        index = f.readline().rstrip()
        if  not index.isnumeric():
            break
        sequence1 = sequence1[0:int(index) + 1] + sequence1 + sequence1[int(index) + 1:]
    
    sequence2 = index
    while True:
        index = f.readline().rstrip()
        if  not index.isnumeric():
            break
        sequence2 = sequence2[0:int(index) + 1] + sequence2 + sequence2[int(index) + 1:]
    f.close()
    return sequence1, sequence2


def getMinimumPenalty(x : str,y : str):
    i , j = 0, 0
    m, n = len(x), len(y)
    opt = [[0 for _ in range (n + 1)] for _ in range(m + 1)]
    for i in range(m + 1):
        opt[i][0] = i * delta
    
    for j in range(n + 1):
        opt[0][j] = j * delta
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if x[i - 1]  == y[j - 1]:
                opt[i][j] = opt[i - 1][j -1]
            
            else:
                opt[i][j] = min(opt[i - 1][j - 1] + alpha[x[i -1]+y[j-1]], opt[i - 1][j] + delta, opt[i][j - 1] + delta)
    l = m + n
    i = m
    j = n
    x_position, y_position = l, l
    x_final = ['']*(l+1)
    y_final = ['']*(l+1)
    
    while not(i == 0 or j == 0):
        if x[i-1] == y[j-1]:
            x_final[x_position] = x[i - 1]
            y_final[y_position] = y[j - 1]
            x_position -= 1
            y_position -= 1
            i -= 1
            j -= 1
        elif opt[i - 1][j - 1] + alpha[x[i-1] + y[j-1]] == opt[i][j]:
            x_final[x_position] = x[i - 1]
            y_final[y_position] = y[j - 1]
            x_position -= 1
            y_position -= 1
            i -= 1
            j -= 1
        elif opt[i-1][j] + delta == opt[i][j]:
            x_final[x_position] = x[i - 1]
            y_final[y_position] = "_"
            x_position -= 1
            y_position -= 1
            i -= 1
        elif opt[i][j - 1]  + delta == opt[i][j]:
            y_final[y_position] = y[j - 1]
            x_final[x_position] = "_"
            x_position -= 1
            y_position -= 1
            j -= 1       
    
    while x_position > 0: 
        if i > 0:
            x_final[x_position] = x[i-1]
            x_position -= 1
            i -=1
        else:
            x_final[x_position] = "_"
            x_position -= 1
        
    while y_position > 0:
        if j > 0:
            y_final[y_position] = y[j-1]
            y_position -= 1
            j -=1
        else:
            y_final[y_position] = "_"
            y_position -= 1
    index = 1
    while x_final[index] == '_' and y_final[index] == '_':
        index = index + 1
    
    s1 = ''.join(x_final[index:])
    s2 = ''.join(y_final[index:])
    return s1, s2, opt[m][n]


def spaceEfficientAlignment(x,y):
    m = len(x)
    n = len(y)
    opt = [[0]*(n+1) for i in range(2)]
    for i in range(n+1):
        opt[0][i] = i * delta
    opt[1][0] = delta 
    for i in range(1,m+1):
        for j in range(1,n+1):
            if x[i-1] == y[j-1]:
                opt[1][j] = opt[0][j-1]
            else:
                opt[1][j] = min(alpha[x[i-1] + y[j-1]] + opt[0][j-1], delta + opt[0][j], delta + opt[1][j-1])
    
        for i in range(n+1):
            opt[0][i] = opt[1][i]
        opt[1][0] = opt[1][0]+delta
    return opt[-1]


def divideAndConquer(x ,y):
    if len(x) <= 2 or len(y)<= 2:
        return getMinimumPenalty(x, y)
    
    mid_x = len(x)//2
    x_left_partition = x[:mid_x]
    x_right_partition = x[mid_x:]
    y_left_penalty = spaceEfficientAlignment(x_left_partition, y)
    y_right_penalty = spaceEfficientAlignment(x_right_partition[::-1],y[::-1])
    pc = []
    for i, j in zip(y_left_penalty, y_right_penalty[::-1]):
        pc.append(i+j)
    
    y_mid = pc.index(min(pc))
    s1x, s1y, cost_left_partition = divideAndConquer(x_left_partition, y[:y_mid])
    s2x, s2y, cost_right_partition = divideAndConquer(x_right_partition, y[y_mid:])
    return (s1x + s2x, s1y + s2y, cost_left_partition + cost_right_partition)
if __name__ == "__main__":
    tracemalloc.start()
    gene1, gene2 = generateString()
    time_start=time.time()
    delta = 30
    alpha = {
        "AC":110,
        "AG":48,
        "AT":94,
        "CA":110,
        "CG":118,
        "CT":48,
        "GA":48,
        "GC":118,
        "GT":110,
        "TA":94,
        "TC":48,
        "TG":110,
        "AA":0,
        "TT":0,
        "GG":0,
        "CC":0
    }
    a,b,c = divideAndConquer(gene1, gene2)
    _, peak = tracemalloc.get_traced_memory()    
    writeToFile(a[:50]+' '+a[-50:]+'\n'+b[:50]+' '+b[-50:]+'\n'+str(c)+'\n'+str(time.time()-time_start)+'\n'+str(peak / 10**3)+'\n')    
    tracemalloc.stop()