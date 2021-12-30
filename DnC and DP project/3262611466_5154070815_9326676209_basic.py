import tracemalloc
import time
import sys
def writeToFile(text):
    f=open("output.txt","w")
    f.write(text)
    f.close()
    del f
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
    f=None
    return sequence1, sequence2
    
            
def getMinimumPenalty(x : str,y : str, alpha: int, delta: int):

    i , j = 0, 0
    n, m = len(x), len(y)

    opt = [[0 for _ in range (m + 1)] for _ in range(n + 1)]
    for i in range(n + 1):
        opt[i][0] = i * delta
    
    for j in range(m + 1):
        opt[0][j] = j * delta

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if x[i - 1]  == y[j - 1]:
                opt[i][j] = opt[i - 1][j -1]
            
            else:
                opt[i][j] = min(opt[i - 1][j - 1] + alpha[x[i -1]+y[j-1]], opt[i - 1][j] + delta, opt[i][j - 1] + delta)

    
    l = m + n
    i = n
    j = m
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

        elif opt[i][j - 1]  + delta == opt[i][j]:
            y_final[y_position] = y[j - 1]
            x_final[x_position] = "_"
            x_position -= 1
            y_position -= 1
            j -= 1       
        
        elif opt[i-1][j] + delta == opt[i][j]:
            x_final[x_position] = x[i - 1]
            y_final[y_position] = "_"
            x_position -= 1
            y_position -= 1
            i -= 1
    
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
    del x_final
    del y_final

    return s1,s2,opt[n][m]



if __name__ == "__main__":
    tracemalloc.start()
    time_start=time.time()
    gene1, gene2 = generateString()
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
    delta = 30   
    s1,s2,opt = getMinimumPenalty(gene1, gene2, alpha, delta)
    _, peak = tracemalloc.get_traced_memory()    
    writeToFile(s1[0:50]+' '+s1[-50:]+'\n'+s2[0:50]+' '+s2[-50:]+'\n'+str(opt)+'\n'+str(float(time.time()-time_start))+'\n'+str(float(peak / 10**3))+'\n')
    tracemalloc.stop()