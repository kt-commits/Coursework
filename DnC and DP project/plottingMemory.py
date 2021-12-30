import numpy as np
import matplotlib.pyplot as plt
x=[1,2,4,8,16,32,64,128,256,512,1024,2048]
memory_efficient=[17.625,17.729,17.737,17.777,17.795,17.825,24.199,39.321,73.369,132.438,265.429,516.811]
memory_basic=[17.655,17.759,17.767,17.777,
17.795,36.179,134.991,516.935,1968.605,7722.941,31209.405,125663.447]

plt.title('Memory vs Problem Size')
plt.ylabel('Memory in KB')
plt.xlabel('Problem Size (length of string in characters)')
plt.plot(x, memory_efficient,label="Efficient")
plt.plot(x,memory_basic,label='Basic')
plt.legend()
plt.show()