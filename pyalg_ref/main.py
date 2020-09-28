import numpy as np

#### Setting 
X_LEN = 1.0
Y_LEN = 1.0
Cell_I = 4
Cell_J = 4
deltaT = 0.01
#### Initial Condition
rhoc = 1.0
k = 0.1
####
I = Cell_I+1
J = Cell_J+1
deltaX = X_LEN / I
deltaY = Y_LEN / J
####

def initialGrid():
    tn = np.zeros([I, J])
    tn = SetBoundary(tn)
    return tn

def advect_ij(T, i, j):
    assert (i != 0) | (j != 0) | (i!= I-1) | (j != J-1), "i/j should not be 0 or LEN_X/Y"
    partial2_x2 = (T[i+1, j] + T[i-1, j] - 2 * T[i, j]) / deltaX**2
    partial2_y2 = (T[i, j-1] + T[i, j+1] - 2 * T[i, j]) / deltaY**2

    Tn_ij = ( k * (partial2_x2 + partial2_y2) * deltaT + T[i,j] ) / rhoc

    return Tn_ij


def advect(T):
    Tn = initialGrid()
    for i in range(1, I-1):
        for j in range(1, J-1):
            Tn[i, j] = advect_ij(T, i, j)
    
    return Tn

def SetBoundary(T):
    ## T[0, j] = 100
    for j in range(T.shape[1]):
        T[0, j] = 100
        T[I-1, j] = 0

    for i in range(1, T.shape[0]):
        T[i, 0] = 0
        T[i, J-1] = 0
    
    return T

def PringGrid(T):
    with np.printoptions(precision=3, suppress=True):
        print(np.flipud(T.T))

def initialA(T):
    size = I * J
    A = np.zeros([size, size])

if __name__ == "__main__":
    T = initialGrid()
    PringGrid(T)
    timestep = 0
    while (timestep <= 10):
        T = advect(T)
        print("\n\nTimestep: {}".format(timestep))
        PringGrid(T)
        timestep += 1


