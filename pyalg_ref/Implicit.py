import numpy as np

#### Setting 
X_LEN = 1.0
Y_LEN = 1.0
Cell_I = 10
Cell_J = 10
deltaT = 0.01
#### Initial Condition
rhoc = 1.0
k = 0.1
####
I = Cell_I+1
J = Cell_J+1
deltaX = X_LEN / Cell_I
deltaY = Y_LEN / Cell_J
####

def initialGrid():
    tn = np.zeros([I, J])
    tn = SetBoundary(tn)
    return tn

def PringGrid(T):
    with np.printoptions(precision=3, suppress=True):
        print(np.flipud(T.T))

def indexA2T(i, j):
    """
    Index from A_25_25 to T_5_5
    """
    pass

def indexT2A(i, j):
    """
    Index from T_5_5 to A_25_25
    """
    index = i*(J)+ j
    return index
    

def initialA():
    """
    A would scale from T00->T0J->T10->T1J->...->TI0->TIJ
    """
    size = I * J
    s = k * deltaT/rhoc

    A = np.zeros([size, size])
    A = setABoundary(A)

    ## Deal with inner voxcel
    for i in range(1, I-1):
        for j in range(1, J-1):
            ind_ij = indexT2A(i, j)
            ind_r = indexT2A(i+1, j)
            ind_l = indexT2A(i-1, j)
            ind_t = indexT2A(i, j+1)
            ind_b = indexT2A(i, j-1)

            A[ind_ij, ind_ij] = 2*s/deltaX**2 +2*s/deltaY**2 + 1

            A[ind_ij, ind_r] = -s/deltaX**2
            A[ind_ij, ind_l] = -s/deltaX**2

            A[ind_ij, ind_t] = -s/deltaY**2
            A[ind_ij, ind_b] = -s/deltaY**2
    
    return A

def initialb():
    size = I * J
    b = np.zeros([size, 1])
    b = setbBoundary(b)
    return b

def setbBoundary(b):
    for j in range(J):
        ind = indexT2A(0, j)
        b[ind, 0] = 100

        ind2 = indexT2A(I-1, j)
        b[ind2, 0] = 0

    for i in range(1, I):
        ind3 = indexT2A(i, 0)
        b[ind3, 0] = 0

        ind4 = indexT2A(i, J-1)
        b[ind4, 0] = 0
    
    return b

def setABoundary(A):
    for j in range(T.shape[1]):
        ind = indexT2A(0, j)
        A[ind, ind] = 1

        ind2 = indexT2A(I-1, j)
        A[ind2, ind2] = 1

    for i in range(1, T.shape[0]):
        ind3 = indexT2A(i, 0)
        A[ind3, ind3] = 1

        ind4 = indexT2A(i, J-1)
        A[ind4, ind4] = 1
    
    return A

def SetBoundary(T):
    for j in range(T.shape[1]):
        T[0, j] = 100
        T[I-1, j] = 0

    for i in range(1, T.shape[0]):
        T[i, 0] = 0
        T[i, J-1] = 0
    
    return T

def print_related(i, j):
    print("i, j: ({}, {}) -> {}\n".format(i, j, indexT2A(i,j)))
    print("i+1, j -> ({})\n".format(indexT2A(i+1, j)))
    print("i-1, j -> ({})\n".format(indexT2A(i-1, j)))
    print("i, j-1 -> ({})\n".format(indexT2A(i, j-1)))
    print("i, j+1 -> ({})\n".format(indexT2A(i, j+1)))

def plot_meshgrid(Z):
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    x = np.linspace(0, 1, I)
    y = np.linspace(0, 1, J)
    X, Y = np.meshgrid(x, y)
    ax.plot_wireframe(X, Y, Z, color='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show()

if __name__ == "__main__":

    T = initialGrid()
    A = initialA()
    b = initialb()

    timestep = 0
    while (timestep <= 20):
        x = np.linalg.solve(A, b)
        # print("\n\nTimestep: {}".format(timestep))
        # PringGrid(x.reshape(I, J))
        b = setbBoundary(x)
        timestep += 1
    
    re = np.flipud(x.reshape(I, J).T)
    plot_meshgrid(re)

