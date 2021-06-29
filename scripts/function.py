import numpy as np



def norm(p:np.array)->float:
    '''
    Compute the L2-norm of a vector
    Parameters
    ---
    param : type
        introduction
    Return
    ---
    return instruction
    '''
    sum = 0
    for i in p:
        sum += i*i
    
    return np.sqrt(sum)

def angle(p1, p2, p3):
    p21 = p1 - p2
    p23 = p3 - p2
    return np.real(np.arccos(np.dot(p21, p23) / norm(p21) / norm(p23)))