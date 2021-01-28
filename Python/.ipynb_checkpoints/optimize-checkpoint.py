import scipy
import numpy as np
import matplotlib
from matplotlib import pyplot

numgenes = 10
numtfs = 5

alpha = np.array([1.0 for i in range(numgenes)])
beta = np.array([0.0 for i in range(numgenes)])
cs = np.random.rand(numgenes, numtfs)
tfa = np.random.rand(numgenes, numtfs)

"""
Calculates the predicted gene expression values for the given parameters, assuming enhancers
"""
def getGeneExpression(alpha, beta, cs, tfa):
    return [beta[i] + alpha[i]*sum([tfa[i][x] / (tfa[i][x] + cs[i][x]) for x in range(len(cs[0]))]) for i in range(len(cs))]

exprs = getGeneExpression(alpha, beta, cs, tfa)
# print(exprs)

"""
Calculates the error in the predicted gene expression values given the true values and the parameters
"""
def getErrorFromParameters(exprs, alpha, beta, cs, tfa):
    return sum(map(lambda x: (x[1] - x[0])**2, zip(getGeneExpression(alpha, beta, cs, tfa), exprs)))

"""
Calculates error given predicted and actual values
"""
def getError(exprs, pred):
    return sum(map(lambda x: (x[1] - x[0])**2, zip(pred, exprs)))

cs2 = np.random.rand(numgenes, numtfs)
tfa2 = np.random.rand(numgenes, numtfs)

# print(getErrorFromParameters(exprs, alpha, beta, cs2, tfa2))


"""
Calculates the predicted gene expression values for the given parameters, assuming enhancers
In addition, this simulates transcription factor at index tfdel is deleted, therefore setting
the activity of that transcription factor in the tfa matrix to zero
"""
def getGeneExpressionWithDeletion(alpha: np.array, beta: np.array, cs: np.array, tfa: np.array, tfdel: int):
    updatedtfa = tfa.copy()
    newcol = np.zeros(len(tfa))
    updatedtfa[:,tfdel] = newcol
    return getGeneExpression(alpha, beta, cs, updatedtfa)

# for i in range(numtfs):
#     print(getGeneExpressionWithDeletion(alpha, beta, cs, tfa, i))


"""
Calculates total error from matrix of gene expressions where column i is the expression vector for
the system with transcription factor i deleted. Thus this is a numgenes * numtfs matrix
"""
def getErrorWithDeletionsFromParameters(expr_matrix, alpha, beta, cs, tfa):
    pred_matrix = np.array([getGeneExpressionWithDeletion(alpha, beta, cs, tfa, i) for i in range(len(tfa[0]))]).transpose()
    return np.square(pred_matrix - expr_matrix).sum()

"""
Calculates total error from matrix of gene expressions where column i is the expression vector for
the system with transcription factor i deleted. Thus this is a numgenes * numtfs matrix
"""
def getErrorWithDeletions(expr_matrix, pred_matrix):
    return np.square(pred_matrix - expr_matrix).sum()


del_exprs = np.random.rand(numgenes, numtfs)
print(getErrorWithDeletionsFromParameters(del_exprs, alpha, beta, cs, tfa))