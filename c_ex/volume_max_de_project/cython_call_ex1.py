from volume_cython_ex1 import DifferentialEvolution

#fittness function (cost function)
def evaluate(designVariablel):
    surface = 80.0
    # if pType is 1, the penality is negative (maximization problem)
    # if pType is 0, the penality is positive (minimization problem)
    penality = -1000
  
    z = (surface-designVariablel[0]*designVariablel[1])/(2.0*(designVariablel[0]\
         +designVariablel[1]))
    volume = designVariablel[0]*designVariablel[1]*z
    
    if(volume <= 0):
        return penality
  
    # box length and width need to be larger than 0
    if(designVariablel[0] <= 0):
        return penality
  
    if(designVariablel[1] <= 0):
        return penality
    return volume
    
volume = DifferentialEvolution(evaluate, 1, 3, 2, 100, 0.6, 0.85, [0, 0], [50, 50], 100, 10)
volume.run()