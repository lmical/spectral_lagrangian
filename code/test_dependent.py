import numpy as np

class DATA_CLASS:
    def __init__(self,test):
        self.test = test
        if test=="Sod": #Sod
            # Extrema
            self.xL=0
            self.xR=1
            # Final time
            self.T=0.231
            # Periodicity of the mesh
            self.periodic=False
        else:
            print("Error in class DATA_CLASS, in test_dependent, test not available")
            quit()