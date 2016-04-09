class EuropeanCallOption(object):   
   def __init__(self):
      self.__rho = 0 
      self.__sigma = 0 
      self.__r = 0
      self.__q = 0
      self.__theta = 0
      self.__U = 0
      pass
   def __del__(self):
      print "Destructor started"

   # Set rho, sigma, r, q, theta
   def SetParameter(self,rho, sigma, r, q, theta):
      self.__rho = rho
      print self.__rho
      print "SetParameter"
   
   # Set boundary condition
   def SetBC(self):
      print "SetBC"
   
   # Set num of grid in s, v, t direction
   def SetNumofGrid(self):
      print "SetNumofGrid"
   
   # Construct matrix, rhs 
   def Build(self):
      print "Build"
   
   # Solve the linear system
   def Solve(self):
      print "Solve"

