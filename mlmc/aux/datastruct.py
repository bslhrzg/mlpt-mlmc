import string

class Complextype:
 def __init__(self,dtyp,coefs,tag,what,whattags,where,value,status):
   self.dtyp=dtyp           # sum/norm
   self.coefs=coefs
   self.tag=tag
   self.what=what
   self.whattags=whattags
   self.where=where
   self.value=value
   self.status=status
   self.labels=None

class Linesearch:
 def __init__(self):
   self.a=[None]         # braketing left
   self.b=[None]         # minimum
   self.c=[None]         # braketing right
   self.direction=None
   self.energies=[None]
   self.steps=[None]
   self.oldgrad=None                 # gradients
   self.newgrad=None
   self.oldpos=None                  # positions
   self.nexpos=None
   self.status='u'
   self.cg=None

class Elementprop:
  def __init__(self):
    self.covalentradii={'Ru': 1.5, 'Re': 1.55, 'Ra': 2.100, 'Rb': 1.670, 'Rn': 0.200, 'Rh': 1.650,\
     'Be': 0.550, 'Ba': 1.540, 'Bi': 1.740, 'Bk': 0.200, 'Br': 1.410, 'H': 0.430, 'P': 1.250, \
     'Os': 1.570, 'Hg': 1.900, 'Ge': 1.370, 'Gd': 1.990, 'Ga': 1.420, 'Pr': 2.020, 'Pt': 1.700, \
     'Pu': 0.200, 'C': 0.900, 'Pb': 1.740, 'Pa': 1.810, 'Pd': 1.700, 'Cd': 1.890, 'Po': 1.880, \
     'Pm': 2.000, 'Ho': 1.940, 'Hf': 1.770, 'K': 1.530, 'He': 0.741, 'Mg': 1.300, 'Mo': 1.670,\
     'Mn': 1.550, 'O': 0.880, 'S': 1.220, 'W': 1.570, 'Zn': 1.650, 'Eu': 2.190, 'Zr': 1.760, \
     'Er': 1.930, 'Ni': 1.700, 'Na': 1.170, 'Nb': 1.680, 'Nd': 2.010, 'Ne': 0.815, 'Np': 1.750, \
     'Fr': 0.200, 'Fe': 1.540, 'B': 1.030, 'F': 0.840, 'Sr': 1.320, 'N': 0.880, 'Kr': 1.069, \
     'Si': 1.000, 'Sn': 1.660, 'Sm': 2.000, 'V': 1.530, 'Sc': 1.640, 'Sb': 1.660, 'Se': 1.420,\
     'Co': 1.530, 'Cm': 0.200, 'Cl': 1.190, 'Ca': 1.190, 'Cf': 1.730, 'Ce': 2.030, 'Xe': 1.750, \
     'Lu': 1.920, 'Cs': 1.870, 'Cr': 1.550, 'Cu': 1.000, 'La': 2.070, 'Li': 0.880, 'Tl': 1.750, \
     'Tm': 1.920, 'Th': 1.990, 'Ti': 1.670, 'Te': 1.670, 'Tb': 1.960, 'Tc': 1.550, 'Ta': 1.630, \
     'Yb': 2.140, 'Dy': 1.950, 'I': 1.600, 'U': 1.780, 'Y': 1.980, 'Ac': 2.080, 'Ag': 1.790, \
     'Ir': 1.520, 'Am': 1.710, 'Al': 1.550, 'As': 1.410, 'Ar': 0.995, 'Au': 1.500, 'At': 0.200, \
     'In': 1.830}
    self.mendelejev={'Ru': 4, 'Re': 5, 'Ra': 6, 'Rb': 4, 'Rn': 5, 'Rh': 4, 'Be': 1, 'Ba': 5, \
    'Bi': 5, 'Bk': 6, 'Br': 3, 'H': 0, 'P': 2, 'Os': 5, 'Hg': 5, 'Ge': 3, 'Gd': 5, 'Ga': 3, \
    'Pr': 5, 'Pt': 5, 'Pu': 6, 'C': 1, 'Pb': 5, 'Pa': 6, 'Pd': 4, 'Cd': 4, 'Po': 5, 'Pm': 5, \
    'Ho': 5, 'Hf': 5, 'K': 3, 'He': 0, 'Mg': 2, 'Mo': 4, 'Mn': 3, 'O': 1, 'S': 2, 'W': 5, \
    'Zn': 3, 'Eu': 5, 'Zr': 4, 'Er': 5, 'Ni': 3, 'Na': 2, 'Nb': 4, 'Nd': 5, 'Ne': 1, 'Np': 6, \
    'Fr': 6, 'Fe': 3, 'B': 1, 'F': 1, 'Sr': 4, 'N': 1, 'Kr': 3, 'Si': 2, 'Sn': 4, 'Sm': 5, \
    'V': 3, 'Sc': 3, 'Sb': 4, 'Se': 3, 'Co': 3, 'Cm': 6, 'Cl': 2, 'Ca': 3, 'Cf': 6, 'Ce': 5, \
    'Xe': 4, 'Lu': 5, 'Cs': 5, 'Cr': 3, 'Cu': 3, 'La': 5, 'Li': 1, 'Tl': 5, 'Tm': 5, 'Th': 6, \
    'Ti': 3, 'Te': 4, 'Tb': 5, 'Tc': 4, 'Ta': 5, 'Yb': 5, 'Dy': 5, 'I': 4, 'U': 6, 'Y': 4, \
    'Ac': 6, 'Ag': 4, 'Ir': 5, 'Am': 6, 'Al': 2, 'As': 3, 'Ar': 2, 'Au': 5, 'At': 5, 'In': 4}
