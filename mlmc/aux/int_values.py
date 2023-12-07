from numpy import *

#from p4vasp.matrix import *
from math import *
import mymath

class Int_values:
  """from cartesians and given reciept (intwhat, intwhere) and
  lattice vectors (lattmat) calculates values of internal
  coordinates
  """

  def __init__(self,cart,coords,lattmat,ircdata):
  # cart must be an array!!!
    for i in range(len(coords)):
     if coords[i].dtyp=='simple':
      if coords[i].tag=='X':
        coords[i].value=self.set_singles(cart,0,coords[i].what[0],coords[i].where[0],lattmat)
      elif coords[i].tag=='Y':
        coords[i].value=self.set_singles(cart,1,coords[i].what[0],coords[i].where[0],lattmat)
      elif coords[i].tag=='Z':
        coords[i].value=self.set_singles(cart,2,coords[i].what[0],coords[i].where[0],lattmat)
      elif coords[i].tag=='Xr':
        coords[i].value=self.set_singles(cart,0,coords[i].what[1],coords[i].where[1],lattmat)
      elif coords[i].tag=='Yr':
        coords[i].value=self.set_singles(cart,1,coords[i].what[1],coords[i].where[1],lattmat)
      elif coords[i].tag=='Zr':
        coords[i].value=self.set_singles(cart,2,coords[i].what[1],coords[i].where[1],lattmat)
      elif coords[i].tag=='fX':
        coords[i].value=self.set_fsingles(cart,0,coords[i].what,lattmat)
      elif coords[i].tag=='fY':
        coords[i].value=self.set_fsingles(cart,1,coords[i].what,lattmat)
      elif coords[i].tag=='fZ':
        coords[i].value=self.set_fsingles(cart,2,coords[i].what,lattmat)
      elif coords[i].tag=='R':
        coords[i].value=self.set_lengths(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='IR1':
        coords[i].value=self.set_inverse_lengths1(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='IR6':
        coords[i].value=self.set_inverse_lengths6(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='M':
        coords[i].value=self.set_midlengths(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='RatioR':
        coords[i].value=self.set_ratior(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='A':
        coords[i].value=self.set_angles(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='AL':
        coords[i].value=self.set_directed_angles(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='acosZR':
        coords[i].value=self.set_acosXR(cart,coords[i].what,coords[i].where,lattmat,2)
      #elif coords[i].tag=='AL':
      #  coords[i].value=self.set_directed_angles(cart,coords[i].what,coords[i].where,lattmat,1)
      #elif coords[i].tag=='AL':
      #  coords[i].value=self.set_directed_angles(cart,coords[i].what,coords[i].where,lattmat,2)
      elif coords[i].tag=='T':
        coords[i].value=self.set_dihs(cart,coords[i].what,coords[i].where,lattmat)
      elif coords[i].tag=='tV':
        coords[i].value=self.set_tetrahedralVol(cart,coords[i].what,lattmat)
      elif coords[i].tag=='LR':
       coords[i].value=self.set_llength(coords[i].what,lattmat)
      elif coords[i].tag=='LA':
       coords[i].value=self.set_langle(coords[i].what,lattmat)
      elif coords[i].tag=='LB':
       coords[i].value=self.set_lbngle(coords[i].what,lattmat)
      elif coords[i].tag=='LV':
       coords[i].value=self.set_lvolume(lattmat)
     elif coords[i].dtyp=='sum':
       coords[i].value=self.set_sum(cart,coords[i],lattmat)
     elif coords[i].dtyp=='tsum':
       coords[i].value=self.set_tsum(cart,coords[i],lattmat)
     elif coords[i].dtyp=='norm':
       coords[i].value=self.set_norm(cart,coords[i],lattmat)
     elif coords[i].dtyp=='cn':
       coords[i].value=self.set_cnum(cart,coords[i],lattmat)
     elif coords[i].dtyp=='is':
       coords[i].value=self.set_is(cart,coords[i],lattmat,ircdata)

  def shortest_dist(self,cartesians,lattmat,atom1,atom2):
    """finds the shortest distance between two atoms
    """
    cart1=cartesians[atom1]
    cart2=cartesians[atom2]
    dists=[]
    what=[]

    for i in [-1,0,1]:
      for j in [-1,0,1]:
        for k in [-1,0,1]:
          trans=i*lattmat[0]+j*lattmat[1]+k*lattmat[2]
          point2=cart2+trans
          dist=(sum((cart1-point2)**2))**0.5
          dists.append(dist)
          what.append([i,j,k])

    dists=array(dists)
    dummy=argmin(dists)
    return [[0,0,0],[what[dummy][0],what[dummy][1],what[dummy][2]]]
	
  def set_lengths(self,cart,what,where,lattmat):
    """Calculates bond lengths
    """
    index1=what[0]
    index2=what[1]
    #where=array(self.shortest_dist(cart,lattmat,what[0],what[1]))
    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    bond=(sum((a-b)**2))**0.5
    return bond

  def set_inverse_lengths1(self,cart,what,where,lattmat):
    """Calculates inverse bond lengths
    """
    index1=what[0]
    index2=what[1]
    #where=array(self.shortest_dist(cart,lattmat,what[0],what[1]))
    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    bond=(sum((a-b)**2))**0.5
    bond=1./bond
    return bond

  def set_inverse_lengths6(self,cart,what,where,lattmat):
    """Calculates inverse bond lengths powered by 6
    """
    index1=what[0]
    index2=what[1]
    #where=array(self.shortest_dist(cart,lattmat,what[0],what[1]))
    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    bond=(sum((a-b)**2))**0.5
    bond=1./bond**6
    return bond



  def set_acosXR(self,cart,what,where,lattmat,mu):
    """acos(X/R); X,R- xth component and length of vector
    """
    index1=what[0]
    index2=what[1]
    #where=array(self.shortest_dist(cart,lattmat,what[0],what[1]))
    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    rvect=a-b
    r=(sum((a-b)**2))**0.5
    #xr=abs(rvect[mu]/r)
    xr=(rvect[mu]/r)
    xr=acos(xr)
    return xr


  def set_midlengths(self,cart,what,where,lattmat):
    """Calculates distance from an atom to midpoint between two atoms.
    """
    index1=what[0]
    index2=what[1]
    index3=what[2]

    #where=array(self.shortest_dist(cart,lattmat,what[0],what[1]))
    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    #where2=array(self.shortest_dist(cart,lattmat,what[1],what[2]))
    #c=cart[index3]+dot(where[1],lattmat)+matrixmultiply(where2[1],lattmat)
    c=cart[index3]+dot(where[2],lattmat)
    m=(b+c)/2
    bond=(sum((a-m)**2))**0.5
    return bond

  def set_ratior(self,cart,what,where,lattmat):
    r1=self.set_lengths(cart,what[:2],where[:2],lattmat)
    r2=self.set_lengths(cart,what[2:],where[2:],lattmat)
    return r1/r2

  def set_angles(self,cart,what,where,lattmat):
    """Calculates angles.
    """
    index1=what[0]
    index2=what[1]
    index3=what[2]
    #delement1=self.shortest_dist(cart,lattmat,what[1],what[0])
    #delement2=self.shortest_dist(cart,lattmat,what[1],what[2])
    #where=array([delement1[1],delement1[0],delement2[1]])

    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    c=cart[index3]+dot(where[2],lattmat)
    vector1=a-b
    vector2=c-b
    size1=mymath.vector_size(vector1)
    size2=mymath.vector_size(vector2)
    angle=sum(vector1*vector2)/(size1*size2)
    if angle>1:angle=1
    elif angle<-1:angle=-1
    angle=abs(acos(angle))
    return angle

  def set_directed_angles(self,cart,what,where,lattmat):
    """Calculates directed angles 
    """
    index1=what[0]
    index2=what[1]
    index3=what[2]

    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    c=cart[index3]+dot(where[2],lattmat)
    d=where[3][0]*lattmat[0]+where[3][1]*lattmat[1]+where[3][2]*lattmat[2]
    vector1=a-b
    vector2=c-b
    vector3=c-d
    size1=mymath.vector_size(vector1)
    size2=mymath.vector_size(vector2)
    angle=sum(vector1*vector2)/(size1*size2)
    if angle>1:angle=1
    elif angle<-1:angle=-1
    angle=abs(acos(angle))

    cross1=mymath.cross(vector1,vector2)
    #cross2=mymath.cross(vector2,vector3)
    #cross1_size=mymath.vector_size(cross1)
    #cross2_size=mymath.vector_size(cross2)
    #fu=sum(cross1*cross2)/(cross1_size*cross2_size)
    #if fu>1: fu=1.0
    #if fu<-1: fu=-1.0
    #dangle=acos(fu)
    if sum(cross1*vector3)>=0: angle=2*pi-angle
    return angle

  def set_dihs(self,cart,what,where,lattmat):
    """Calculates torsions.
    """
    index1=what[0]
    index2=what[1]
    index3=what[2]
    index4=what[3]
    #delement1=self.shortest_dist(cart,lattmat,what[1],what[0])
    #delement2=self.shortest_dist(cart,lattmat,what[1],what[2])
    #delement3=self.shortest_dist(cart,lattmat,what[2],what[3])
    #where=array([delement1[1],delement1[0],delement2[1],array(delement3[1])+array(delement2[1])])

    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    c=cart[index3]+dot(where[2],lattmat)
    d=cart[index4]+dot(where[3],lattmat)
    vector1=a-b
    vector2=b-c
    vector3=c-d
    cross1=mymath.cross(vector1,vector2)
    cross2=mymath.cross(vector2,vector3)
    cross1_size=mymath.vector_size(cross1)
    cross2_size=mymath.vector_size(cross2)
    fuck=sum(cross1*cross2)/(cross1_size*cross2_size)
    if fuck>1: fuck=1.0
    if fuck<-1: fuck=-1.0
    dangle=acos(fuck)
    if sum(cross1*vector3)>=0: dangle=-dangle
    return dangle
 
  def set_tetrahedralVol(self,cart,what,lattmat):
    """Calculates torsions.
    """
    index1=what[0]
    index2=what[1]
    index3=what[2]
    index4=what[3]

    delement1=self.shortest_dist(cart,lattmat,what[0],what[1])
    delement2=self.shortest_dist(cart,lattmat,what[0],what[2])
    delement3=self.shortest_dist(cart,lattmat,what[0],what[3])
    where=[delement1[0],delement1[1],delement2[1],delement3[1]]

    a=cart[index1]+dot(where[0],lattmat)
    b=cart[index2]+dot(where[1],lattmat)
    c=cart[index3]+dot(where[2],lattmat)
    d=cart[index4]+dot(where[3],lattmat)

    vector1=b-a
    vector2=c-a
    vector3=d-a
    #print 'vector1',vector1
    #print 'vector2',vector2
    #print 'vector3',vector3
    cross1=mymath.cross(vector1,vector2)
    #tv=sum(cross1*vector3)/2
    #volume must be positive
    tv=abs(sum(cross1*vector3)/6)
    return tv


  def set_singles(self,cart,xyz,what,where,lattmat):
    index=what
    #a=cart[index]+dot(where,lattmat)
    a=cart[index]
    internal=a[xyz]
    #internal=cart[index][xyz]
    return internal

  def set_fsingles(self,cart,xyz,what,lattmat):
    index=what[0]
    x=cart[index]
    x=dot(x,linalg.inv(lattmat))
    internal=x[xyz]
    return internal

  def set_tsum(self,cart,coord,lattmat):
    complexcoord=0.0
    dummyp=0.
    dummyq=0.
    dummyr=0.
    dummys=0.
    dummyt=0.
    for i in range(len(coord.tag)):
      if abs(coord.coefs[i])>0.:
        if coord.tag[i]=='T':
          dist=self.set_dihs(cart,coord.what[i],coord.where[i],lattmat)
          if coord.coefs[i]>0.:
            dummyp=mymath.min_image_cyclic(dummyp,dist)
            dummyq+=coord.coefs[i]*dist
            dummyp=dist
            dummys+=coord.coefs[i]
          elif coord.coefs[i]<0.:
            dummyp=mymath.min_image_cyclic(dummyp,dist)
            dummyr+=coord.coefs[i]*dist
            dummyp=dist
            dummyt+=coord.coefs[i]
    if dummys>0.:
      dummyq=dummyq/dummys
      dummyq=mymath.min_image_cyclic(0.,dummyq)
    if dummyt<0.:
      dummyr=dummyr/dummyt
      dummyq=mymath.min_image_cyclic(dummyq,dummyr)

    complexcoord=dummyq+dummyr
    return complexcoord



  def set_sum(self,cart,coord,lattmat):
    complexcoord=0.0
    for i in range(len(coord.tag)):
      if abs(coord.coefs[i])>0.:
        if coord.tag[i]=='X':
          dist=self.set_singles(cart,0,coord.what[i][0],coord.where[i][0],lattmat)
        if coord.tag[i]=='Y':
          dist=self.set_singles(cart,1,coord.what[i][0],coord.where[i][0],lattmat)
        if coord.tag[i]=='Z':
          dist=self.set_singles(cart,2,coord.what[i][0],coord.where[i][0],lattmat)
        if coord.tag[i]=='Xr':
          dist=self.set_singles(cart,0,coord.what[i][1],coord.where[i][1],lattmat)
        if coord.tag[i]=='Yr':
          dist=self.set_singles(cart,1,coord.what[i][1],coord.where[i][1],lattmat)
        if coord.tag[i]=='Zr':
          dist=self.set_singles(cart,2,coord.what[i][1],coord.where[i][1],lattmat)
        if coord.tag[i]=='fX':
          dist=self.set_fsingles(cart,0,coord.what[i],lattmat)
        if coord.tag[i]=='fY':
          dist=self.set_fsingles(cart,1,coord.what[i],lattmat)
        if coord.tag[i]=='fZ':
          dist=self.set_fsingles(cart,2,coord.what[i],lattmat)
        if coord.tag[i]=='R':
          dist=self.set_lengths(cart,coord.what[i],coord.where[i],lattmat)
        if coord.tag[i]=='IR1':
          dist=self.set_inverse_lengths1(cart,coord.what[i],coord.where[i],lattmat)
        if coord.tag[i]=='IR6':
          dist=self.set_inverse_lengths6(cart,coord.what[i],coord.where[i],lattmat)
        if coord.tag[i]=='A':
          dist=self.set_angles(cart,coord.what[i],coord.where[i],lattmat)
        if coord.tag[i]=='T':
          dist=self.set_dihs(cart,coord.what[i],coord.where[i],lattmat)
        if coord.tag[i]=='tV':
          dist=self.set_tetrahedralVol(cart,coord.what[i],lattmat)
        if coord.tag[i]=='RatioR':
          dist=self.set_ratior(cart,coord.what[i],coord.where[i],lattmat)
      
        complexcoord+=coord.coefs[i]*dist
    return complexcoord

  def set_norm(self,cart,coord,lattmat):
    complexcoord=0.0
    for i in range(len(coord.tag)):
      if coord.tag[i]=='X':
        dist=self.set_singles(cart,0,coord.what[i][0],coord.where[i][0],lattmat)
      if coord.tag[i]=='Y':
        dist=self.set_singles(cart,1,coord.what[i][0],coord.where[i][0],lattmat)
      if coord.tag[i]=='Z':
        dist=self.set_singles(cart,2,coord.what[i][0],coord.where[i][0],lattmat)
      if coord.tag[i]=='Xr':
        dist=self.set_singles(cart,0,coord.what[i][1],coord.where[i][1],lattmat)
      if coord.tag[i]=='Yr':
        dist=self.set_singles(cart,1,coord.what[i][1],coord.where[i][1],lattmat)
      if coord.tag[i]=='Zr':
        dist=self.set_singles(cart,2,coord.what[i][1],coord.where[i][1],lattmat)
      if coord.tag[i]=='fX':
        dist=self.set_fsingles(cart,0,coord.what[i],lattmat)
      if coord.tag[i]=='fY':
        dist=self.set_fsingles(cart,1,coord.what[i],lattmat)
      if coord.tag[i]=='fZ':
        dist=self.set_fsingles(cart,2,coord.what[i],lattmat)
      if coord.tag[i]=='R':
        dist=self.set_lengths(cart,coord.what[i],coord.where[i],lattmat)
      if coord.tag[i]=='A':
        dist=self.set_angles(cart,coord.what[i],coord.where[i],lattmat)
      if coord.tag[i]=='T':
        dist=self.set_dihs(cart,coord.what[i],coord.where[i],lattmat)
      if coord.tag[i]=='RatioR':
        dist=self.set_ratior(cart,coord.what[i],coord.where[i],lattmat)
      complexcoord+=(coord.coefs[i]*dist)**2
    complexcoord=complexcoord**0.5
    return complexcoord
    
  def set_cnum(self,cart,coord,lattmat):
    complexcoord=0.0
    for i in range(len(coord.tag)):
      if coord.tag[i]=='R':
        if abs(coord.coefs[i])>1e-4:
          #print coord.coefs[i]
          #dist=coord.value
          dist=self.set_lengths(cart,coord.what[i],coord.where[i],lattmat)
          dummyq=abs(dist/coord.coefs[i])
          if abs(dummyq-1.0)<1e-4: dummyq=1.0001
          if coord.coefs[i] > 0.:
            complexcoord=complexcoord+(1.0-dummyq**9.)/(1.0-dummyq**14.)
          else:
            complexcoord=complexcoord-(1.0-dummyq**9.)/(1.0-dummyq**14.)
          #print complexcoord,dist,coord.coefs[i],i,'i'
    return complexcoord

  def set_is(self,cart,coord,lattmat,ircdata):
    complexcoord=0.0
    xx=[]
    for i in range(len(coord.tag)):
      if coord.tag[i]=='X':
        dist=self.set_singles(cart,0,coord.what[i][0],coord.where[i][0],lattmat)
      if coord.tag[i]=='Y':
        dist=self.set_singles(cart,1,coord.what[i][0],coord.where[i][0],lattmat)
      if coord.tag[i]=='Z':
        dist=self.set_singles(cart,2,coord.what[i][0],coord.where[i][0],lattmat)
      if coord.tag[i]=='Xr':
        dist=self.set_singles(cart,0,coord.what[i][1],coord.where[i][1],lattmat)
      if coord.tag[i]=='Yr':
        dist=self.set_singles(cart,1,coord.what[i][1],coord.where[i][1],lattmat)
      if coord.tag[i]=='Zr':
        dist=self.set_singles(cart,2,coord.what[i][1],coord.where[i][1],lattmat)
      if coord.tag[i]=='fX':
        dist=self.set_fsingles(cart,0,coord.what[i],lattmat)
      if coord.tag[i]=='fY':
        dist=self.set_fsingles(cart,1,coord.what[i],lattmat)
      if coord.tag[i]=='fZ':
        dist=self.set_fsingles(cart,2,coord.what[i],lattmat)
      if coord.tag[i]=='R':
        dist=self.set_lengths(cart,coord.what[i],coord.where[i],lattmat)
      if coord.tag[i]=='A':
        dist=self.set_angles(cart,coord.what[i],coord.where[i],lattmat)
      if coord.tag[i]=='T':
        dist=self.set_dihs(cart,coord.what[i],coord.where[i],lattmat)
      if coord.tag[i]=='tV':
        dist=self.set_tetrahedralVol(cart,coord.what[i],lattmat)
      if coord.tag[i]=='RatioR':
        dist=self.set_ratior(cart,coord.what[i],coord.where[i],lattmat)
      xx.append(dist)

    dummyp=0.
    dummyq=0.
    for k in range(len(ircdata)):
      dummyr=0.
      for j in range(len(xx)):
        if coord.coefs[j]>0.:
          dummyr+=coord.coefs[j]*(xx[j]-ircdata[k,j])**2
      dummyr=e**(-dummyr)
      dummyp=dummyp+(k)*dummyr
      dummyq=dummyq+dummyr
    if dummyq>0.:
      isVal=dummyp/(dummyq*(len(ircdata)-1))
    return isVal


  def set_llength(self,intwhat,lattmat):
    index1=intwhat[0]
    llength=(sum(lattmat[index1]*lattmat[index1]))**0.5
    return llength

  def set_langle(self,intwhat,lattmat):
    a=intwhat[0]
    b=intwhat[1]
    diffav=lattmat[a]
    diffbv=lattmat[b]
    d1=mymath.vector_size(diffav)
    d2=mymath.vector_size(diffbv)
    cosalpha=(sum(diffav*diffbv)/(d1*d2))
    alpha=acos(cosalpha)
    return alpha

  def set_lbngle(self,intwhat,lattmat):
    a=intwhat[0]
    b=intwhat[1]
    c=intwhat[2]
    av=lattmat[a]
    bv=lattmat[b]
    cv=lattmat[c]

    v1=mymath.cross(av,bv)
    v2=mymath.cross(av,cv)

    v1=v1/sum(v1**2)**0.5
    v2=v2/sum(v2**2)**0.5
    alpha=sum(v1*v2)
    alpha=acos(alpha)
    return alpha

  def set_lvolume(self,lattmat):
    l1=lattmat[0]
    l2=lattmat[1]
    l3=lattmat[2]
    volume=mymath.cross(l1,l2)
    volume=sum(volume*l3)
    return volume

  def set_ratiollength(self,intwhat,lattmat):
    index1=intwhat[0]
    llengtha=(sum(lattmat[index1]*lattmat[index1]))**0.5
    index2=intwhat[1]
    llengthb=(sum(lattmat[index2]*lattmat[index2]))**0.5
    ratio=llengtha/llengthb
    return ratio

    
    
    
    
    
    
    
    
    
    
    






