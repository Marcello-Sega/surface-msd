import numpy as np
import MDAnalysis as mda
import pytim
from pytim.observables import Observable, Position, Mass

class ResidueWeightedObservable(Observable):
    
    def __init__(self, obs, weights, normalized, wrap):
        self.obs, self.weights, self.normalized, self.wrap = obs, weights, normalized, wrap

    def compute(self, inp=None, kargs=None):
        obs,weights,normalized, wrap = self.obs, self.weights, self.normalized, self.wrap
        vals = obs.compute(inp)
        if weights is None :
            w = np.ones(vals.shape[0])
        else:
            w = weights.compute(inp)
        box = inp.universe.dimensions[obs.dirmask]
        labels, gi, first_idx = np.unique(inp.resindices, return_inverse=True, return_index=True)
        G, d = labels.size, inp.n_atoms
        if wrap:
            # reconstruct starting from the 1st atom in the residiue
            # this works also if the set of residues is non homogeneous
            dv = vals - vals[first_idx]
            # minimal-image shift to bring each atom near its reference
            # round(dX/L) is the integer number of boxes separating atom from ref
            newvals  = vals  - box * np.floor(dv/box + 0.5) # +0.5 for stability
        else:
            newvals = vals
                
        try:
            num = np.zeros((G, newvals.shape[1]), dtype=type(newvals))
            np.add.at(num, first_idx, w[:,None] * newvals)
        except:
            num = np.zeros(G, dtype=type(newvals))
            np.add.at(num, first_idx, w * newvals)
            
        if normalized:
            den = np.zeros(G, dtype=type(w))
            np.add.at(den, first_idx, w)
            try: 
                return num / den[:, None]
            except:
                return num/den

        else: 
            return num

class ResidueCenterOfMassPosition(ResidueWeightedObservable):
    def __init__(self,arg='xyz'):
        self._obs = Position(arg)
        self._weights = Mass()
        self._normalized,self._wrap = True, True
        ResidueWeightedObservable.__init__(self, self._obs, self._weights, self._normalized, self._wrap)
        

class ResidueCenterOfMassVelocity(ResidueWeightedObservable):
    def __init__(self,arg='xyz'):
        self._obs = Velocity(arg)
        self._weights = Mass()
        self._normalized,self._wrap = True, False
        ResidueWeightedObservable.__init__(self, self._obs, self._weights, self._normalized, self._wrap)
        

class LayerMSD(object):
    def __init__(self, refgroup, direction='xy',npoints=100):
        self.data = np.zeros(npoints+1)
        self.norm = np.zeros(npoints+1)
        self.refgroup = refgroup
        self.npoints = npoints
        self.ind0 = refgroup.residues.resindices[0]
        self.com = ResidueCenterOfMassPosition(direction)
        self.all_coms = []
           
    def sample(self,inp):
        s1,s2 = len(self.refgroup),len(self.com.obs.dirmask)
        a = np.empty((s1,s2))
        a[:] = np.nan
        self.a=a  
        resind = inp.residues.resindices
        a[resind - self.ind0 ] = self.com.compute(inp)
        self.all_coms.append(a)
        nan_mask = np.zeros(a.shape) * a
        self.norm[0] += len(resind)
        for dt in range(1, min(len(self.all_coms),self.npoints+1)):
            end = len(self.all_coms) - 1
            nan_mask *= self.all_coms[end-dt]
            dr = (self.all_coms[end] - self.all_coms[end-dt]) + nan_mask
            dr = np.abs(dr)
            dr[np.greater(dr,inp.universe.dimensions[0]/2.)]-=inp.universe.dimensions[0] # PBC
            self.data[dt] += np.nansum(dr**2)
            self.norm[dt] += np.count_nonzero(~np.isnan(dr))/len(dr[0])
            self.dr=dr
            if np.count_nonzero(np.isnan(nan_mask))/s2 == s1:
                break

        if len(self.all_coms)>self.npoints: self.all_coms.pop(0)

