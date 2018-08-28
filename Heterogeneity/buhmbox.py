import numpy as np
from numpy.linalg import inv
from scipy.stats import norm, chi2
import sys

class BuhmBox:
    def bb(self, cases, controls):
        self.latest = {}
        self.index = 0
        self.numer = 0.0
        self.denom = 0.0
        self.num_snps = cases.shape[1]
        self.z = np.zeros(self.num_snps)

    def get_full_mat(self):
        return self.R

    def get_values(self, clist):
        numer = 0.0
        denom = 0.0
        num_snps = len(clist)
        for i in range(num_snps):
            for j in range(i+1,num_snps):
                wij = self.z[clist[i]] * self.z[clist[j]]
                yij = self.Y[clist[i],clist[j]]
                numer += wij * yij
                denom += wij * wij
        SBB = numer / np.sqrt(denom)
        return SBB

    def get_weights(self):
        return self.z

class BuhmBox_GWAS(BuhmBox):
    def bb(self, cases, controls):
        super().bb(cases, controls)
        N = float(len(cases))
        Np = float(len(controls))
        self.R = np.corrcoef(cases.T)
        self.Rp = np.corrcoef(controls.T)
        self.Y = np.sqrt(N*Np/(N+Np)) * (self.R-self.Rp)
        for i in range(self.num_snps):
            p_case = np.sum(cases[:,i]) / (2*cases.shape[0])
            p_control = np.sum(controls[:,i]) / (2*controls.shape[0])
            gamma = p_case/(1-p_case) / (p_control/(1-p_control))
            self.z[i] = np.sqrt(p_control*(1-p_control)) * (gamma-1) / (p_control*(gamma-1) + 1)

class BuhmBox_TWAS(BuhmBox):
    def bb(self, cases, controls):
        super().bb(cases, controls)
        N = float(len(cases))
        Np = float(len(controls))
        self.R = np.corrcoef(cases.T)
        self.Rp = np.corrcoef(controls.T)
        self.Y = np.sqrt(N*Np/(N+Np)) * (self.R-self.Rp)
        n = controls.shape[0]
        for i in range(self.num_snps):
            mui_control = np.mean(controls[:,i])
            mui_case = np.mean(cases[:,i])
            stdi_control = np.std(controls[:,i])
            self.z[i] = (mui_case - mui_control)/stdi_control

