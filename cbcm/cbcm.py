def _q_from_m1_m2(m1,m2):
    return m2/m1

def _eta_from_q(q):
    return q/(1.+q)**2

def _mtot_from_m1_m2(m1,m2):
    return m1+m2

def _mc_from_mtot_eta(mtot,eta):
    return eta**(3./5.)*mtot

def _eta_from_mtot_mc(mtot,mc):
    return (mc/mtot)**(5./3.)

def _m1_from_mtot_q(mtot,q):
    return mtot/(1.+q)

def _mtot_from_mc_eta(mc,eta):
    return mc/eta**(3./5.)

def _m2_from_m1_eta(m1,eta):
    return m1*(1.-2.*eta - (1.-4.*eta)**0.5)/(2.*eta)

def _m1_from_m2_eta(m2,eta):
    return m2*(1.-2.*eta + (1.-4.*eta)**0.5)/(2.*eta)

def _q_from_eta(eta):
    return (1.-2.*eta - (1-4*eta)**0.5)/(2.*eta)

def _q_from_mc_m1(mc,m1):
    
    # Follow Cardano's formula, e.g. https://en.wikipedia.org/wiki/Cubic_equation
    r = (mc/m1)**5.
    delta = r**2/4. - r**3/27.
    u1 = r/2. + delta**0.5
    u2 = r/2. - delta**0.5

    return (u1**(1./3.) + u2**(1./3.)).real

def _q_from_mc_m2(mc,m2):

    # Follow Cardano's formula for depressed polynomial to get q_tilde = q + 1/3
    c1 = -1./3.
    c0 = (2.-27.*(m2/mc)**5)/27.
    delta = c0**2/4. + c1**3./27.
    u1 = -c0/2. + delta**0.5
    u2 = -c0/2. - delta**0.5
    q_tilde = u1**(1./3.) + u2**(1./3.)
    return (q_tilde - 1./3.).real

class cbc():

    """
    Class to help me save time transforming between mass parameters
    """

    def __init__(self,**kwargs):

        # We need exactly two variables to describe a set of cbc masses
        if len(kwargs)!=2:
            raise Exception("Requires exactly two arguments")

        # Verify that we have known parameters
        dimensionful_params = ["total_mass","chirp_mass","mass_1","mass_2"]
        dimensionless_params = ["mass_ratio","symmetric_mass_ratio"] 
        if any(kw not in dimensionful_params+dimensionless_params for kw in kwargs):
            raise Exception("Parameters must be among the following:",dimensionful_params,dimensionless_params)

        # Make sure at least one parameter is unitful
        if not any(param in kwargs for param in dimensionful_params):
            raise Exception("Requires at least one dimensionful parameter",dimensionful_params)

        # Assign attributes
        self.__dict__.update(kwargs)
        self._check_params()
        self._populate()

    def _check_params(self):

        all_params = ["total_mass","chirp_mass",'mass_1',"mass_2","mass_ratio","symmetric_mass_ratio"]
        for param in all_params:
            if hasattr(self,param):
                if self.__dict__[param]<0:
                    raise Exception("All mass parameters must be positive")

        if hasattr(self,'mass_ratio'):
            if self.mass_ratio>1:
                raise Exception("Mass ratio must be less than unity")

        if hasattr(self,'symmetric_mass_ratio'):
            if self.symmetric_mass_ratio>0.25:
                raise Exception("Symmetric mass ratio must be less than 0.25")

        if hasattr(self,'mass_1'):
            if hasattr(self,'mass_2'):
                if self.mass_2>self.mass_1:
                    raise Exception("Mass 2 must be less than or equal to mass 1")
            elif hasattr(self,'total_mass'):
                if self.total_mass<self.mass_1:
                    raise Exception("Total mass must be greater than mass 1")
            elif hasattr(self,'chirp_mass'):
                if self.chirp_mass>self.mass_1/2**(1./5.):
                    raise Exception("Unphysical chirp mass; chirp mass cannot exceed m1/2**(1/5)")

        elif hasattr(self,'mass_2'):
            if hasattr(self,'total_mass'):
                if self.total_mass<self.mass_2:
                    raise Exception("Total mass must be greater than mass 2")
            elif hasattr(self,'chirp_mass'):
                if self.chirp_mass>self.mass_2/2**(1./5.):
                    raise Exception("Unphysical chirp mass; chirp mass cannot exceed m2/2**(1/5)")

        elif hasattr(self,'total_mass'):
            if hasattr(self,'chirp_mass'):
                if self.chirp_mass>self.total_mass/2**(6./5.):
                    raise Exception("Unphysical chirp mass; chirp mass cannot exceed total_mass/2**(6/5)")

    def _populate(self):

        """
        Method to fully populate class attributes with self-consistent set of mass parameters.
        """

        # Consider various cases of input source parameters
        if hasattr(self,'mass_1'):

            if hasattr(self,'mass_2'):
                self.__dict__.update(mass_ratio = _q_from_m1_m2(self.mass_1,self.mass_2))
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(total_mass = _mtot_from_m1_m2(self.mass_1,self.mass_2))
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

            elif hasattr(self,'total_mass'):
                self.__dict__.update(mass_2 = self.total_mass - self.mass_1)
                self.__dict__.update(mass_ratio = _q_from_m1_m2(self.mass_1,self.mass_2))
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

            elif hasattr(self,'chirp_mass'):
                self.__dict__.update(mass_ratio = _q_from_mc_m1(self.chirp_mass,self.mass_1))
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)
                self.__dict__.update(total_mass = _mtot_from_m1_m2(self.mass_1,self.mass_2))

            elif hasattr(self,'mass_ratio'):
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(total_mass = _mtot_from_m1_m2(self.mass_1,self.mass_2))
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

            elif hasattr(self,'symmetric_mass_ratio'):
                self.__dict__.update(mass_2 = _m2_from_m1_eta(self.mass_1,self.symmetric_mass_ratio))
                self.__dict__.update(total_mass = self.mass_1 + self.mass_2)
                self.__dict__.update(mass_ratio = self.mass_2/self.mass_1)
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

        elif hasattr(self,'mass_2'):

            if hasattr(self,'total_mass'):
                self.__dict__.update(mass_1 = self.total_mass - self.mass_2)
                self.__dict__.update(mass_ratio = _q_from_m1_m2(self.mass_1,self.mass_2))
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

            elif hasattr(self,'chirp_mass'):
                self.__dict__.update(mass_ratio = _q_from_mc_m2(self.chirp_mass,self.mass_2))
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(mass_1 = self.mass_2/self.mass_ratio)
                self.__dict__.update(total_mass = _mtot_from_m1_m2(self.mass_1,self.mass_2))

            elif hasattr(self,'mass_ratio'):
                self.__dict__.update(mass_1 = self.mass_2/self.mass_ratio)
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(total_mass = _mtot_from_m1_m2(self.mass_1,self.mass_2))
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

            elif hasattr(self,'symmetric_mass_ratio'):
                self.__dict__.update(mass_1 = _m1_from_m2_eta(self.mass_2,self.symmetric_mass_ratio))
                self.__dict__.update(mass_ratio = self.mass_2/self.mass_1)
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))
                self.__dict__.update(total_mass = self.mass_1 + self.mass_2)

        elif hasattr(self,'total_mass'):

            if hasattr(self,'chirp_mass'):
                self.__dict__.update(symmetric_mass_ratio = _eta_from_mtot_mc(self.total_mass,self.chirp_mass))
                self.__dict__.update(mass_ratio = _q_from_eta(self.symmetric_mass_ratio))
                self.__dict__.update(mass_1 = _m1_from_mtot_q(self.total_mass,self.mass_ratio))
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)

            elif hasattr(self,'mass_ratio'):
                self.__dict__.update(mass_1 = _m1_from_mtot_q(self.total_mass,self.mass_ratio))
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

            elif hasattr(self,'symmetric_mass_ratio'):
                self.__dict__.update(mass_ratio = _q_from_eta(self.symmetric_mass_ratio))
                self.__dict__.update(mass_1 = _m1_from_mtot_q(self.total_mass,self.mass_ratio))
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)
                self.__dict__.update(chirp_mass = _mc_from_mtot_eta(self.total_mass,self.symmetric_mass_ratio))

        elif hasattr(self,'chirp_mass'):
            
            if hasattr(self,'mass_ratio'):
                self.__dict__.update(symmetric_mass_ratio = _eta_from_q(self.mass_ratio))
                self.__dict__.update(total_mass = _mtot_from_mc_eta(self.chirp_mass,self.symmetric_mass_ratio))
                self.__dict__.update(mass_1 = _m1_from_mtot_q(self.total_mass,self.mass_ratio))
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)

            elif hasattr(self,'symmetric_mass_ratio'):
                self.__dict__.update(mass_ratio = _q_from_eta(self.symmetric_mass_ratio))
                self.__dict__.update(total_mass = _mtot_from_mc_eta(self.chirp_mass,self.symmetric_mass_ratio))
                self.__dict__.update(mass_1 = _m1_from_mtot_q(self.total_mass,self.mass_ratio))
                self.__dict__.update(mass_2 = self.mass_ratio*self.mass_1)
