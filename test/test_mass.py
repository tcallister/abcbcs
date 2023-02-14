import unittest
from abcbcs import mass

class TestMassConversion(unittest.TestCase):

    # Set up by preparing dictionary of reference values
    def setUp(self):

        m1 = 10.
        m2 = 3.
        q = m2/m1
        eta = q/(1.+q)**2
        mtot = m1+m2
        mc = eta**(3./5.)*mtot
        self.reference_dictionary = {\
            'mass_1':m1,
            'mass_2':m2,
            'mass_ratio':q,
            'symmetric_mass_ratio':eta,
            'total_mass':mtot,
            'chirp_mass':mc}
        
    def tearDown(self):
        del self.reference_dictionary

    def test_conversion(self):

        # Loop across all combinations of parameters
        for i,keyA in enumerate(list(self.reference_dictionary.keys())[:-1]):
            for keyB in list(self.reference_dictionary.keys())[i+1:]:

                # The one combination we *can't* do is q and eta, so skip this one
                if 'mass_ratio' not in (keyA,keyB) and 'symmetric_mass_ratio' not in (keyA,keyB):

                    # Pass combinations of parameters to mass translator
                    kwargs = {keyA:self.reference_dictionary[keyA],keyB:self.reference_dictionary[keyB]}
                    masses = mass.masses(**kwargs)

                    # Test derived params against precomputed values
                    self.assertAlmostEqual(masses.mass_1,self.reference_dictionary['mass_1'])
                    self.assertAlmostEqual(masses.mass_2,self.reference_dictionary['mass_2'])
                    self.assertAlmostEqual(masses.mass_ratio,self.reference_dictionary['mass_ratio'])
                    self.assertAlmostEqual(masses.symmetric_mass_ratio,self.reference_dictionary['symmetric_mass_ratio'])
                    self.assertAlmostEqual(masses.total_mass,self.reference_dictionary['total_mass'])
                    self.assertAlmostEqual(masses.chirp_mass,self.reference_dictionary['chirp_mass'])
