import unittest
from rdkit import Chem

class ElectrostaticPotential(unittest.TestCase):

    ''' The tests cover all functional groups, except nitro alkane.'''

    def setUp(self):
        
        from electrostatic_potential import ElectrostaticPotential
        self.ep = ElectrostaticPotential()

    def test_standardisation(self):
    
        smiles = 'C(C)C(=O)[O-]'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles)
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'carboxylic_acid': 1, 'alcohol': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'carboxylic_acid': 1})

        
        smiles = 'c1cccc[nH+]1'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles)
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'pyridine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'pyridine': 1})


        smiles = 'c1ccccc1[NH3+]'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles) 
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'amine': 1, 'benzene': 1, 'aniline': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'aniline': 1})
        

        smiles = 'c1ccccc1[O-]'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles) 
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'benzene': 1, 'phenol': 1, 'alcohol': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'phenol': 1})


        smiles = 'C[N-]S(=O)(=O)C'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles) 
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'sulfonamide': 1, 'amine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'sulfonamide': 1})


        smiles = 'C[N-]C=C'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles) 
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'alkene': 1, 'amine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'alkene': 1, 'amine': 1})


        smiles = 'c1ccc[n-]1'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles) 
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'pyrrole': 1, 'amine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'pyrrole': 1})


        smiles = 'CC[N-]C(=O)CC'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        neutral_smiles = self.ep.neutral_smiles()
        self.assertNotEqual(smiles, neutral_smiles) 
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'amide': 1, 'amine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'amide': 1})

    def test_adjusted_matches_partA(self):
    
        smiles = 'C1=CC=CC(=C1)[S](C2=CC=CC=C2)(C3=CC=CC(=C3)C4=CC(=CC=C4)[S](C5=CC=CC=C5)(C6=CC(=CC=C6)O)C7=CC(=C(C=C7)O)O)C8=CC=CC=C8'
        mol = Chem.MolFromSmiles(smiles)
        
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches,{'thioether': 12, 'benzene': 8, 'phenol': 3, 'alcohol': 3})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches,{'thioether': 2, 'benzene': 6, 'phenol': 2})

        self.assertAlmostEqual(ep_values.MaxAlpha, 3.8)
        self.assertAlmostEqual(ep_values.MaxBeta,3.6)
        self.assertAlmostEqual(ep_values.MinAlpha,1.0)
        self.assertAlmostEqual(ep_values.MinBeta,2.2)
        self.assertAlmostEqual(ep_values.TotalAlpha,15.6)
        self.assertAlmostEqual(ep_values.TotalBeta,25.8)
        self.assertAlmostEqual(ep_values.AverageAlpha,1.6)
        self.assertAlmostEqual(ep_values.AverageBeta,2.6)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised,0.021)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised,0.035)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised ,0.05)
        self.assertAlmostEqual(ep_values.BetaVSANormalised ,0.082)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised,1.186)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised,1.962)
       
        self.assertAlmostEqual(ep_values.Alpha_VSA0,0)
        self.assertAlmostEqual(ep_values.Alpha_VSA1,0.115)
        self.assertAlmostEqual(ep_values.Alpha_VSA2,0)
        self.assertAlmostEqual(ep_values.Alpha_VSA3,0)
        self.assertAlmostEqual(ep_values.Alpha_VSA4,0.032)

        self.assertAlmostEqual(ep_values.Beta_VSA0,0)
        self.assertAlmostEqual(ep_values.Beta_VSA1,0.176)
        self.assertAlmostEqual(ep_values.Beta_VSA2,0)
        self.assertAlmostEqual(ep_values.Beta_VSA3,0)
        self.assertAlmostEqual(ep_values.Beta_VSA4,0)    
       
    def test_groups_0(self):
    
        smiles = 'CC(=C)C(=O)OC(c1ccccc1)(c1ccccc1)c1ccncc1'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'alkene': 1, 'ether': 1, 'ester': 1, 'benzene': 2, 'pyridine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'alkene': 1, 'ester': 1, 'benzene': 2, 'pyridine': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 1.5)
        self.assertAlmostEqual(ep_values.MaxBeta, 7.0)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.7)
        self.assertAlmostEqual(ep_values.MinBeta, 1.1)
        self.assertAlmostEqual(ep_values.TotalAlpha, 5.6)
        self.assertAlmostEqual(ep_values.TotalBeta, 17.8)
        self.assertAlmostEqual(ep_values.AverageAlpha, 1.100)
        self.assertAlmostEqual(ep_values.AverageBeta, 3.600)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.017)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.054)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.038)
        self.assertAlmostEqual(ep_values.BetaVSANormalised , 0.121)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 1.246)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 3.962)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.292)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.066)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.038)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.076)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.041)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.084)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.000)    
    
    def test_groups_1(self):

        smiles = 'C(C(CCCCC(=O)O)C(NC)=S)=O'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'alkane': 2, 'aldehyde': 1, 'thioamide': 1, 'carboxylic_acid': 1, 'amine': 1, 'alcohol': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'alkane': 1, 'aldehyde': 1, 'thioamide': 1, 'carboxylic_acid': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 3.6)
        self.assertAlmostEqual(ep_values.MaxBeta, 5.8)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.4)
        self.assertAlmostEqual(ep_values.MinBeta, 0.3)
        self.assertAlmostEqual(ep_values.TotalAlpha, 8.9)
        self.assertAlmostEqual(ep_values.TotalBeta, 16.1)
        self.assertAlmostEqual(ep_values.AverageAlpha, 2.200)
        self.assertAlmostEqual(ep_values.AverageBeta, 4.000)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.041)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.074)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.100)
        self.assertAlmostEqual(ep_values.BetaVSANormalised , 0.182)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 8.960)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 16.209)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.073)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.054)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.118)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.073)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.196)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.000)      
 
    def test_groups_2(self):

        smiles = 'C(=O)(OCOCC1=CC=C(C=C1)S)N'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'carbamate': 1, 'ether': 2, 'thiophenol': 1, 'amine': 1, 'benzene': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'carbamate': 1, 'ether': 1, 'thiophenol': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 2.8)
        self.assertAlmostEqual(ep_values.MaxBeta, 7.3)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.9)
        self.assertAlmostEqual(ep_values.MinBeta, 2.2)
        self.assertAlmostEqual(ep_values.TotalAlpha, 5.5)
        self.assertAlmostEqual(ep_values.TotalBeta, 14.8)
        self.assertAlmostEqual(ep_values.AverageAlpha, 1.800)
        self.assertAlmostEqual(ep_values.AverageBeta, 4.900)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.026)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.069)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.063)
        self.assertAlmostEqual(ep_values.BetaVSANormalised , 0.169)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 3.561)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 9.581)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.055)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.145)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.055)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.056)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.078)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.070)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.000)
        
    def test_groups_3(self):

        smiles = 'C(CC(=O)CSC)(=O)N'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'ketone': 1, 'amide': 1, 'thioether': 1, 'amine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'ketone': 1, 'amide': 1, 'thioether': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 2.9)
        self.assertAlmostEqual(ep_values.MaxBeta, 8.3)
        self.assertAlmostEqual(ep_values.MinAlpha, 1.0)
        self.assertAlmostEqual(ep_values.MinBeta, 3.6)
        self.assertAlmostEqual(ep_values.TotalAlpha, 5.4)
        self.assertAlmostEqual(ep_values.TotalBeta, 17.7)
        self.assertAlmostEqual(ep_values.AverageAlpha, 1.800)
        self.assertAlmostEqual(ep_values.AverageBeta, 5.900)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.037)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.120)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.093)
        self.assertAlmostEqual(ep_values.BetaVSANormalised , 0.305)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 0.000)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 0.000)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.206)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.084)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.100)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.101)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.101)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.103)

    def test_groups_4(self):

        smiles = 'C(C(C(C1=CC(=CC=C1)Cl)C2=CC=C(C=C2)Br)C3=CC=C(C=C3)I)CC(S)C4=CC=CC(=C4)F'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)

        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'thiol': 1, 'aryl_chloride': 1, 'aryl_fluoride': 1, 'aryl_bromide': 1, 'aryl_iodide': 1, 'benzene': 4})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'thiol': 1, 'aryl_chloride': 1, 'aryl_fluoride': 1, 'aryl_bromide': 1, 'aryl_iodide': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 1.7)
        self.assertAlmostEqual(ep_values.MaxBeta, 2.7)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 1.6)
        self.assertAlmostEqual(ep_values.TotalAlpha, 4.4)
        self.assertAlmostEqual(ep_values.TotalBeta, 9.1)
        self.assertAlmostEqual(ep_values.AverageAlpha, 0.900)
        self.assertAlmostEqual(ep_values.AverageBeta, 1.800)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.007)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.014)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.019)
        self.assertAlmostEqual(ep_values.BetaVSANormalised , 0.039)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 0.431)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 0.891)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.052)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.122)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.015)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.022)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.000)
    
    def test_groups_5(self):

        smiles = 'C(C(C(Cl)Br)I)(C(C#C)F)N'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'alkyl_iodide': 1, 'alkyl_fluoride': 1, 'alkyl_chloride': 1, 'alkyl_bromide': 1, 'alkyne': 1, 'amine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'alkyl_iodide': 1, 'alkyl_fluoride': 1, 'alkyl_chloride': 1, 'alkyl_bromide': 1, 'alkyne': 1, 'amine': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 1.9)
        self.assertAlmostEqual(ep_values.MaxBeta, 7.8)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 2.2)
        self.assertAlmostEqual(ep_values.TotalAlpha, 7.2)
        self.assertAlmostEqual(ep_values.TotalBeta, 19.8)
        self.assertAlmostEqual(ep_values.AverageAlpha, 1.200)
        self.assertAlmostEqual(ep_values.AverageBeta, 3.300)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.020)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.056)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.078)
        self.assertAlmostEqual(ep_values.BetaVSANormalised , 0.214)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 3.515)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 9.666)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.243)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.349)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.110)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.062)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.000)

    def test_groups_6(self):

        smiles = 'C1(=CC=C(C=C1)C(C#N)NC(N)=O)O'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'amine': 2, 'nitrile': 1, 'benzene': 1, 'urea': 1, 'phenol': 1, 'alcohol': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'nitrile': 1, 'urea': 1, 'phenol': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 3.8)
        self.assertAlmostEqual(ep_values.MaxBeta, 8.3)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 2.7)
        self.assertAlmostEqual(ep_values.TotalAlpha, 6.8)
        self.assertAlmostEqual(ep_values.TotalBeta, 15.7)
        self.assertAlmostEqual(ep_values.AverageAlpha, 2.300)
        self.assertAlmostEqual(ep_values.AverageBeta, 5.200)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.036)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.082)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.084)
        self.assertAlmostEqual(ep_values.BetaVSANormalised, 0.194)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 10.877)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 25.113)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.066)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.072)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.064)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.072)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.076)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.075)

    def test_groups_7(self):

        smiles = 'C(NCC(C(C1=CC=C(C=C1)O)O)C#N)(N)=O'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)
        
        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'amine': 2, 'nitrile': 1, 'benzene': 1, 'urea': 1, 'phenol': 1, 'alcohol': 2})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'nitrile': 1, 'urea': 1, 'phenol': 1, 'alcohol': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 3.8)
        self.assertAlmostEqual(ep_values.MaxBeta, 8.3)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 2.7)
        self.assertAlmostEqual(ep_values.TotalAlpha, 9.5)
        self.assertAlmostEqual(ep_values.TotalBeta, 21.5)
        self.assertAlmostEqual(ep_values.AverageAlpha, 2.400)
        self.assertAlmostEqual(ep_values.AverageBeta, 5.400)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.040)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.091)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.097)
        self.assertAlmostEqual(ep_values.BetaVSANormalised, 0.219)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 40.654)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 92.006)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.054)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.111)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.052)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.059)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.125)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.062)        

    def test_groups_8(self):

        smiles = 'C[S](NCC(C1=COC=C1)=NCSSCNC(C)=N)(=O)=O'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)    

        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'sulfonamide': 1, 'furan': 1, 'amine': 2, 'imine': 2, 'disulfide': 1, 'amidine': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'sulfonamide': 1, 'furan': 1, 'imine': 1, 'disulfide': 1, 'amidine': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 2.8)
        self.assertAlmostEqual(ep_values.MaxBeta, 8.9)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 2.2)
        self.assertAlmostEqual(ep_values.TotalAlpha, 4.5)
        self.assertAlmostEqual(ep_values.TotalBeta, 26.4)
        self.assertAlmostEqual(ep_values.AverageAlpha, 0.900)
        self.assertAlmostEqual(ep_values.AverageBeta, 5.300)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.013)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.075)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.034)
        self.assertAlmostEqual(ep_values.BetaVSANormalised, 0.200)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 2.998)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 17.587)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.040)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.033)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.064)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.092)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.043)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.044)        

    def test_groups_9(self):

        smiles = 'C1=CC(=C[N]1[H])CC(CCCC2=CC=C(C=C2)N)N[S](C)=O'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)

        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'alkane': 1, 'pyrrole': 1, 'sulfinamide': 1, 'amine': 3, 'benzene': 1, 'aniline': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'alkane': 1, 'pyrrole': 1, 'sulfinamide': 1, 'aniline': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 3.2)
        self.assertAlmostEqual(ep_values.MaxBeta, 8.3)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.4)
        self.assertAlmostEqual(ep_values.MinBeta, 0.3)
        self.assertAlmostEqual(ep_values.TotalAlpha, 8.7)
        self.assertAlmostEqual(ep_values.TotalBeta, 18.0)
        self.assertAlmostEqual(ep_values.AverageAlpha, 2.200)
        self.assertAlmostEqual(ep_values.AverageBeta, 4.500)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.028)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.059)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.068)
        self.assertAlmostEqual(ep_values.BetaVSANormalised, 0.140)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 3.604)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 7.457)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.050)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.084)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.033)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.050)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.093)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.086)      

    def test_groups_10(self):

        smiles = 'C(C(C[S](C)(=O)=O)SC#N)(CC1=C[N](C=N1)[H])N=C=S'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)

        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'thioether': 2, 'amine': 1, 'sulfone': 1, 'nitrile': 1, 'isothiocyanate': 1, 'thiocyanate': 1, 'imidazole': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'sulfone': 1, 'isothiocyanate': 1, 'thiocyanate': 1, 'imidazole': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 3.7)
        self.assertAlmostEqual(ep_values.MaxBeta, 6.3)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 0.0)
        self.assertAlmostEqual(ep_values.TotalAlpha, 3.7)
        self.assertAlmostEqual(ep_values.TotalBeta, 11.7)
        self.assertAlmostEqual(ep_values.AverageAlpha, 0.900)
        self.assertAlmostEqual(ep_values.AverageBeta, 2.900)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.012)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.037)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.031)
        self.assertAlmostEqual(ep_values.BetaVSANormalised, 0.097)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 3.520)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 11.131)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.044)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.041)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.052)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.045)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.081)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.000)     

    def test_groups_11(self):

        smiles = 'C([P](OC)(OC)=O)C(C(C[S](=O)(OC)OC)C[P](C)(C)=O)C[S](C)=O'
        mol = Chem.MolFromSmiles(smiles)
        ep_values = self.ep.calc_electrostatic_potential_descriptors(mol)

        all_matches = self.ep.all_matches()
        self.assertEqual(all_matches, {'thioether': 1, 'phosphinate_diester': 1, 'sulfate_diester': 1, 'sulfoxide': 1, 'phosphine_oxide': 1})
        adjusted_matches = self.ep.adjusted_matches()
        self.assertEqual(adjusted_matches, {'phosphinate_diester': 1, 'sulfate_diester': 1, 'sulfoxide': 1, 'phosphine_oxide': 1})

        self.assertAlmostEqual(ep_values.MaxAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MaxBeta, 9.9)
        self.assertAlmostEqual(ep_values.MinAlpha, 0.0)
        self.assertAlmostEqual(ep_values.MinBeta, 4.7)
        self.assertAlmostEqual(ep_values.TotalAlpha, 0.0)
        self.assertAlmostEqual(ep_values.TotalBeta, 32.4)
        self.assertAlmostEqual(ep_values.AverageAlpha, 0.000)
        self.assertAlmostEqual(ep_values.AverageBeta, 8.100)
        self.assertAlmostEqual(ep_values.AlphaMolWtNormalised, 0.000)
        self.assertAlmostEqual(ep_values.BetaMolWtNormalised, 0.073)
        self.assertAlmostEqual(ep_values.AlphaVSANormalised , 0.000)
        self.assertAlmostEqual(ep_values.BetaVSANormalised, 0.210)
        self.assertAlmostEqual(ep_values.AlphaLogPNormalised, 0.000)
        self.assertAlmostEqual(ep_values.BetaLogPNormalised, 15.333)

        self.assertAlmostEqual(ep_values.Alpha_VSA0, 0.030)
        self.assertAlmostEqual(ep_values.Alpha_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA2, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Alpha_VSA4, 0.000)

        self.assertAlmostEqual(ep_values.Beta_VSA0, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA1, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA2, 0.068)
        self.assertAlmostEqual(ep_values.Beta_VSA3, 0.000)
        self.assertAlmostEqual(ep_values.Beta_VSA4, 0.083)        
        
if __name__ == '__main__':
    unittest.main()
    
