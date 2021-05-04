# This is a library of 24 electrostatic potential based molecular descriptors 
# based on Hunter, C. A. (2004). Angewandte Chemie - International Edition, 43(40), 5310–5324
# https://doi.org/10.1002/anie.200301739
# This code was created by Andre Frade (with support of Patrick McCabe and Richard Cooper)

try:
    import os
    import pickle
    import numpy as np
    import pandas as pd
    from collections import namedtuple
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.ML.Descriptors import MoleculeDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError as e:
            print('Error {}'.format(e))

class ElectrostaticPotential:

    def __init__(self):
        """
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        
        ep = ElectrostaticPotential()
        ep_values = ep.calc_electrostatic_potential_descriptors(mol)
        
        ep_values is a named tuple containing: MaxAlpha MaxBeta MinAlpha MinBeta TotalAlpha TotalBeta 
        AverageAlpha AverageBeta AlphaMolWtNormalised BetaMolWtNormalised AlphaVSANormalised BetaVSANormalised
        AlphaLogPNormalised BetaLogPNormalised Alpha_VSA0 Alpha_VSA1 Alpha_VSA2 Alpha_VSA3 Alpha_VSA4 Beta_VSA0
        Beta_VSA1 Beta_VSA2 Beta_VSA3 Beta_VSA4')
    
        :param ht_filename: absolute filename of csv file containing alpha, beta group values
        Alpha and beta group values taken from https://doi.org/10.1002/anie.200301739
        :param complexgroup_filename: absolute filename of metadata dict
        A complex group is a supergroup that contains subgroups. For example phenol contains a benzene and an alcohol
        """
        
        
        import pkg_resources
        ht_filepath = pkg_resources.resource_filename('electrostatic_potential', 'data/ABvalues.csv')
        cg_filepath = pkg_resources.resource_filename('electrostatic_potential', 'data/complex_groups')
        
        with open(complexgroup_filepath, "rb") as complexgroup_file:
            self._complex_groups = pickle.load(cg_filepath)
            
        self._ht_df = pd.read_csv(ht_filepath)    
        self._smarts_patterns = self._ht_df.smarts.values
        self._neutral_smiles = ''
        self._all_matches = {}
        self._all_indices = {}
        self._adjusted_matches = {}
        self._adjusted_indices = {}

    def _uncharge_molecule(self, mol):
    
        '''Neutralizes molecular charges for functional group identification.
        It may fail.'''
    
        smiles = Chem.MolToSmiles(mol)

        if '+' in smiles or '-' in smiles:
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)           
            
        self._neutral_smiles = Chem.MolToSmiles(mol)
        return mol

    def _find_group(self, smarts, mol):

        indices = []

        try:
            patt = Chem.MolFromSmarts(smarts)
            indices = [list(x) for x in mol.GetSubstructMatches(patt)]
        except RuntimeError as e:
            print('Error {} {}'.format(smarts, e))

        return indices

    def _find_matches(self, mol):

        self._all_matches = {}
        self._all_indices = {} 
        
        for group in self._smarts_patterns:
            indices = self._find_group(group, mol)

            if len(indices) > 0:
                group = self._ht_df.group.loc[self._ht_df.smarts == group].values[0]
                self._all_matches[group] = len(indices)
                self._all_indices[group] = indices
               
    def _adjust_matches(self):
    
        adjusted_indices = dict(self._all_indices)
        
        # [1] ADJUST GROUP OCCURENCE COUNT - check if no atom belongs to more than one group occurence
        # =======================================================================
        
        #for i in adjusted_indices: print(i, adjusted_indices[i])
        for group in adjusted_indices:
            group_events = adjusted_indices[group]
            remove_ =[]
            for i in range(len(group_events)-1):
                if set(group_events[i]) & set(group_events[i+1]):
                    remove_.append(group_events[i])
            for r in remove_:
                adjusted_indices[group].remove(r)
   
        # [2] ADJUST COMPLEX GROUP COUNT - remove subgroup counts
        # =======================================================================
        
        # get indices of cg and their potential sg identified in molecule
        # -----------------------------------------------------------------------
        for cg in self._complex_groups: # for each cg found in mol
        
            if not cg in adjusted_indices:
                continue
           
            for cg_indices in adjusted_indices[cg]: # get its indices
                temp = [] # accumulates sg indices overlaping each cg
                    
                for sg in self._complex_groups[cg]: # and for the corresponding subgroups, get their indices too
                
                    if not sg in adjusted_indices:
                        continue
                           
                    # check whether any of the matches of a group is a sg (search by indices overlap)
                    # -----------------------------------------------------------------------                
                    for sg_indices in adjusted_indices[sg]:
                    
                        if set(cg_indices) & set(sg_indices): # check if sg and cg indices overlap

                            if sg_indices not in temp: # if so, store the sg indices once
                                temp.append(sg_indices)
                
                    # keep only cg and simpler ones that dont act as sg
                    # -----------------------------------------------------------------------    
                    if sg in adjusted_indices: # keep simpler groups that dont act as sg
                        keep = [x for x in adjusted_indices[sg] if x not in temp]
                        adjusted_indices[sg] = keep

        self._adjusted_indices = { k:v for k, v in adjusted_indices.items() if len(v) > 0}
        self._adjusted_matches = { k:len(v) for k, v in adjusted_indices.items() if len(v) > 0}
        
    def _max_min_alpha_beta(self):

        max_alpha = 0
        max_beta = 0
        min_alpha = 0
        min_beta = 0

        groups = [group for group in self._adjusted_matches]
        alphas = self._ht_df.alpha.loc[self._ht_df.group.isin(groups)].values
        betas = self._ht_df.beta.loc[self._ht_df.group.isin(groups)].values

        if len(alphas) > 0:
            max_alpha = max(alphas)
            min_alpha = min(alphas)

        if len(betas) > 0:
            max_beta = max(betas)
            min_beta = min(betas)

        return max_alpha, max_beta, min_alpha, min_beta
        
    def _average_total_alpha_beta(self):
    
        total_alpha = 0
        total_beta = 0
        average_alpha = 0
        average_beta = 0
        
        groups = []
        for group in self._adjusted_matches:
            for i in range (self._adjusted_matches[group]):    
                groups.append(group)
                
        alphas = [self._ht_df.alpha.loc[self._ht_df.group == group].values[0] for group in groups]
        betas = [self._ht_df.beta.loc[self._ht_df.group == group].values[0] for group in groups]
        
        if len(alphas) > 0:
            total_alpha = round(np.sum(alphas), 1)
            average_alpha = round(np.mean(alphas), 1)

        if len(betas) > 0:
            total_beta = round(np.sum(betas), 1)
            average_beta = round(np.mean(betas), 1)
            
        return total_alpha, total_beta, average_alpha, average_beta
     
    def _normalised_alpha_beta(self, mol, total_alpha, total_beta):
        
        alpha_mw_norm, beta_mw_norm = 0, 0
        alpha_vsa_norm, beta_vsa_norm = 0, 0
        alpha_logp_norm, beta_logp_norm = 0, 0
   
        # by MolWt
        '''The average molecular weight of the molecule'''
        obj = MoleculeDescriptors.MolecularDescriptorCalculator(['MolWt'])
        res = obj.CalcDescriptors(mol)
        mw = res[0]
        if mw > 0:
            alpha_mw_norm = round(total_alpha/mw,3)
            beta_mw_norm = round(total_beta/mw,3)
        
        # by VSA
        '''Labute. P. J. Mol. Graph. Mod. _18_ 464-477 (2000)'''
        vsa = Chem.MolSurf.LabuteASA(mol)
        if vsa > 0:
            alpha_vsa_norm = round(total_alpha/vsa,3)
            beta_vsa_norm = round(total_beta/vsa,3)    
        
        # by LogP
        '''S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)'''
        logp = Chem.Crippen.MolLogP(mol)
        if logp > 0:        
            alpha_logp_norm = round(total_alpha/logp,3)
            beta_logp_norm = round(total_beta/logp,3)        
        
        return alpha_mw_norm, beta_mw_norm, alpha_vsa_norm, beta_vsa_norm, alpha_logp_norm, beta_logp_norm

    def _VSA_alpha_beta(self, mol):
    
        ''' calculates Labute's Approximate Surface Area.
        Definition from P. Labute's article in the Journal of the Chemical Computing Group
        and J. Mol. Graph. Mod.  _18_ 464-477 (2000)
        '''
        
        # calculate the per-atom contributions to the surface area
        (ats, hs)= Chem.rdMolDescriptors._CalcLabuteASAContribs(mol, includeHs = True)
        atomic_vsa = [round(vsa, 3) for vsa in ats]
        total_vsa = np.sum(atomic_vsa)
        
        # calculate the per-atom contributions to the Gasteiger partial charge
        # this dictates the most positive (donor) and negative (acceptor) atom of the group
        AllChem.ComputeGasteigerCharges(mol)
        atomic_pcharge = [round(mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge'), 3) for i in range(mol.GetNumAtoms())]

        res_alpha = {}
        res_beta = {}
        
        for group in self._adjusted_indices: 
        
            vsa_alpha = 0
            vsa_beta = 0
            
            # for each group event, get vsa of 'donnor' and 'acceptor' atoms
            for event in self._adjusted_indices[group]:
            
                # get p charge for elements in group
                local_pcharge = [atomic_pcharge[atom] for atom in event]
                local_vsa = [atomic_vsa[atom] for atom in event]
                
                # identify alpha (min p charge) and beta (max p charge) atoms, and retrieve their vsa contribution
                vsa_alpha += np.sum([local_vsa[i] for i in range(len(event)) if local_pcharge[i] == min(local_pcharge)])
                vsa_beta += np.sum([local_vsa[i] for i in range(len(event)) if local_pcharge[i] == max(local_pcharge)])
        
            # get group's alpha and beta value
            alpha = self._ht_df.alpha.loc[self._ht_df.group == group].values[0]
            beta = self._ht_df.beta.loc[self._ht_df.group == group].values[0]
            
            # store normalised vsa's for those alphas and betas
            res_alpha[alpha]= vsa_alpha/total_vsa
            res_beta[beta]= vsa_beta/total_vsa
        
        a_VSA0, a_VSA1, a_VSA2, a_VSA3, a_VSA4 = self._alpha_VSA(res_alpha)
        b_VSA0, b_VSA1, b_VSA2, b_VSA3, b_VSA4 = self._beta_VSA(res_beta)
        
        return a_VSA0, a_VSA1, a_VSA2, a_VSA3, a_VSA4, b_VSA0, b_VSA1, b_VSA2, b_VSA3, b_VSA4
        
    def _alpha_VSA(self, dic):
    
        ''' alpha: 5 bins, 1 unit interval [0], ]0, 1] ]1, 2] ]2, 3], ]3, inf['''

        alpha_VSA0 = 0
        alpha_VSA1 = 0
        alpha_VSA2 = 0
        alpha_VSA3 = 0
        alpha_VSA4 = 0

        for alpha in dic:

            if alpha == 0:
                alpha_VSA0 += dic[alpha]

            elif 0 < alpha <=1:
                alpha_VSA1 += dic[alpha]

            elif 1 < alpha <=2:
                alpha_VSA2 += dic[alpha]

            elif 2 < alpha <=3:
                alpha_VSA3 += dic[alpha]

            elif 3 < alpha:
                alpha_VSA4 += dic[alpha]

        return round(alpha_VSA0,3), round(alpha_VSA1,3), round(alpha_VSA2,3), round(alpha_VSA3,3), round(alpha_VSA4,3)
    
    def _beta_VSA(self, dic):
    
        ''' beta: 5 bins, 2 units interval [0, 2], ]2, 4] ]4, 6], ]6, 8], ]8, inf['''
    
        beta_VSA0 = 0
        beta_VSA1 = 0
        beta_VSA2 = 0
        beta_VSA3 = 0
        beta_VSA4 = 0

        for beta in dic:

            if 0 <= beta <=2:
                beta_VSA0 += dic[beta]

            elif 2 < beta <=4:
                beta_VSA1 += dic[beta]

            elif 4 < beta <=6:
                beta_VSA2 += dic[beta]

            elif 6 < beta <=8:
                beta_VSA3 += dic[beta]

            elif 8 < beta:
                beta_VSA4 += dic[beta]

        return round(beta_VSA0,3), round(beta_VSA1,3), round(beta_VSA2,3), round(beta_VSA3,3), round(beta_VSA4,3)
        
    def neutral_smiles(self):
        
        return self._neutral_smiles        

    def adjusted_matches(self):
        
        return self._adjusted_matches
   
    def all_matches(self):

        return self._all_matches
          
    def calc_electrostatic_potential_descriptors(self, mol):
    
        '''
        :param mol: a mol object
        :return ep_values: a namedtuple containing 
                           MaxAlpha MaxBeta MinAlpha MinBeta 
                           TotalAlpha TotalBeta AverageAlpha AverageBeta 
                           AlphaMolWtNormalised BetaMolWtNormalised 
                           AlphaVSANormalised BetaVSANormalised 
                           AlphaLogPNormalised BetaLogPNormalised 
                           Alpha_VSA0 Alpha_VSA1 Alpha_VSA2 Alpha_VSA3 Alpha_VSA4 
                           Beta_VSA0 Beta_VSA1 Beta_VSA2 Beta_VSA3 Beta_VSA4
    .
        
            - MaxAlpha, MaxBeta: maximum alpha and beta values in the molecule
            - MinAlpha, MinBeta: minimum alpha and beta values in the molecule
            - TotalAlpha, TotalBeta: total sum of all alpha and beta values in a molecule
            - AverageAlpha, AverageBeta: total sum of all alpha and beta values in a molecule, averaged by nr of identified functional groups
            - AlphaMolWtNormalised, BetaMolWtNormalised: total sum of all alpha and beta values in a molecule normalised by MolWt
            - AlphaVSANormalised, BetaVSANormalised: total sum of all alpha and beta values in a molecule normalised by van der Waals surface area
            - AlphaLogPNormalised, BetaLogPNormalised: total sum of all alpha and beta values in a molecule normalised by LogP
            - Alpha_VSA0, Alpha_VSA1, Alpha_VSA2, Alpha_VSA3, Alpha_VSA4: fraction of surface area with a value of alpha between x and y
            - Beta_VSA0, Beta_VSA1, Beta_VSA2, Beta_VSA3, Beta_VSA4: fraction of surface area with a value of beta between x and y
        '''
    
        ep_values = namedtuple('ep_values', 'MaxAlpha MaxBeta MinAlpha MinBeta TotalAlpha TotalBeta AverageAlpha AverageBeta AlphaMolWtNormalised BetaMolWtNormalised AlphaVSANormalised BetaVSANormalised AlphaLogPNormalised BetaLogPNormalised Alpha_VSA0 Alpha_VSA1 Alpha_VSA2 Alpha_VSA3 Alpha_VSA4 Beta_VSA0 Beta_VSA1 Beta_VSA2 Beta_VSA3 Beta_VSA4')
    
        mol = self._uncharge_molecule(mol)

        self._find_matches(mol)

        self._adjust_matches()

        ep_values.MaxAlpha, ep_values.MaxBeta, ep_values.MinAlpha, ep_values.MinBeta = self._max_min_alpha_beta()
        
        ep_values.TotalAlpha, ep_values.TotalBeta, ep_values.AverageAlpha, ep_values.AverageBeta = self._average_total_alpha_beta() 
        
        ep_values.AlphaMolWtNormalised, ep_values.BetaMolWtNormalised, ep_values.AlphaVSANormalised, ep_values.BetaVSANormalised, ep_values.AlphaLogPNormalised, ep_values.BetaLogPNormalised = self._normalised_alpha_beta(mol, ep_values.TotalAlpha, ep_values.TotalBeta)

        ep_values.Alpha_VSA0, ep_values.Alpha_VSA1, ep_values.Alpha_VSA2, ep_values.Alpha_VSA3, ep_values.Alpha_VSA4, ep_values.Beta_VSA0, ep_values.Beta_VSA1, ep_values.Beta_VSA2, ep_values.Beta_VSA3, ep_values.Beta_VSA4 = self._VSA_alpha_beta(mol) 

        return ep_values
        
        
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        