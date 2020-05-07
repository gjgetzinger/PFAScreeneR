import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
import molvs
import sys
import os
import subprocess
import tempfile
import fnmatch
import xlrd
import sqlite3
from molvs.tautomer import TautomerTransform
from rdkit import RDLogger
from progress.bar import Bar
from time import sleep
from time import process_time


TAUTOMER_TRANSFORMS = (
    TautomerTransform('1,3 (thio)keto/enol f', '[CX4!H0]-[C]=[O,S,Se,Te;X1]'),
    TautomerTransform('1,3 (thio)keto/enol r', '[O,S,Se,Te;X2!H0]-[C]=[C]'),
    TautomerTransform('1,5 (thio)keto/enol f', '[CX4,NX3;!H0]-[C]=[C][CH0]=[O,S,Se,Te;X1]'),
    TautomerTransform('1,5 (thio)keto/enol r', '[O,S,Se,Te;X2!H0]-[CH0]=[C]-[C]=[C,N]'),
    TautomerTransform('aliphatic imine f', '[CX4!H0]-[C]=[NX2]'),
    TautomerTransform('aliphatic imine r', '[NX3!H0]-[C]=[CX3]'),
    TautomerTransform('special imine f', '[N!H0]-[C]=[CX3R0]'),
    TautomerTransform('special imine r', '[CX4!H0]-[c]=[n]'),
    TautomerTransform('1,3 aromatic heteroatom H shift f', '[#7!H0]-[#6R1]=[O,#7X2]'),
    TautomerTransform('1,3 aromatic heteroatom H shift r', '[O,#7;!H0]-[#6R1]=[#7X2]'),
    TautomerTransform('1,3 heteroatom H shift', '[#7,S,O,Se,Te;!H0]-[#7X2,#6,#15]=[#7,#16,#8,Se,Te]'),
    TautomerTransform('1,5 aromatic heteroatom H shift', '[#7,#16,#8;!H0]-[#6,#7]=[#6]-[#6,#7]=[#7,#16,#8;H0]'),
    TautomerTransform('1,5 aromatic heteroatom H shift f', '[#7,#16,#8,Se,Te;!H0]-[#6,nX2]=[#6,nX2]-[#6,#7X2]=[#7X2,S,O,Se,Te]'),
    TautomerTransform('1,5 aromatic heteroatom H shift r', '[#7,S,O,Se,Te;!H0]-[#6,#7X2]=[#6,nX2]-[#6,nX2]=[#7,#16,#8,Se,Te]'),
    TautomerTransform('1,7 aromatic heteroatom H shift f', '[#7,#8,#16,Se,Te;!H0]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6]-[#6,#7X2]=[#7X2,S,O,Se,Te,CX3]'),
    TautomerTransform('1,7 aromatic heteroatom H shift r', '[#7,S,O,Se,Te,CX4;!H0]-[#6,#7X2]=[#6]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[NX2,S,O,Se,Te]'),
    TautomerTransform('1,9 aromatic heteroatom H shift f', '[#7,O;!H0]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#7,O]'),
    TautomerTransform('1,11 aromatic heteroatom H shift f', '[#7,O;!H0]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#7X2,O]'),
    TautomerTransform('furanone f', '[O,S,N;!H0]-[#6r5]=[#6X3r5;$([#6]([#6r5])=[#6r5])]'),
    TautomerTransform('furanone r', '[#6r5!H0;$([#6]([#6r5])[#6r5])]-[#6r5]=[O,S,N]'),
    TautomerTransform('keten/ynol f', '[C!H0]=[C]=[O,S,Se,Te;X1]', bonds='#-'),
    TautomerTransform('keten/ynol r', '[O,S,Se,Te;!H0X2]-[C]#[C]', bonds='=='),
    TautomerTransform('ionic nitro/aci-nitro f', '[C!H0]-[N+;$([N][O-])]=[O]'),
    TautomerTransform('ionic nitro/aci-nitro r', '[O!H0]-[N+;$([N][O-])]=[C]'),
    TautomerTransform('oxim/nitroso f', '[O!H0]-[N]=[C]'),
    TautomerTransform('oxim/nitroso r', '[C!H0]-[N]=[O]'),
    TautomerTransform('oxim/nitroso via phenol f', '[O!H0]-[N]=[C]-[C]=[C]-[C]=[OH0]'),
    TautomerTransform('oxim/nitroso via phenol r', '[O!H0]-[c]=[c]-[c]=[c]-[N]=[OH0]'),
    TautomerTransform('cyano/iso-cyanic acid f', '[O!H0]-[C]#[N]', bonds='=='),
    TautomerTransform('cyano/iso-cyanic acid r', '[N!H0]=[C]=[O]', bonds='#-'),
    # TautomerTransform('formamidinesulfinic acid f', '[O,N;!H0]-[C]=[S,Se,Te]=[O]', bonds='=--'),  # TODO: WAT!?
    # TautomerTransform('formamidinesulfinic acid r', '[O!H0]-[S,Se,Te]-[C]=[O,N]', bonds='=--'),
    TautomerTransform('isocyanide f', '[C-0!H0]#[N+0]', bonds='#', charges='-+'),
    TautomerTransform('isocyanide r', '[N+!H0]#[C-]', bonds='#', charges='-+'),
    TautomerTransform('phosphinic acid f', '[OH]-[PD3X3H0]', bonds='='),
    TautomerTransform('phosphinic acid r', '[PD3X3H1]=[O]', bonds='-'),
)

s = molvs.Standardizer(tautomer_transforms=TAUTOMER_TRANSFORMS, prefer_organic=True)

HYD_RXNS = dict({
'hyd_Xali_SN_A':Chem.AllChem.ReactionFromSmarts('[$([CX4]);$(C[Cl,Br,I]);!$(C([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]);!$(C([F,Cl,Br,I])[F,Cl,Br,I]):1][Cl,Br,I:2]>>[*:1][OH1].[*:2]'),
'hyd_Xali_SN_B':Chem.AllChem.ReactionFromSmarts('[$([CX4]);$(C[Cl,Br,I]);$(CC[Cl,Br,I]);!$(C([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]);!$(C([F,Cl,Br,I])[F,Cl,Br,I]):1][Cl,Br,I:2]>>[*:1][OH1].[*:2]'),
'hyd_Xali_SN_C':Chem.AllChem.ReactionFromSmarts('[Cl,Br,I:5]-[CX4:1]([F,Cl,Br,I:3])([*:4])([*:2])>>[*:5].[OH1][A:1]([*:3])([*:4])([*:2])'),
'hyd_Xali_E':Chem.AllChem.ReactionFromSmarts('[F,Cl,Br,I:3]-[CX4:1]-[C;X4;H1,H2:2]-[#1:4]>>[C:1]=[C:2].[*:4][*:3]'),
'hyd_epox':Chem.AllChem.ReactionFromSmarts('([OX2r3]1[#6r3:1][#6r3:2]1)>>([OX2r3][*:1][*:2][OH1])'),
'hyd_OPE_1':Chem.AllChem.ReactionFromSmarts('[PX4D4:1](=[!#6,O,S:2])([N,O,S;H0:3])([N,O,S;H0:4])[N,O,S;H0:5]>>[*:1]([OH1])(=[*:2])([*:3])[*:4].[H][*:5]'),
'hyd_ester':Chem.AllChem.ReactionFromSmarts('[CX3;$([#6]),$([H1]):1](=[OX1:4])[OX2:2][#6;!$(C=[O,N,S]):3]>>[*:1](=[*:4])[OH1].[*H1:2][*:3]'),
'hyd_lactone':Chem.AllChem.ReactionFromSmarts('([C][CX3R:1](=[OX1])[#8X2:2][C;!$(C=[O,N,S])])>>([*:1][OH1].[H][*:2])'),
'hyd_carbonate':Chem.AllChem.ReactionFromSmarts('[#6;!$(C=[O,N,S])][#8X2:1][#6X3:2](=[OX1])[#8X2:3][#6;!$(C=[O,N,S])]>>[*H1:1].[*:2].[H][*:3]'),
'hyd_anhydride':Chem.AllChem.ReactionFromSmarts('[O:1]([C](=O)[C])([C:2](=O)[C])>>[OH1:1]([C](=O)[C]).[OH1][C:2](=O)[C]'),
'hyd_cyanhydride':Chem.AllChem.ReactionFromSmarts('([$([CR]=O),$([CR]=[OX1]):1][#8X2R:2][$([CR]=O),$([CR]=[OX1]):3])>>([*:1][OH1].[OH1][*:3]).[*:2]'),
'hyd_amide':Chem.AllChem.ReactionFromSmarts('[CX3;$([R0][#6]),$([H1R0]):1](=[OX1:3])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]):2]>>[OH1][*:1](=[*:3]).[H][*:2]'),
'hyd_lactam':Chem.AllChem.ReactionFromSmarts('([CR:4][CX3R:1](=[OX1])[#7X3;$([H1][C;!$(C=[O,N,S])]),$([H0]([C;!$(C=[O,N,S])])[C;!$(C=[O,N,S])]):2][CR:3])>>([*:4][*:1](=O)[OH1].[*:2][*:3])'),
'hyd_carbamate':Chem.AllChem.ReactionFromSmarts('[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4]):1][C:2](=O)[O:3][$([#6]);!$(C=[O,S,N])]>>[H][*:1].[*:2].[H][*:3]'),
'hyd_thiocarbamate':Chem.AllChem.ReactionFromSmarts('[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4]):1][C:2](=O)[SX2:3][$([#6]);!$(C=[O,S,N])]>>[H][*:1].[*:2].[H][*:3]'),
'hyd_urea':Chem.AllChem.ReactionFromSmarts('[#7X3;!$([#7][!#6]):1][CX3:2](=[OX1])[#7X3;!$([#7][!#6]):3]>>[H][*:1].[*:2].[H][*:3]'),
'hyd_sulfonylurea':Chem.AllChem.ReactionFromSmarts('[SX4;$([H1]),$([H0][#6])](=[!#6])(=[!#6])[#7X3;!$([#7][!#6]):1][#6X3:2](=[OX1])[#7X3;!$([#7][!#6]):3]>>[H][*:1].[*:2].[H][*:3]'),
'hyd_nitrile':Chem.AllChem.ReactionFromSmarts('[NX1:1]#[CX2:2]>>[H][*:1]([H])[*:2](=O)'),
'hyd_NS':Chem.AllChem.ReactionFromSmarts('[N:1][SX2:2][C,N,O:3]>>[H][*:1].[OH1][*:2][*:3]'),
'hyd_imdide':Chem.AllChem.ReactionFromSmarts('([#6X3;$([H0][#6]),$([H1]):1](=[OX1])[#7X3:2]([H,C,N])[#6X3;$([H0][#6]),$([H1])](=[OX1]))>>([*:1][OH1].[H][*:2])'),
'hyd_acidhalide':Chem.AllChem.ReactionFromSmarts('[$([S,C]=[O,S]):1][F,Br,Cl,I:2]>>[*:1][OH1].[*:2]'),
'hyd_gemdiol':Chem.AllChem.ReactionFromSmarts('[OX2H1:1][CX4;!$(C([OX2H])([OX2H])[O,S,#7,#15]):2][OX2H1:3]>>[*X1H0:1]=[*:2].[*:3]')
})

def hydrolyze(mol = None,nsteps = 1,hyd_rxns = None):
  """
  Hydrolyze a molecule. 
  
  Subject an input molecule to hydrolyze. Multistep reactions (nsteps>1) 
  hydrolyze products until the specified number of reaction steps is achieved. 
    
  Parameters
  ----------
  mol: rdkit.Chem.rdchem.Mol
    An rdkit mol object. 
  nsteps: int
    The number of reaction steps to apply. 
  filename: str
    A filename where hydrolysis reactions are stored as SMARTS.
    
  Returns
  -------
  rdkit.Chem.rdchem.Mol 
    A list of rdkit mol objects.
  """
  
  if not isinstance(mol, Chem.rdchem.Mol):
    sys.exit('Input mol is not of type rdkit.Chem.rdchem.Mol')
  if mol is None:
    sys.exit('Missing input molecule.')
  if not isinstance(hyd_rxns, dict):
    sys.exit('hyd_rxns should be a dictionary')
  if not np.all([isinstance(hyd_rxns[key],rdkit.Chem.rdChemReactions.ChemicalReaction) for key in list(hyd_rxns.keys())]):
    sys.exit('Each element of hyd_rxns should be type rdkit.Chem.rdChemReaction.ChemicalReaction')
  n = 1
  products = {}
  while n <= nsteps:
    if n == 1:
      mols = [mol]
    else:
      if len(products)>0:
        mols = list(products.values())
      else:
        break 
    for m in mols:
      prods = {k:list(np.concatenate(v)) for k,v in {k:hyd_rxns[k].RunReactants((m,)) for k in hyd_rxns}.items() if len(v)>0}
      for rxn,molprod in prods.items():
        for mp in molprod:
          mp = s.super_parent(mp) 
          inchikey = Chem.MolToInchiKey(mp)
          mp.SetProp('Reaction', '\t'.join(map(str,[Chem.MolToInchiKey(m),inchikey,rxn, n])))
          if inchikey in [key for key in products]:
            products[inchikey] = MergeMolProps(products[inchikey],mp,True)
          else: 
            products[inchikey] = mp
    n+=1
  return(products)
  
def biotransform(sd_in=None, sd_out=None, nsteps = 1,jar_file=None):
  """
  Predict biotransformations.
  
  Use the biotransformer algorithm to predict environmental microbial 
    biotransformation products. 
  
  """
  sdfFile = open(sd_out, 'a')
  sdfFile_w = Chem.SDWriter(sdfFile)
  supp = Chem.SDMolSupplier(sd_in)
  n = 5000
  x = [range(len(supp))[i * n:(i + 1) * n] for i in range((len(range(len(supp))) + n - 1) // n )]
  for i in x:
    tmp_sdf_in = tempfile.NamedTemporaryFile("w+b",suffix=".sdf",delete=True).name
    tmp_sdf_out = tempfile.NamedTemporaryFile("w+b",suffix=".sdf",delete=True).name
    w = Chem.SDWriter(tmp_sdf_in)
    [w.write(supp[mol]) for mol in i]
    w.close()
    os.system('java -jar '+jar_file+' -b env -isdf '+tmp_sdf_in+' -k pred -osdf '+tmp_sdf_out+' -s '+str(nsteps))
    supp2 = Chem.SDMolSupplier(tmp_sdf_out)
    for m in supp2:
      if m is not None:
        try:
          sdfFile_w.write(m)
        except ValueError:
          continue
      else:
        continue
  print('Total metabolites stored:', sdfFile_w.NumMols())
  sdfFile_w.close()
  sdfFile.close()
    
def MergeMolProps(mol1,mol2,inchikeymatch = True):
  '''
  Merge properties of mols into one new mol.
  
  Example
  -------
  mol1 = Chem.MolFromSmiles('c1ccccc1')
  mol2 = Chem.MolFromSmiles('c1ccccc1')
  mol1.SetProp('name', 'benzene')
  mol2.SetProp('name', 'cyclohexatriene')
  mol1.SetProp('id', 'abcd')
  mol2.SetProp('ID', 'xyz')
  mol3 = MergeMolProps(mol1, mol2)
  mol3.GetPropsAsDict()
  '''
  if inchikeymatch and not Chem.MolToInchiKey(mol1)[:14] == Chem.MolToInchiKey(mol2)[:14]:
    sys.exit('InChIKeys do not match')
  else:
    if not mol1.GetPropsAsDict() == mol2.GetPropsAsDict():
      props = mergeDict(mol1.GetPropsAsDict(),mol2.GetPropsAsDict(),True)
      for prop in props:
        mol1.SetProp(prop, str(props[prop]))
  return mol1

def mergeDict(dict1, dict2,join=True):
  ''' Merge dictionaries and keep values of common keys in list'''
  dict3 = {**dict1, **dict2}
  v_types = {k:type(v) for k,v in dict3.items()}
  for key, value in dict3.items():
    if key in dict1 and key in dict2:
      if join:
        dict3[key] = ",".join(set([str(value) , str(dict1[key])]))
      else:
        dict3[key] = [str(value), str(dict1[key])]
  return dict3
  
def make_pfas_masslist(mollist_path, biotrans_jarfile, hyd_sdf, biotrans_in_sdf, biotrans_out_sdf, result_sdf, result_db):
  """
  Make a PFAS molecular databse
  
  Take input molecule lists, apply hydrolysis and biotransformation, 
  calculate properties, write results to sdf and sql database. 
  
  Parameters
  ----------
  mollist_path: path to a directory containing excel sheets exported from 
    EPA Chemicals Dasboard 
  biotrans_jarfile: Full path to a the biotransformer jar file 
  hyd_sdf: A file name for the sdf file that will store results of hydrolysis
    prediction. 
  biotrans_in_sdf: A filename for the input sdf file for biotransformation prediction. 
  biotrans_out_sdf: A filename for the outpur sdf file from biotransformation. 
  result_sdf: A file name for writing the resultant databsae to sdf file 
  result_db: A database filename for storing the results as a SQLite database. 
  """
  # define substructure patterns 
  carbon = Chem.MolFromSmarts('[#6]')
  CF2 = Chem.MolFromSmarts('[#6](F)F')
  fluorine = Chem.MolFromSmarts('[F]')
  
  # get a list of input files
  xls_files = {}
  for file in os.listdir(mollist_path):
    if fnmatch.fnmatch(file, '*.xls'):
      key = str.replace(file, '.xls', '')
      xls_files[key] = mollist_path+file 

  mols = []  
  for key in xls_files:
    # read the mol list 
    xls = xlrd.open_workbook(xls_files[key],encoding_override="gb2312")
    df = pd.read_excel(xls,na_values='-')
    # find the column with InChI string 
    inchi_col = df.columns[df.columns.to_series().str.contains('INCHI_STRING|InChI')]
    df = df.dropna(subset=inchi_col).reset_index()
    for inchi in df[inchi_col[0]]:
      mol = Chem.MolFromInchi(inchi, sanitize=False) # make mols from inchi 
      # split any non-covalent bonds 
      frags = Chem.GetMolFrags(mol, asMols = True, sanitizeFrags = False) 
      for f in frags:
        # get fragment, charge, isotope, stereochemistry and tautomer parent
        try:
          f_san = s.super_parent(f)
        except ValueError:
          continue 
        # keep carbon and CF2 containing molecules 
        if f_san.HasSubstructMatch(carbon) and f_san.HasSubstructMatch(CF2):
          f_san.SetProp('List', key) # set the input list 
          f_san.SetProp('InChIKey', Chem.MolToInchiKey(f_san))
          mols.append(f_san) # append to mols list 

  #make mols a dictionary keyed by inchikey 
  mols_dict = {}
  for mol in mols:
    inchikey = mol.GetProp('InChIKey')[:14]
    if inchikey in [key for key in mols_dict]:
      mols_dict[inchikey] = MergeMolProps(mols_dict[inchikey], mol)
    else:
      mols_dict[inchikey] = mol

# see which molecules have hits for hydrolyzable moieties 
  hyd_substruc = pd.DataFrame({}, columns = [key for key in HYD_RXNS])
  for inchikey in mols_dict:
    df = pd.DataFrame({key:mols_dict[inchikey].HasSubstructMatch(HYD_RXNS[key].GetReactantTemplate(0)) for key in HYD_RXNS}, index = [inchikey])
    hyd_substruc = hyd_substruc.append(df)

  hydrolyzable_inchikeys = hyd_substruc[hyd_substruc.any(axis = 'columns')].index.tolist()

  #predict hydrolysis products for each mol
  
  w = Chem.SDWriter(hyd_sdf)
  for inchikey in hydrolyzable_inchikeys:
    [w.write(v) for v in hydrolyze(mol=mols_dict[inchikey], nsteps=2, hyd_rxns=HYD_RXNS).values()]
  w.close()
 
  mols_rxn_hyd = pd.DataFrame({}, columns=['Precursor_InChIKey', 'Product_InChIKey', 'Reaction', 'RxnStep'])
  supp = Chem.SDMolSupplier(hyd_sdf)
  for mol in supp:
    if mol is None:
      continue 
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][F]')):
      continue 
    else:
      rxns = mol.GetProp('Reaction').split(',')
      rxn = [r.split('\t') for r in rxns]
      rxn_entry = pd.DataFrame({
        'Precursor_InChIKey':[r[0] for r in rxn],
        'Product_InChIKey':[r[1] for r in rxn], 
        'Reaction':[r[2] for r in rxn],
        'RxnStep':[r[3] for r in rxn]
      })
      mol.SetProp('List', 'HYD')
      inchikey = Chem.MolToInchiKey(mol)
      mol.SetProp('InChIKey', inchikey)
      mols_rxn_hyd = mols_rxn_hyd.append(rxn_entry)
      if inchikey[:14] in [key for key in mols_dict]:
        mols_dict[inchikey[:14]] = MergeMolProps(mols_dict[inchikey[:14]],mol)
      else:
        mols_dict[inchikey[:14]] = mol

  # predict biotrans products 
  w = Chem.SDWriter(biotrans_in_sdf)
  [w.write(mol) for mol in mols_dict.values()]
  w.close()

  biotransform(sd_in=biotrans_in_sdf, sd_out=biotrans_out_sdf, nsteps = 2,jar_file=biotrans_jarfile)

  # add biotransformation products to list
  mols_rxn_bio = []
  supp = Chem.SDMolSupplier(biotrans_out_sdf)
  for mol in supp:
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][F]')):
      continue 
    mol = MergeMolProps(s.super_parent(mol), mol,False)
    inchikey = Chem.MolToInchiKey(mol)
    inchikey_parent = mol.GetProp('Precursor_InChIKey')
    rxn = mol.GetProp('Reaction')
    for key in mol.GetPropsAsDict():
      if key not in ['Reaction','Precursor_InChIKey','InChIKey']:
        mol.ClearProp(key)
    mol.SetProp('InChIKey', inchikey)
    mol.SetProp('List', 'BIOTRANS')
    mol.SetProp('BioTrans_Reaction', rxn)
    mol.SetProp('Reaction', '\t'.join(map(str,[inchikey_parent,inchikey,rxn])))
    mols_rxn_bio.append(mol.GetPropsAsDict())
    if inchikey[:14] in [key for key in mols_dict]:
      mols_dict[inchikey[:14]] = MergeMolProps(mols_dict[inchikey[:14]],mol,True)
    else:
      mols_dict[inchikey[:14]] = mol

  mols_rxn_bio = pd.DataFrame(mols_rxn_bio)

  # clean up mol props 
  for inchikey in mols_dict:
    props = mols_dict[inchikey].GetPropsAsDict()
    for prop in props:
      mols_dict[inchikey].SetProp(prop, ",".join(set(str.split(props[prop],','))))

  # keep only one inchikey per entry 
  for inchikey in mols_dict:
    x = mols_dict[inchikey].GetProp('InChIKey')
    mols_dict[inchikey].SetProp('InChIKey_syn',x)
    mols_dict[inchikey].SetProp('InChIKey', x.split(',')[0])
  
  mol_props = []
  for inchikey in mols_dict:
    m = mols_dict[inchikey]
    mol_props.append({'InChIKey':inchikey,'MolForm':Chem.rdMolDescriptors.CalcMolFormula(m),'ExactMass':Chem.rdMolDescriptors.CalcExactMolWt(m),'FormalCharge':Chem.GetFormalCharge(m),'HBD':Chem.rdMolDescriptors.CalcNumHBD(m),'HBA':Chem.rdMolDescriptors.CalcNumHBA(m)})
  
  mol_props = pd.DataFrame(mol_props)

  # write molecules table 
  mols_tab = []
  list_tab = []
  for inchikey in mols_dict:
    mols_tab.append({'ID': inchikey, 
      'InChIKey':Chem.MolToInchiKey(mols_dict[inchikey]),
      'InChI':Chem.MolToInchi(mols_dict[inchikey]),
      'SMILES':Chem.MolToSmiles(mols_dict[inchikey]),
      'MolStructure':Chem.MolToMolBlock(mols_dict[inchikey])})
    lists = mols_dict[inchikey].GetProp('List').split(',')
    list_tab.append({
      'InChIKey': [Chem.MolToInchiKey(mols_dict[inchikey]) for l in lists], 
      'List': lists 
      })
  
  mols_tab_df = pd.DataFrame(mols_tab)
  list_tab_df = pd.DataFrame({}, columns = ['InChIKey','List'])
  for i in list_tab:
    list_tab_df = list_tab_df.append(pd.DataFrame(i))
  
  # write sdf 
  w = Chem.SDWriter(result_sdf)
  [w.write(mol) for mol in mols_dict.values()]
  w.close()

  # make a database 
  conn = sqlite3.connect(result_db)
  mols_tab_df.set_index(['ID','InChIKey']).to_sql("molecules", conn, if_exists='replace')
  mol_props.set_index('InChIKey').to_sql("mol_props", conn, if_exists='replace')
  list_tab_df.set_index('InChIKey').to_sql("list_membership", conn, if_exists='replace')
  mols_rxn_hyd.set_index(['Precursor_InChIKey','Product_InChIKey']).to_sql("hyd_rxns", conn, if_exists='replace')
  mols_rxn_bio.set_index(['Precursor_InChIKey','InChIKey']).to_sql("bio_rxns", conn, if_exists='replace')
  conn.close()
