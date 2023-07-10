# Import libraries
import pandas as pd
import numpy as np
import streamlit as st
import pickle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors

## Custom Function
# Calculate molecular descriptiors
def AromaticProportion(m):
    aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
    aa_count = []
    for i in aromatic atoms:
        if i == True:
            aa_count.append(1)
    AromaticAtom = sum(aa_count)
    HeavyAtom = Descriptors.HeavyAtomCount(m)
    AR = AromaticAtom/HeavyAtom
    return AR

def generate(smiles, verbose = False):
    
    moldata = []
    for elem in smiles:
        mol = chem.MolFromSmiles(elem)
        moldata.append(mol)
    
    baseData = np.arange(1,1)
    i = 0
    for mol in moldata:
        
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        desc_Aromatic Proportion = Aromatic Proportion(mol)
        
        row = np.array([desc_MolLogP,
                        desc_MolWt,
                        desc_NumRotatableBonds,
                        desc_AromaticProportion])
        
        if(i==0):
            baseData = row
        else:
            baseData=np.vstack([baseData, row])
        i = i+1
        
    columnNames = ['MolLogP','MolWt','NumRotatableBonds','AromaticProportion']
    descriptors = pd.DataFrame(data=baseData, columns = columnNames)
    
    return descriptors

# Page Title

image = Image.open("solubility-logo.jpg")

st.image(image, use_column_width = True)

st.write(""" 
# Molecular Solubility Prediction App

This app predicts the ** Solubility (LogS) ** values of molecules!

Data obtained from the John S. Delaney. [ESOL: Estimating Aqueous Solubility Directly from molecular Structure]

""")


# Input molecules (Side Panel)

st.sidebar.header('User Input features')

## Read Smiles Input

SMILES_input = "CCCCC\nCCC\nCN"
SMILES = st.sidebar.text_area("SMILES input", SMILES_input)
SMILES = "C\n" + SMILES # Adds C as a dummy, first item
SMILES = SMILES.split('\n')

st.header('Input SMILES')
SMILES[1:] # skips the dummy first item

## Calculate molecular descriptors
st.header('Computed molecular descriptors')
X = generate(SMILES)
X[1:] # skips the dummy first item

# Pre-built model

# Reads in saved Model
load_model = pickle.load(open('solubility_model.pkl', 'rb'))

# Apply model to make Predictions
prediction = load_model.predict(X)

st.header("Predicted LogS values")

prediction




