#! usr/bin/env python

# a collection of pre_processsing(), pred_power_calculator() and evalLogicGate*() functions for calculation of results of N-gene logicGates.

# Define global parameter
specificity_cutoff = 0.9 # 0.8 0.9 should use 0.9 to match average safety score of clinical targets!!!

def pre_processing(fnIn1, fnIn2, fnGene):
    # Read data
    print('  Reading in candidate gene list ...')
    gene_candi = pd.read_csv(fnGene, sep='\t', header=None).iloc[:, 0].tolist()
    print('  Reading in TME raw expression matrix ...')
    data1 = pd.read_csv(fnIn1, sep='\t', skiprows=[0, 1, 2], index_col=0, header=None)
    gene_list1 = data1.index.tolist()
    print('  Reading in normal raw expression matrix ...')
    data2 = pd.read_csv(fnIn2, sep='\t', skiprows=[0, 1], index_col=0, header=None)
    gene_list2 = data2.index.tolist()
    set2 = set(gene_list1)
    set3 = set(gene_list2)
    gene_candi = [gene for gene in gene_candi if gene in set2 and gene in set3]
    data1 = data1.loc[gene_candi]
    data2 = data2.loc[gene_candi]
    del gene_list1, gene_list2, set2, set3  # Destroy unnecessary variables
    data = pd.concat([data1, data2], axis=1)
    del data1, data2  # Destroy unnecessary variables
    # Binarize
    print('  Binarizing the expression matrix ...')
    data = (data > 0)
    data_new = ~data
    data_new.index = [f"{idx}_not" for idx in data.index]
    data = pd.concat([data, data_new])
    del data_new  # Destroy unnecessary variable
    # Transpose
    data = data.T
    # Read labels
    labels1 = pd.read_csv(fnIn1, sep='\t', nrows=3, header=None, index_col=0)
    labels2 = pd.read_csv(fnIn2, sep='\t', nrows=2, header=None, index_col=0)
    labels2.iloc[[0, 1]] = labels2.iloc[[1, 0]].values
    new_row = pd.DataFrame([["Normal"] * labels2.shape[1]], index=["NewIndex"])
    new_row.columns = range(1, labels2.shape[1] + 1)
    labels2 = pd.concat([new_row, labels2])
    labels1.reset_index(drop=True, inplace=True)
    labels2.reset_index(drop=True, inplace=True)
    labels = pd.concat([labels1, labels2], axis=1)
    del labels1, labels2, new_row  # Destroy unnecessary variables
    # Set labels
    cell_label = list(labels.iloc[2])
    cell_label = [False if ('Malignant' not in c) and ('Tumor' not in c) else True for c in cell_label]
    patient_label = list(labels.iloc[1])
    cancerType_label = list(labels.iloc[0])
    del labels  # Destroy unnecessary variable
    # Add labels to data
    data['cell_type'] = cell_label
    data['patient'] = patient_label
    data['cancer_type'] = cancerType_label
    print('Total cell number %d, tumor cell number %d' % (len(cell_label), sum(cell_label)))
    del cell_label, patient_label, cancerType_label, gene_candi  # Destroy unnecessary variables
    return data

# def pre_processing(fnIn1, fnIn2, fnGene):
#     # Read data
#     print('  Reading in candidate gene list ...')
#     gene_candi = pd.read_csv(fnGene, sep='\t', header=None)
#     gene_candi = gene_candi.iloc[:, 0].tolist()
#     print('  Reading in TME raw expression matrix ...')
#     data1 = pd.read_csv(fnIn1, sep='\t', skiprows=[0, 1, 2], index_col=0, header=None)
#     gene_list1 = data1.index.tolist()
#     print('  Reading in normal raw expression matrix ...')
#     data2 = pd.read_csv(fnIn2, sep='\t', skiprows=[0, 1], index_col=0, header=None)
#     gene_list2 = data2.index.tolist()
#     set2 = set(gene_list1)
#     set3 = set(gene_list2)
#     gene_candi = [gene for gene in gene_candi if gene in set2 and gene in set3]
#
#     data1 = data1.loc[gene_candi]
#     data2 = data2.loc[gene_candi]
#     data = pd.concat([data1, data2], axis=1)
#
#     # Binarize
#     print('  Binarizing the expression matrix ...')
#     data = (data>0)
#     data_new = ~data
#     data_new.index = [f"{idx}_not" for idx in data.index]
#     data = pd.concat([data, data_new])
#     # Transpose
#     data = data.T
#     # Read labels
#     labels1 = pd.read_csv(fnIn1, sep='\t', nrows=3, header=None, index_col=0)
#     labels2 = pd.read_csv(fnIn2, sep='\t', nrows=2, header=None, index_col=0)
#     labels2.iloc[[0, 1]] = labels2.iloc[[1, 0]].values
#     new_row = pd.DataFrame([["Normal"] * labels2.shape[1]], index=["NewIndex"])
#     new_row.columns = range(1, labels2.shape[1] + 1)
#     labels2 = pd.concat([new_row, labels2])
#     labels1.reset_index(drop=True, inplace=True)
#     labels2.reset_index(drop=True, inplace=True)
#     labels = pd.concat([labels1, labels2], axis=1)
#     # Set labels
#     cell_label = list(labels.iloc[2])
#     cell_label = [False if ('Malignant' not in c) and ('Tumor' not in c) else True for c in cell_label]
#     patient_label = list(labels.iloc[1])
#     cancerType_label = list(labels.iloc[0])
#     # Add labels to data
#     data['cell_type'] = cell_label
#     data['patient'] = patient_label
#     data['cancer_type'] = cancerType_label
#     print('Total cell number %d, tumor cell number %d'%(len(cell_label),sum(cell_label)))
#     return data

def pred_power_calculator(selected_genes):
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power

############ 1-gene logic combinations
def evalLogicGate_zero(individual):
    selected_genes = df.iloc[:,individual[0]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

############ 2-gene logic combinations
def evalLogicGateA(individual):
    selected_genes = df.iloc[:,individual[0]] & df.iloc[:,individual[1]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO(individual):
    selected_genes = df.iloc[:,individual[0]] | df.iloc[:,individual[1]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

############ 3-gene logic combinations
def evalLogicGateAA(individual):
    selected_genes = (df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO(individual):
    selected_genes = (df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO(individual):
    selected_genes = df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

############ 4-gene logic combinations
def evalLogicGateAA_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_A(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAA_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_O(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateA_A_O(individual):
    selected_genes = (df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & (df.iloc[:,individual[2]] | df.iloc[:,individual[3]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateA_O_A(individual):
    selected_genes = (df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | (df.iloc[:,individual[2]] & df.iloc[:,individual[3]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO_A_O(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & (df.iloc[:,individual[2]] | df.iloc[:,individual[3]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO_O_A(individual):
    selected_genes = df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | (df.iloc[:,individual[2]] & df.iloc[:,individual[3]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,






############ 5-gene logic combinations
def evalLogicGateAA_A_A(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_A_A(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_A_A(individual):
    selected_genes = (((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_A_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAA_O_A(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_O_A(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_O_A(individual):
    selected_genes = (((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_O_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateA_A_O_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & (df.iloc[:,individual[2]] | df.iloc[:,individual[3]])) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateA_O_A_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | (df.iloc[:,individual[2]] & df.iloc[:,individual[3]])) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO_A_O_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & (df.iloc[:,individual[2]] | df.iloc[:,individual[3]])) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO_O_A_A(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | (df.iloc[:,individual[2]] & df.iloc[:,individual[3]])) & df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,



def evalLogicGateAA_A_O(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_A_O(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_A_O(individual):
    selected_genes = (((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_A_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) & df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAA_O_O(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_O_O(individual):
    selected_genes = (((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_O_O(individual):
    selected_genes = (((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_O_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) | df.iloc[:,individual[3]]) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateA_A_O_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & (df.iloc[:,individual[2]] | df.iloc[:,individual[3]])) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateA_O_A_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | (df.iloc[:,individual[2]] & df.iloc[:,individual[3]])) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO_A_O_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & (df.iloc[:,individual[2]] | df.iloc[:,individual[3]])) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateO_O_A_O(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | (df.iloc[:,individual[2]] & df.iloc[:,individual[3]])) | df.iloc[:,individual[4]]
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,



def evalLogicGateAA_a_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_a_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_a_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_a_A(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,



def evalLogicGateAA_o_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_o_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_o_A(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_o_A(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] & df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,



def evalLogicGateAA_a_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_a_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_a_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_a_O(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) & (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,



def evalLogicGateAA_o_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateAO_o_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] & df.iloc[:,individual[1]]) | df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOA_o_O(individual):
    selected_genes = ((df.iloc[:,individual[0]] | df.iloc[:,individual[1]]) & df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,

def evalLogicGateOO_o_O(individual):
    selected_genes = (df.iloc[:,individual[0]] | df.iloc[:,individual[1]] | df.iloc[:,individual[2]]) | (df.iloc[:,individual[3]] | df.iloc[:,individual[4]])
    not_selected_genes = (selected_genes == False)
    tn = np.sum(not_true_labels & not_selected_genes)
    fp = np.sum(not_true_labels & selected_genes)
    fn = np.sum(true_labels & not_selected_genes)
    tp = np.sum(true_labels & selected_genes)
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0
    pred_power = sensitivity * (specificity >= specificity_cutoff)
    return pred_power,
