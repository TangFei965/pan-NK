import scanpy as sc
import numpy as np
import pandas as pd
import os
from collections import Counter

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import permutation_test_score

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')

cd56 = sc.read_h5ad('cd56.h5ad')
cd56 = cd56[cd56.obs['meta_tissue_new']!='Other tissue']
t2t = {'Blood':1, 'Normal':0, 'Tumor':0}

choose56 = np.random.choice(cd56.obs['cellID'],2000,replace=False)

cd56_choose = cd56[cd56.obs['cellID'].isin(choose56)]
tissue56=[t2t[i] for i in cd56_choose.obs['meta_tissue_new'].tolist()]

clf = SVC(kernel="rbf", random_state=7)
cv = StratifiedKFold(2, shuffle=True, random_state=0)

score_rgs1, perm_scores_rgs1, pvalue_rgs1 = permutation_test_score(
clf, cd56_choose[:,['RGS1']].X.todense(), tissue56, scoring="accuracy", cv=cv,n_permutations=1000
)

score_cd69, perm_scores_cd69, pvalue_cd69 = permutation_test_score(
clf, cd56_choose[:,['CD69']].X.todense(), tissue56, scoring="accuracy", cv=cv,n_permutations=1000
)


fig, ax = plt.subplots()
ax.hist(perm_scores_rgs1, bins=20, density=True)
ax.axvline(score_rgs1, ls="--", color="r")
score_label = f"Score on original\ndata: {score_rgs1:.2f}\n(p-value: {pvalue_rgs1:.3f})"
ax.text(0.56, 92, score_label, fontsize=12)
ax.set_xlabel("Accuracy score (RGS1)")
ax.set_ylabel("Probability")
plt.grid(None)

plt.savefig('CD56_RGS1.pdf')
plt.show()


fig, ax = plt.subplots()
ax.hist(perm_scores_cd69, bins=20, density=True,color='orange')
ax.axvline(score_cd69, ls="--", color="r")
score_label = f"Score on original\ndata: {score_cd69:.2f}\n(p-value: {pvalue_cd69:.3f})"
ax.text(0.55, 82, score_label, fontsize=12)
ax.set_xlabel("Accuracy score (CD69)")
ax.set_ylabel("Probability")
plt.grid(None)

plt.savefig('CD56_CD69.pdf')
plt.show()