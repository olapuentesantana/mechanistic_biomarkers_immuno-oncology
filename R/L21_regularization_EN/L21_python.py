
import pandas as pd
import numpy as np
import pprint
import os

from sklearn import linear_model
from sklearn.model_selection import train_test_split

# ***************
# Working directory
print(os.getcwd())
pprint = pprint.PrettyPrinter(indent=4)

# ***************
# declare variables
#l1_values = [.1, .5, .7, .9, .95, .99, 1]
l1_values = [1]
estimates = dict.fromkeys(range(0,100))
#pprint.pprint(estimates) 


# ***************
# data 
pathway_activity = pd.read_csv('./python_tmp/pathways_SKCM.csv', sep=',', header=0, index_col=0)
drug_response = pd.read_csv("./python_tmp/drugs_SKCM.csv",sep=',', header=0, index_col=0)

# X: explanatory data, Y: response data
X = pathway_activity
features = X.columns.values
Y = drug_response
drugs = Y.columns.values


# splitting data into training and test set
for i in range(0,1):
    estimates[i] = dict.fromkeys(l1_values)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.20, random_state=i)

    for j in l1_values:

        clf = linear_model.MultiTaskElasticNet(alpha=0,l1_ratio=j)
        clf.fit(X_train, Y_train)

        df = pd.DataFrame(clf.coef_.tolist(), columns = features, index = drugs)
        df.insert(loc = 0, column = 'intercept', value = clf.intercept_.tolist())
        #print.pprint(df)
        
        estimates[i][j] = df

pprint.pprint(estimates) 

