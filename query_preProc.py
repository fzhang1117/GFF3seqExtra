## query_preProc.py ##
## pre processing the query file ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2022-08-18 ##

import sys
import pandas as pd

fl_query = sys.argv[1]
query = pd.read_table(fl_query, sep = "\t", header = None)
clusterList = pd.unique(query[1])
clusternum = len(clusterList)

for cluster in clusterList:
    geneList = query[query[1] == cluster][0].tolist()
    geneNum = len(geneList)
    fl_output = './tmp/' + cluster + '_' + str(geneNum) + '_genes.txt'
    with open(fl_output, 'w') as fh_output:
        fh_output.write('\n'.join(geneList))
