import pandas as pd
from py2neo import Graph

df = pd.read_csv('/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization/intermediate/volunteer_wise_analysis/cemitool/carriage_sex/Tables/module.tsv', delimiter='\t')
df.head()

graph = Graph("bolt://localhost:7687", user="neo4j", password="cemitool")

query = '''
        MERGE (e:Gene {name: $gene})
        MERGE (p:Module {name: $mod})
        CREATE (e) - [:MEMBER_OF] -> (p)
        '''
tx = graph.begin()
for index, row in df.iterrows():
    tx.run(query, parameters = {'gene':row['genes'], 'mod':row['modules']})
tx.commit()
