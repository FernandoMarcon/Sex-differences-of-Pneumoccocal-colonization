import pandas as pd
from py2neo import Graph

df = pd.read_csv('/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization/intermediate/volunteer_wise_analysis/cemitool/carriage_sex/Tables/module.tsv', delimiter='\t')
df.head()

graph = Graph("bolt://localhost:7687", user="neo4j", password="cemitool")

tx = graph.begin()
for gene in df.genes.unique():
    tx.run('CREATE (n:Gene {name: $gene})', parameters = {'gene': gene})
for mod in df.modules.unique():
    tx.run('CREATE (n:Module {name: $module})', parameters = {'module': mod})
for index, row in df.iterrows():
    tx.run('''
       MATCH (a:Gene {name:$gene}), (b:Module {name:$module})
       MERGE (a)-[r:MEMBER_OF]->(b)
       ''', parameters = {'gene': row['genes'], 'module': row['modules']})
tx.commit()

query = '''
        MERGE (e:Gene), (p:Module)
        ON CREATE SET e.name = $gene, p.name = $mod
        ON MATCH SET e.name = $gene, p.name = $mod
        CREATE (e) - [:MEMBER_OF] -> (p)
        '''
tx = graph.begin()
for index, row in df.iterrows():
    tx.run(query, parameters = {'gene':row['genes'], 'mod':row['modules']})
tx.commit()

# for index, row in df.iterrows():
#     tx.run('''
#         CREATE (p1:Gene { name: $gene })
#         CREATE (p2:Module { name: $module })
#         CREATE (p1)-[:MEMBER_OF]->(p2)
#        ''', parameters = {'gene': row['genes'], 'module': row['modules']})
# tx.commit()
