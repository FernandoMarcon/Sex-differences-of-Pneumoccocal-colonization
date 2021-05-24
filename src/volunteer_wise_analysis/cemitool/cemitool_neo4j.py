import pandas as pd
from py2neo import Graph

graph = Graph("bolt://localhost:7687", user="neo4j", password="cemitool")

basedir = 'intermediate/volunteer_wise_analysis/cemitool/carriage_sex/Tables/'
# df = pd.read_csv('/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization/intermediate/volunteer_wise_analysis/cemitool/carriage_sex/Tables/module.tsv', delimiter='\t')
interactions = pd.read_csv(basedir+'interactions.tsv', delimiter='\t')
interactions.head()

query = '''
        MERGE (g1:Gene {name: $gene1})
        MERGE (g2:Gene {name: $gene2})
        MERGE (m:Module {name: $mod})
        CREATE (g1) - [:MEMBER_OF] -> (m)
        CREATE (g2) - [:MEMBER_OF] -> (m)
        CREATE (g1) - [:PPI] -> (g2)
        '''
tx = graph.begin()
for index, row in interactions.iterrows():
    tx.run(query, parameters = {'gene1':row['Gene1'],'gene2':row['Gene2'], 'mod':row['Module']})
tx.commit()


nes = pd.read_csv(basedir+'enrichment_nes.tsv',delimiter='\t')
nes.head()
nes_long = pd.melt(nes, id_vars=['pathway'], var_name='group', value_name='nes')
nes_long = nes_long[nes_long['pathway'] != 'Not.Correlated']
nes_long.head()
query = '''
        MERGE (g:Module {name: $mod})
        MERGE (h:Group {name: $group})
        CREATE (h) -[:NES {value: $nes}]-> (g)
        '''
tx = graph.begin()
for index, row in nes_long.iterrows():
    tx.run(query, parameters = {'mod':row['pathway'],'group':row['group'], 'nes':row['nes']})
tx.commit()


ora = pd.read_csv(basedir+'ora.tsv', delimiter='\t')
ora = ora[ora['p.adjust'] < 0.01]
ora.Module.unique()
ora.head()

tx = graph.begin()
for index, row in ora.iterrows():
    tx.run('''
            MERGE (m:Module {name:$mod})
            MERGE (p:Pathway {name:$path})
            CREATE (p)-[:ENRICHED]->(m)
            ''', parameters = {'mod':row['Module'], 'path':row['ID']})
tx.commit()
# a
