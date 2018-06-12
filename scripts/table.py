import json
from tabulate import tabulate
import pandas as pd

fmt="psql" #"pipe" #"psql"

with open('out.json') as f:
    d = json.load(f) # --> dict

    print("# Moves\n")
    for _d in d['moves']:
        df = pd.DataFrame(_d)
        a = tabulate(df,headers="keys", tablefmt=fmt)
        print(a, end="\n\n")

    for _d in d['energy'][0]['hamiltonian'][1]["nonbonded"]:
        df = pd.DataFrame(_d)
        a = tabulate(df,headers='keys', tablefmt=fmt)
        print(a, end="\n\n")


    print("# Analysis\n")
    for _d in d['analysis']:
        df = pd.DataFrame(_d)
        a = tabulate(df,headers="keys", tablefmt=fmt)
        print(a, end="\n\n")

    print("# Groups\n")
    _d = d['groups']
    df = pd.DataFrame(_d)
    a = tabulate(df,headers="keys", tablefmt=fmt)
    print(a, end="\n\n")
