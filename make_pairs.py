import sys
import json
import pandas as pd
from pathlib import Path

def main(tree_path):
    path = Path(tree_path)

    bin_tree = json.load(open(path))

    # Convert json tree into more friendly format
    leaves = {} # gvid : accession
    for o in bin_tree["objects"]:
        leaves[o["_gvid"]] = o["name"]

    # [(parent, child), ...]
    pairs = []
    for edge in bin_tree["edges"]:
        pairs.append({
            "parent": leaves[edge["tail"]],
            "child": leaves[edge["head"]]
        })    

    df_pairs = pd.DataFrame.from_dict(pairs)
    df_pairs.to_csv(path.with_name("pairs.csv"), index=False)

if __name__ == "__main__":

    ## Convert tree.dot to tree.json
    # dot -Txdot_json -o tree.json tree.dot

    # sys.argv[1] = /path/to/tree.json
    main(sys.argv[1])