import json
import re

import pandas as pd


def process_afle_json(path):
    subject = re.search(r"sub-([a-zA-Z\d]+)_", path).group(1)
    with open(path, "r", encoding="utf-8") as json_file:
        afle_dict = json.load(json_file)
    return subject, afle_dict


def main():
    afle_paths = snakemake.input["afles"]
    afle_info = {subject: afle_dict for subject, afle_dict in [process_afle_json(path) for path in afle_paths]}
    columns = list(process_afle_json(afle_paths[0])[1].keys())
    pd.DataFrame(dict({column: [afle_dict[column] for afle_dict in afle_info.values()] for column in columns}, subjects=list(afle_info.keys()))).to_csv(snakemake.output["table"], sep="\t", index=False)


if __name__ == "__main__":
    main()
