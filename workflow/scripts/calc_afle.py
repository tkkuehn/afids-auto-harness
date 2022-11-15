import csv
import json

import numpy as np
from numpy.linalg import norm


FCSV_FIELDS = (
    "id",
    "x",
    "y",
    "z",
    "ow",
    "ox",
    "oy",
    "oz",
    "vis",
    "sel",
    "lock",
    "label",
    "desc",
    "associatedNodeID",
)


def read_fcsv(fcsv_path):
    with open(fcsv_path, "r", encoding="utf-8") as fcsv_file:
        return {
            row["desc"].lower(): np.array(
                [float(row["x"]), float(row["y"]), float(row["z"])]
            )
            for row in csv.DictReader(fcsv_file, fieldnames=FCSV_FIELDS)
            if row["desc"] not in ["desc", None]
        }


def calc_afles(ground_truth_path, model_path):
    ground_truth_afids = read_fcsv(ground_truth_path)
    model_afids = read_fcsv(model_path)

    return {
        key: norm(ground_truth_afids[key] - model_afids[key])
        for key in ground_truth_afids.keys()
    }


def main():
    afles = calc_afles(
        snakemake.input["ground_truth_fcsv"], snakemake.input["model_fcsv"]
    )
    afles["mean"] = np.mean(np.array(list(afles.values())))
    with open(snakemake.output["afle"], "w", encoding="utf-8") as afle_file:
        json.dump(afles, afle_file)


if __name__ == "__main__":
    main()
