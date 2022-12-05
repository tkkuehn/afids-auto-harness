import math
from os import chmod, PathLike
from pathlib import Path
from shutil import copy
import stat
from typing import List

import numpy as np
from sklearn.model_selection import train_test_split

train_images, test_images, train_groundtruth, test_groundtruth = train_test_split(
    np.array(
        snakemake.input["hcp_images"]
        + snakemake.input["oasis_images"]
        + snakemake.input["snsx_images"]
    ),
    np.array(
        snakemake.input["hcp_groundtruth"]
        + snakemake.input["oasis_groundtruth"]
        + snakemake.input["snsx_groundtruth"]
    ),
    random_state=1,
    stratify=np.array(
        ["hcp" for _ in snakemake.input["hcp_images"]]
        + ["oasis" for _ in snakemake.input["oasis_images"]]
        + ["snsx" for _ in snakemake.input["snsx_images"]]
    ),
)

train_out = Path("output/train")
train_out.mkdir(parents=True, exist_ok=True)
test_out = Path("output/test")
test_out.mkdir(parents=True, exist_ok=True)


def copy_subset(files: List[str], out_root: PathLike):
    for file_ in files:
        if snakemake.config["hcp_dir"] in file_:
            root_dir = snakemake.config["hcp_dir"]
            ds_name = "hcp"
        elif snakemake.config["snsx_dir"] in file_:
            root_dir = snakemake.config["snsx_dir"]
            ds_name = "snsx"
        elif snakemake.config["oasis_dir"] in file_:
            root_dir = snakemake.config["oasis_dir"]
            ds_name = "oasis"
        else:
            raise ValueError(f"Unknown root for {fileset}")

        src_path = Path(file_)
        dest_path = out_root / ds_name / src_path.relative_to(root_dir)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        copy(src_path, dest_path, follow_symlinks=True)
        chmod(
            dest_path,
            stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH,
        )


for fileset, out_dir in zip(
    [train_images, train_groundtruth, test_images, test_groundtruth],
    [train_out, train_out, test_out, test_out],
):
    copy_subset(fileset, out_dir)
