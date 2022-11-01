import math
from os import chmod, PathLike
from pathlib import Path
from shutil import copy
import stat
from typing import List
import random

random.seed(snakemake.params["seed"])

train_indices = random.sample(
    range(len(snakemake.input["input_images"])),
    math.floor(
        len(snakemake.input["input_images"]) * snakemake.params["train_proportion"]
    ),
)
test_indices = list(
    set(range(len(snakemake.input["input_images"]))) - set(train_indices)
)

train_images = [snakemake.input["input_images"][idx] for idx in train_indices]
train_fcsvs = [snakemake.input["ground_truth"][idx] for idx in train_indices]
test_images = [snakemake.input["input_images"][idx] for idx in test_indices]
test_fcsvs = [snakemake.input["ground_truth"][idx] for idx in test_indices]

train_out = Path(snakemake.output["training_data"])
train_out.mkdir(parents=True)
test_out = Path(snakemake.output["testing_data"])
test_out.mkdir(parents=True)


def copy_subset(files: List[str], src_root: PathLike, out_root: PathLike):
    for file_ in files:
        src_path = Path(file_)
        dest_path = out_root / src_path.relative_to(src_root)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        copy(src_path, dest_path, follow_symlinks=True)
        chmod(
            dest_path,
            stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH,
        )


for fileset, out_dir in zip(
    [train_images, train_fcsvs, test_images, test_fcsvs],
    [train_out, train_out, test_out, test_out],
):
    copy_subset(fileset, snakemake.config["bids_dir"], out_dir)
