import math
from pathlib import Path
from shutil import copy
import random

random.seed(snakemake.params["seed"])

train_images = random.sample(snakemake.input["input_images"], math.floor(len(snakemake.input["input_images"]) * snakemake.params["train_proportion"]))
test_images = list(set(snakemake.input["input_images"]) - set(train_images))

train_out = Path(snakemake.output["training_data"])
train_out.mkdir(parents=True)
for image in train_images:
    src_path = Path(image)
    dest_path = train_out / src_path.relative_to(snakemake.config["bids_dir"])
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    copy(src_path, dest_path, follow_symlinks=True)

test_out = Path(snakemake.output["testing_data"])
test_out.mkdir(parents=True)
for image in test_images:
    src_path = Path(image)
    dest_path = test_out / src_path.relative_to(snakemake.config["bids_dir"])
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    copy(src_path, dest_path, follow_symlinks=True)
