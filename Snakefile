import os
import pathlib
import pandas as pd


configfile: "files/config.yaml"


OUTDIR = config["output"]["dir"]
DATA_DIR = pathlib.Path(config["input"]["datadir"])

# assumed unique row identify linking to embryo name
samples = pd.read_csv(config["input"]["logfile"])
assert "file" in samples.columns
samples = samples.set_index(
    samples.apply(lambda x: x.file.replace("/", "_").split(".")[0], axis=1)
)
for emb in samples.file:
    emb_path = DATA_DIR.joinpath(emb)
    if not emb_path.exists():
        print(f"{emb_path} not found")

print("Samples:\n\t" + "\n\t".join(samples.index) + "\n")


rule all:
    input:
        os.path.join(OUTDIR, "final", "counts.csv"),


def get_embryo_param(wc, col):
    return samples.at[wc.embryo, col]


def get_image(wc):
    return DATA_DIR.joinpath(samples.loc[wc.embryo, "file"])


rule normalize_pmc_stains:
    input:
        image=lambda wc: get_image(wc),
    params:
        channel_name="pmc",
        channels=lambda wc: get_embryo_param(wc, "channel_order"),
        z_start=lambda wc: get_embryo_param(wc, "z_start"),
        z_end=lambda wc: get_embryo_param(wc, "z_end"),
    output:
        h5=temp(
            os.path.join(OUTDIR, "pmc_norm", "{embryo}.h5"),
        ),
    conda:
        "envs/hcr_quant.yaml"
    script:
        "scripts/normalize_pmc_stain.py"


rule predict_pmcs:
    input:
        image=os.path.join(OUTDIR, "pmc_norm", "{embryo}.h5"),
        model=config["ilastik"]["model"],
    params:
        ilastik_loc=config["ilastik"]["loc"],
    output:
        temp(os.path.join(OUTDIR, "pmc_probs", "{embryo}.h5")),
    log:
        os.path.join(OUTDIR, "logs", "prediction", "{embryo}.log"),
    shell:
        "({params.ilastik_loc} --headless "
        "--project={input.model} "
        "--output_format=hdf5 "
        "--output_filename_format={output} "
        "{input.image}) 2> {log}"


rule label_pmcs:
    input:
        stain=os.path.join(OUTDIR, "pmc_norm", "{embryo}.h5"),
        probs=os.path.join(OUTDIR, "pmc_probs", "{embryo}.h5"),
    output:
        labels=os.path.join(OUTDIR, "labels", "{embryo}_pmc_labels.h5"),
    log:
        log=os.path.join("logs", "labels", "{embryo}.log"),
    conda:
        "envs/segmentation.yaml"
    script:
        "scripts/label_pmcs.py"


rule quantify_expression:
    input:
        image=lambda wc: get_image(wc),
        labels=os.path.join(OUTDIR, "labels", "{embryo}_pmc_labels.h5"),
    params:
        gene_params=config["quant"]["genes"],
        channels=lambda wc: get_embryo_param(wc, "channel_order"),
        z_start=lambda wc: get_embryo_param(wc, "z_start"),
        z_end=lambda wc: get_embryo_param(wc, "z_end"),
        crop_image=True,
        quant_method="both",
    output:
        image=os.path.join(OUTDIR, "expression", "{embryo}.nc"),
        csv=os.path.join(OUTDIR, "counts", "{embryo}.csv"),
    log:
        "logs/quant/{embryo}.log",
    conda:
        "envs/hcr_quant.yaml"
    script:
        "scripts/count_spots.py"


rule combine_counts:
    input:
        expand(os.path.join(OUTDIR, "counts", "{embryo}.csv"), embryo=samples.index),
    output:
        os.path.join(OUTDIR, "final", "counts.csv"),
    script:
        "scripts/combine_counts.py"
