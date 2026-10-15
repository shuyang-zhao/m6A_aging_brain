# Orthrus embeddings and RiboNN translation-efficiency predictions of the canonical mRNAs
# Needs a CUDA GPU. The stored products of the data release were made on one RTX 6000 Ada (torch 2.2.2+cu118) by an
# earlier version of this script that kept the model, torch and GPU names as separate fields of the npz
# reads the prepared transcript tables; writes transcripts_<species>.parquet (inputs), orthrus_<species>.npz,
# ribonn_input_<species>.txt and ribonn_<species>.tsv
#
#   embed_transcripts.py inputs  SPECIES TRANSCRIPTS EXONS SEQUENCES OUT   (CPU, writes the model input)
#   embed_transcripts.py orthrus SPECIES STAGED OUT                        (GPU)
#   embed_transcripts.py ribonn  SPECIES STAGED OUT RIBONN_DIR             (GPU; RiboNN checkout, weights in models/)
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODEL = "antichronology/orthrus-6-track"
MAX_TOKENS = 262144  # per batch
MAX_BATCH = 64
STOPS = {"TAA", "TAG", "TGA"}
STAGED = ["gene_id", "tx_id", "seq", "cds_start", "cds_end", "cds_end_nostop", "exon_lens"]
STAGED += ["start_ATG", "stop_ok", "frame_ok", "utr5_len", "utr3_len"]
ONE_HOT = np.full(256, 4)
ONE_HOT[[ord(base) for base in "ACGT"]] = range(4)


def stage(transcripts, exons, sequences):
    # join the transcript tables into the input of the GPU steps
    staged = transcripts.merge(sequences[["tx_id", "seq"]], on="tx_id", validate="1:1")
    assert (staged.seq.str.len() == staged.tx_len).all()
    lengths = exons.sort_values(["tx_id", "exon_no"]).groupby("tx_id")["len"]
    staged["exon_lens"] = staged.tx_id.map(lengths.agg(lambda sizes: ",".join(map(str, sizes))))
    staged["start_ATG"] = [seq[at : at + 3] == "ATG" for seq, at in zip(staged.seq, staged.cds_start, strict=True)]
    staged["stop_ok"] = [seq[end - 3 : end] in STOPS for seq, end in zip(staged.seq, staged.cds_end, strict=True)]
    staged["frame_ok"] = (staged.cds_end - staged.cds_start) % 3 == 0
    # the Orthrus CDS track excludes the stop codon
    staged["cds_end_nostop"] = staged.cds_end - 3 * staged.stop_ok
    return staged[STAGED].reset_index(drop=True)


def six_track(row):
    # Orthrus input: one-hot ACGT, codon starts of the CDS and exon ends
    track = np.zeros((len(row.seq), 6), dtype=np.float32)
    index = ONE_HOT[np.frombuffer(row.seq.encode(), dtype=np.uint8)]
    track[np.flatnonzero(index < 4), index[index < 4]] = 1.0
    track[np.arange(row.cds_start, row.cds_end_nostop, 3), 4] = 1.0
    track[np.cumsum([int(size) for size in row.exon_lens.split(",")]) - 1, 5] = 1.0
    return track


def embed(model, items):
    import torch  # imported in the GPU steps only

    # mean-pooled embeddings, longest first, in batches of MAX_TOKENS
    items = sorted(items, key=lambda item: -len(item[1]))
    embeddings, start = {}, 0
    while start < len(items):
        size = max(1, min(MAX_BATCH, MAX_TOKENS // len(items[start][1])))
        keys, tracks = zip(*items[start : start + size], strict=True)
        batch = np.zeros((len(tracks), len(tracks[0]), 6), dtype=np.float32)
        for row, track in enumerate(tracks):
            batch[row, : len(track)] = track
        lengths = torch.tensor([len(track) for track in tracks], device="cuda")
        with torch.no_grad():
            output = model.representation(torch.from_numpy(batch).cuda(), lengths, channel_last=True)
        embeddings.update(zip(keys, output.float().cpu().numpy(), strict=True))
        start += len(keys)
    return embeddings


def run_orthrus(species, staged, out):
    import torch  # imported in the GPU steps only
    from transformers import AutoModel

    table = pd.read_parquet(staged)
    model = AutoModel.from_pretrained(MODEL, trust_remote_code=True).to("cuda").eval()
    embeddings = embed(model, [(row.tx_id, six_track(row)) for row in table.itertuples()])
    np.savez_compressed(
        Path(out) / f"orthrus_{species}.npz",
        gene_id=table.gene_id.to_numpy(str),
        tx_id=table.tx_id.to_numpy(str),
        emb=np.stack([embeddings[tx] for tx in table.tx_id]).astype(np.float32),
        software=np.array([MODEL, torch.__version__, torch.cuda.get_device_name(0)]),
    )


def run_ribonn(species, staged, out, ribonn_dir):
    staged, out = Path(staged).resolve(), Path(out).resolve()
    sys.path.insert(0, str(ribonn_dir))
    os.chdir(ribonn_dir)
    from src.predict import predict_using_nested_cross_validation_models

    # transcripts RiboNN can score: complete CDS, both UTRs, unambiguous bases
    table = pd.read_parquet(staged)
    table = table[table.start_ATG & table.stop_ok & table.frame_ok]
    table = table[(table.utr5_len > 0) & (table.utr3_len > 0)]
    table = table[table.seq.str.fullmatch("[ACGT]+")]
    parts = zip(table.seq, table.cds_start, table.cds_end, strict=True)
    regions = [(seq[:start], seq[start:end], seq[end:]) for seq, start, end in parts]
    inputs = pd.DataFrame(regions, columns=["utr5_sequence", "cds_sequence", "utr3_sequence"])
    inputs.insert(0, "tx_id", table.tx_id.to_numpy())
    inputs.to_csv(out / f"ribonn_input_{species}.txt", sep="\t", index=False)
    # pretrained nested cross-validation models, five per fold, averaged
    runs = pd.read_csv(f"models/{species}/runs.csv")
    path = str(out / f"ribonn_input_{species}.txt")
    predictions = predict_using_nested_cross_validation_models(path, species, runs, 5, batch_size=64, num_workers=4)
    columns = [name for name in predictions.columns if name.startswith("predicted_")]
    mean = predictions.groupby("tx_id", as_index=False)[columns].mean()
    mean["mean_predicted_TE"] = mean[columns].mean(axis=1)
    mean.to_csv(out / f"ribonn_{species}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    step, species, *paths = sys.argv[1:]
    if step == "inputs":
        transcripts, exons, sequences, out = paths
        tables = stage(*(pd.read_parquet(path) for path in (transcripts, exons, sequences)))
        tables.to_parquet(Path(out) / f"transcripts_{species}.parquet", index=False)
    else:
        {"orthrus": run_orthrus, "ribonn": run_ribonn}[step](species, *paths)
