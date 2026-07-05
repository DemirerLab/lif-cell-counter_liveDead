# LIF Cell Counter

A Python command-line workflow for counting live and dead *Arabidopsis thaliana* cells from Leica `.lif` confocal microscopy files.

Expected channel convention:

- `C=0`: fluorescently labeled cytosol of live cells
- `C=1`: fluorescent nuclei/puncta of dead cells
- `C=2`: optional overlay/reference channel, not analyzed by default

For each image series, the tool outputs an Excel spreadsheet and QC preview PNGs with numeric annotations.

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/lif-cell-counter.git
cd lif-cell-counter
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e .
```

## Basic usage

```bash
lif-cell-counter path/to/file.lif \
  --out counts.xlsx \
  --preview-dir previews \
  --zip-previews previews.zip
```

## Reference command

```bash
lif-cell-counter c1c2-1.lif \
  --out c1c2-1_counts_REFINED.xlsx \
  --preview-dir c1c2-1_previews_REFINED \
  --zip-previews c1c2-1_previews_REFINED.zip \
  --channel-live 0 \
  --channel-dead 1 \
  --live-bg-sigma 25 \
  --live-smooth-sigma 1.0 \
  --live-sauvola-k 0.22 \
  --live-min-area 300 \
  --live-max-area 30000 \
  --live-hmax-h 2.0 \
  --dead-bg-sigma 10 \
  --dead-smooth-sigma 0.8 \
  --dead-tophat-radius 10 \
  --dead-min-area 30 \
  --dead-max-area 4000 \
  --dead-min-circularity 0.55
```

## Tuning

For more live-cell separation, lower:

```bash
--live-hmax-h 1.0
```

For stricter nuclei-only dead-cell counting, increase:

```bash
--dead-min-circularity 0.70 --dead-min-area 80 --dead-tophat-radius 14
```

## Notes

`.lif`, `.xlsx`, `.zip`, and preview folders are ignored by `.gitignore` so large microscopy files and generated outputs are not accidentally committed.

## License

MIT
