# ORCO 0.9 Demo

ORCO (Open Residue Classification Oracle) is a lightweight prototype for assigning amino acid residue types based on sparse NMR chemical shift data.

This `orco-0.9-demo` release is a minimal, fully functional scaffold that demonstrates:
- Cascaded confidence-based residue prediction using CA, CB, CG1, CG2, and CE
- Probability-based classification with uncertainty thresholds
- Export in ASAP format

## Included Files
- `main.py` – Core logic for running predictions and exporting results
- `orco_input.csv` – Example spin system input file
- `orco_output.txt` – Output file formatted for ASAP
- `res_seq.txt` – Example residue sequence file
- `stats.csv` – Average chemical shift statistics for each amino acid
- `requirements.txt` – Python dependencies
- `run_in_colab.ipynb` – Optional notebook for running ORCO in Google Colab

## Usage
```bash
python main.py
```

## Notes
- Uncertain predictions list multiple residue types until the confidence threshold is reached.
- This version is designed for educational and demonstration purposes.

