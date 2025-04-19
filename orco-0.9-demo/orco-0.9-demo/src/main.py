import pandas as pd
import numpy as np

# Configurable parameters
TOP_K_SUM_THRESHOLD = 0.9
LOCK_IN_THRESHOLD = 0.90
MISSING_VAL = 1e6
STAGES = [
    ["CA"],
    ["CA", "CB"],
    ["CA", "CB", "CG1"],
    ["CA", "CB", "CG1", "CG2"],
    ["CA", "CB", "CG1", "CG2", "CE"]
]

def load_residue_stats(path):
    df = pd.read_csv(path)
    stats = {}
    for _, row in df.iterrows():
        res = row["Residue"]
        stats[res] = {
            "CA": (row["CA_mean"], 1),
            "CB": (row["CB_mean"], 1),
            "CG1": (row["CG1_mean"], 1),
            "CG2": (row["CG2_mean"], 1),
            "CE": (row["CE_mean"], 1)
        }
    return stats

def load_sequence(path):
    with open(path) as f:
        return list(f.read().replace("\n", "").strip())

def normalize_probs(distances):
    if not distances:
        return {}
    inv = {k: np.exp(-v) for k, v in distances.items()}
    total = sum(inv.values())
    return {k: round(v / total, 5) for k, v in inv.items()}

def oracle_score(spin, features, stats, idx, level, coefa=1.0, coefb=0.7):
    scores = {}
    for res, values in stats.items():
        dist = 0
        valid = 0
        for feat in features:
            obs = spin.get(feat, MISSING_VAL)
            mean, _ = values.get(feat, (MISSING_VAL, 1))
            if obs >= 100000 or mean >= 100000:
                continue
            weight = coefa if feat == "CA" else (coefb if feat == "CB" else 1.0)
            dist += weight * (obs - mean) ** 2
            valid += 1
        if valid > 0:
            scores[res] = dist / valid
    return normalize_probs(scores)

def cascade_predict(spin_df, stats, res_list):
    assigned = {}
    unassigned = spin_df.copy()
    all_probs = {}

    for stage_idx, features in enumerate(STAGES):
        new_unassigned = []
        for idx, spin in unassigned.iterrows():
            spin_dict = spin.to_dict()
            probs = oracle_score(spin_dict, features, stats, idx, f"Oracle{stage_idx+1}")
            if not probs:
                continue
            all_probs[idx] = probs
            top_res = max(probs, key=probs.get)
            top_conf = probs[top_res]
            if top_conf >= LOCK_IN_THRESHOLD:
                assigned[idx] = (top_res, top_conf, stage_idx + 1)
            else:
                new_unassigned.append(spin)
        unassigned = pd.DataFrame(new_unassigned)

    return assigned, all_probs

def output_asap_format(spin_df, assignments, all_probs, output_path):
    output = []
    for idx, row in spin_df.iterrows():
        spin = row.to_dict()
        if idx in assignments:
            residue = assignments[idx][0]
        else:
            probs = all_probs.get(idx, {})
            if probs:
                sorted_probs = sorted(probs.items(), key=lambda x: -x[1])
                residue_list = []
                cumulative = 0.0
                for res, prob in sorted_probs:
                    residue_list.append(res)
                    cumulative += prob
                    if cumulative >= TOP_K_SUM_THRESHOLD:
                        break
                residue = "/".join(residue_list)
            else:
                residue = "UNK"
        line = f"{spin['N']} {spin['CA']} {spin['C']} {spin['CB']} {spin['CG1']} {spin['CG2']} {spin['CE']} 0.3 0.3 0.3 0.3 0.001 0.001 0.001 1 {residue}"
        output.append(line)
    with open(output_path, "w") as f:
        f.write(f"{len(output)} 7\n")
        f.write("\n".join(output))

def main():
    spin_df = pd.read_csv("data/orco_input.csv")
    stats = load_residue_stats("data/stats.csv")
    sequence = load_sequence("data/res_seq.txt")

    assignments, all_probs = cascade_predict(spin_df, stats, sequence)
    output_asap_format(spin_df, assignments, all_probs, "outputs/orco_output.txt")
    print("✅ Finished! Wrote predictions to outputs/orco_output.txt")

if __name__ == "__main__":
    main()
