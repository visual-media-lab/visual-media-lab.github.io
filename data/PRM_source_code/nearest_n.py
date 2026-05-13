import pandas as pd
import numpy as np

MIGRATION_XLSX = "State_to_State_Migration_Table_2022.xlsx"
VARIANT_CSV = "Moran_rev.csv"

N_PERM = 999
SEED = 12345

STATE_ABBR = {
    "Alabama":"AL","Alaska":"AK","Arizona":"AZ","Arkansas":"AR",
    "California":"CA","Colorado":"CO","Connecticut":"CT","Delaware":"DE",
    "District of Columbia":"DC","Florida":"FL","Georgia":"GA","Hawaii":"HI",
    "Idaho":"ID","Illinois":"IL","Indiana":"IN","Iowa":"IA","Kansas":"KS",
    "Kentucky":"KY","Louisiana":"LA","Maine":"ME","Maryland":"MD",
    "Massachusetts":"MA","Michigan":"MI","Minnesota":"MN","Mississippi":"MS",
    "Missouri":"MO","Montana":"MT","Nebraska":"NE","Nevada":"NV",
    "New Hampshire":"NH","New Jersey":"NJ","New Mexico":"NM","New York":"NY",
    "North Carolina":"NC","North Dakota":"ND","Ohio":"OH","Oklahoma":"OK",
    "Oregon":"OR","Pennsylvania":"PA","Rhode Island":"RI",
    "South Carolina":"SC","South Dakota":"SD","Tennessee":"TN","Texas":"TX",
    "Utah":"UT","Vermont":"VT","Virginia":"VA","Washington":"WA",
    "West Virginia":"WV","Wisconsin":"WI","Wyoming":"WY",
    "Puerto Rico":"PR",
    "U.S. Island Area":"VI",
    "Virgin Islands":"VI",
    "United States Virgin Islands":"VI"
}

VALID_STATE_NAMES = set(STATE_ABBR.keys())

def clean_name(x):
    if pd.isna(x):
        return None

    s = str(x).replace("\xa0", " ").strip()
    s = " ".join(s.split())

    if s.startswith("United States"):
        return None

    return s

def read_migration_matrix(path):
    raw = pd.read_excel(
        path,
        sheet_name="Table",
        header=None,
        engine="openpyxl"
    )

    row_label_cols = [0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 121]
    header_rows = [6, 45]

    origin_cols_dict = {}

    for header_row in header_rows:
        for c in range(raw.shape[1]):
            state_name = clean_name(raw.iat[header_row, c])
            subheader = clean_name(raw.iat[header_row + 1, c])

            if state_name in VALID_STATE_NAMES and subheader == "Estimate":
                origin_cols_dict[c] = STATE_ABBR[state_name]

    origin_cols = list(origin_cols_dict.items())

    flows = []

    for r in range(raw.shape[0]):
        destination = None

        for c in row_label_cols:
            if c < raw.shape[1]:
                nm = clean_name(raw.iat[r, c])

                if nm in VALID_STATE_NAMES:
                    destination = STATE_ABBR[nm]
                    break

        if destination is None:
            continue

        for c, origin in origin_cols:
            val = raw.iat[r, c]

            if pd.isna(val):
                continue

            try:
                val = float(val)
            except Exception:
                continue

            if val <= 0 or origin == destination:
                continue

            flows.append((origin, destination, val))

    flow_df = pd.DataFrame(
        flows,
        columns=["origin", "destination", "flow"]
    )

    flow_df = flow_df.groupby(
        ["origin", "destination"],
        as_index=False
    )["flow"].sum()

    return flow_df

def read_variant_counts(path):
    df = pd.read_csv(
        path,
        header=None,
        names=["variant", "state", "count"],
        dtype={
            "variant": str,
            "state": str,
            "count": float
        }
    )

    df = df.dropna(subset=["variant", "state", "count"])
    df["state"] = df["state"].str.strip().str.upper()

    return df.groupby(
        ["variant", "state"],
        as_index=False
    )["count"].sum()

def make_mobility_distance_matrix(flow_df, states):
    idx = {s: i for i, s in enumerate(states)}
    n_states = len(states)

    # Directed[i, j] = flow from state i to state j
    Directed = np.zeros((n_states, n_states), dtype=float)

    for _, row in flow_df.iterrows():
        origin = row["origin"]
        destination = row["destination"]
        flow = float(row["flow"])

        if origin in idx and destination in idx and origin != destination:
            i = idx[origin]
            j = idx[destination]
            Directed[i, j] += flow

    np.fill_diagonal(Directed, 0)

    # Symmetric flow matrix:
    # F[i, j] = Directed[i, j] + Directed[j, i]
    F = Directed + Directed.T
    np.fill_diagonal(F, 0)

    total_flow = F.sum()

    if total_flow == 0:
        raise ValueError("Total migration flow is zero.")

    # P[i, j] is the globally normalized symmetric migration flow.
    P = F / total_flow

    positive = P[P > 0]

    if len(positive) == 0:
        raise ValueError("No positive migration probabilities.")

    min_positive = positive.min()

    # Replace zero probabilities with a small positive value
    # so that -log(P) is finite.
    P_safe = np.where(
        P > 0,
        P,
        min_positive * 0.5
    )

    # Mobility distance:
    # larger flow -> larger P -> smaller distance
    D = -np.log(P_safe)

    # Same-state distance is defined as zero.
    np.fill_diagonal(D, 0)

    return Directed, F, P, D

def expand_samples(sub):
    samples = []

    for _, row in sub.iterrows():
        state = row["state"]
        count = int(round(float(row["count"])))

        for _ in range(count):
            samples.append(state)

    return samples

def mean_nearest_neighbor_distance(samples, D, states):
    idx = {s: i for i, s in enumerate(states)}

    sample_indices = [
        idx[s]
        for s in samples
        if s in idx
    ]

    if len(sample_indices) < 2:
        return np.nan

    nn_distances = []

    for k, i in enumerate(sample_indices):
        distances = []

        for l, j in enumerate(sample_indices):
            if k == l:
                continue

            distances.append(D[i, j])

        nn_distances.append(min(distances))

    return float(np.mean(nn_distances))

def permutation_p(
    samples,
    D,
    states,
    observed,
    n_perm=999,
    seed=12345
):
    if np.isnan(observed):
        return np.nan

    rng = np.random.default_rng(seed)
    n = len(samples)

    sims = []

    for _ in range(n_perm):
        random_samples = rng.choice(
            states,
            size=n,
            replace=True
        )

        sim = mean_nearest_neighbor_distance(
            random_samples,
            D,
            states
        )

        sims.append(sim)

    sims = np.array(sims)

    # Smaller mean nearest-neighbor distance means stronger clustering.
    p_lower = (
        np.sum(sims <= observed) + 1
    ) / (n_perm + 1)

    return p_lower

# ==========================================================
# Read data
# ==========================================================

flow_df = read_migration_matrix(MIGRATION_XLSX)
variant_df = read_variant_counts(VARIANT_CSV)

# ==========================================================
# Jurisdictions to exclude
# ==========================================================

EXCLUDE = {"GU", "AS", "MP"}

states = sorted(
    set(flow_df["origin"])
    .union(flow_df["destination"])
    .union(variant_df["state"])
)

states = [
    s for s in states
    if s not in EXCLUDE
]

flow_df = flow_df[
    flow_df["origin"].isin(states)
    &
    flow_df["destination"].isin(states)
]

variant_df = variant_df[
    variant_df["state"].isin(states)
]

# ==========================================================
# Build matrices
# ==========================================================

Directed, F, P, D = make_mobility_distance_matrix(
    flow_df,
    states
)

# ==========================================================
# Debug outputs
# ==========================================================

directed_flow_matrix_df = pd.DataFrame(
    Directed,
    index=states,
    columns=states
)

directed_flow_matrix_df.to_csv(
    "debug_directed_migration_flow_matrix.csv"
)

symmetric_flow_matrix_df = pd.DataFrame(
    F,
    index=states,
    columns=states
)

symmetric_flow_matrix_df.to_csv(
    "debug_symmetric_migration_flow_matrix.csv"
)

probability_matrix_df = pd.DataFrame(
    P,
    index=states,
    columns=states
)

probability_matrix_df.to_csv(
    "debug_migration_probability_matrix.csv"
)

distance_matrix_df = pd.DataFrame(
    D,
    index=states,
    columns=states
)

distance_matrix_df.to_csv(
    "debug_mobility_distance_matrix.csv"
)

# ==========================================================
# Main analysis
# ==========================================================

results = []

for variant, sub in variant_df.groupby("variant"):
    samples = expand_samples(sub)

    observed = mean_nearest_neighbor_distance(
        samples,
        D,
        states
    )

    p_value = permutation_p(
        samples,
        D,
        states,
        observed,
        n_perm=N_PERM,
        seed=SEED
    )

    if p_value <= 1 / (N_PERM + 1):
        p_text = f"< {1/(N_PERM+1):.3f}"
    else:
        p_text = f"{p_value:.3f}"

    results.append({
        "variant": variant,
        "sample_count": len(samples),
        "occupied_states": sub["state"].nunique(),
        "mean_mobility_NN_distance": observed,
        "p_value_clustered": p_text
    })

out = pd.DataFrame(results).sort_values("variant")

print(out)

out.to_csv(
    "mobility_nearest_neighbor_results.csv",
    index=False
)

#print("\nDEBUG FILES GENERATED:")
#print("  debug_directed_migration_flow_matrix.csv")
#print("  debug_symmetric_migration_flow_matrix.csv")
#print("  debug_migration_probability_matrix.csv")
#print("  debug_mobility_distance_matrix.csv")
print("  mobility_nearest_neighbor_results.csv")