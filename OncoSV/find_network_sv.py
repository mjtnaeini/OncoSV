import pandas as pd
import numpy as np
import networkx as nx
import community as community_louvain
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pyvis.network import Network
import random

# ─────────────────────────────────────────
# 1.  Build SV graph (shared-read weights)
# ─────────────────────────────────────────
def build_sv_graph(df, min_read_count: int = 2) -> nx.Graph:
    """
    For every row, connect all SV IDs that co-occur in that read pattern.
    Edge weights = Read_Count (can be filtered by min_read_count).
    """
    G = nx.Graph()
    for _, row in df.iterrows():
        svs = str(row["ID"]).split(",")
        for i in range(len(svs)):
            for j in range(i + 1, len(svs)):
                if row["Read_Count"] >= min_read_count:
                    G.add_edge(svs[i], svs[j], weight=row["Read_Count"])
    return G


# ─────────────────────────────────────────
# 2.  Louvain → SV_Clusters
# ─────────────────────────────────────────
def label_sv_clusters(df: pd.DataFrame, G: nx.Graph):
    """
    Add a column Cluster_number = smallest Louvain community ID for the SV IDs in row.ID.
    """
    partition = community_louvain.best_partition(G)
    df = df.copy()
    df["Cluster_number"] = df["ID"].apply(
        lambda ids: min(partition.get(sv, -1) for sv in ids.split(","))
    )
    return df, partition

# ─────────────────────────────────────────
# 3.  AF summary
# ─────────────────────────────────────────
def add_af_stats(df: pd.DataFrame):
    af_lists = (
        df["AF"]
        .astype(str)
        .str.split(";", expand=False)
        .apply(lambda xs: [float(x) for x in xs if x])
    )
    df = df.copy()
    df["mean_AF"] = af_lists.apply(np.mean)
    df["sd_AF"] = af_lists.apply(np.std)
    return df

# ─────────────────────────────────────────
# 4.  Clone assignment  (contiguous subseq rule)
# ─────────────────────────────────────────
def _is_contig_subseq(short: list, long: list) -> bool:
    """short appears contiguously and in order inside long."""
    n, m = len(short), len(long)
    if n > m:
        return False
    for i in range(m - n + 1):
        if long[i : i + n] == short:
            return True
    return False


def assign_clone_ids_per_cluster(df: pd.DataFrame):
    df = df.copy()
    df["ID_list"] = df["ID"].str.split(",")
    df["ID_len"] = df["ID_list"].apply(len)
    df["Clone_ID"] = None

    for clust in sorted(df["Cluster_number"].unique()):
        idx_mask = df["Cluster_number"] == clust
        ordered = (
            df.loc[idx_mask]
            .sort_values(["ID_len", "mean_AF"], ascending=[False, False])
            .index
        )

        bases = []          # list of (label, id_list)
        counter = 1

        for idx in ordered:
            ids = df.at[idx, "ID_list"]
            label = None

            for base_lbl, base_ids in bases:
                if (_is_contig_subseq(ids, base_ids) or
                    _is_contig_subseq(base_ids, ids)):
                    label = (
                        base_lbl if ids == base_ids else f"nested {base_lbl}"
                    )
                    break

            if label is None:
                label = f"SC{clust}.{counter}"
                bases.append((label, ids))
                counter += 1

            df.at[idx, "Clone_ID"] = label

    return df.drop(columns=["ID_list", "ID_len"])

# ─────────────────────────────────────────
# 5.  Build & cluster only rows meeting the read-count threshold
# ─────────────────────────────────────────
def process_with_modularity(df_in: pd.DataFrame, min_read_count: int = 2):
    df = df_in.copy()

    # Filter out weak-support patterns *before* building the graph
    df = df[df["Read_Count"] >= min_read_count].copy()

    # Build shared-read graph
    G = build_sv_graph(df, min_read_count)

    # Louvain → partition  (0-based)
    partition0 = community_louvain.best_partition(G)

    # Shift cluster IDs so smallest is 1
    partition = {sv: cid + 1 for sv, cid in partition0.items()}

    # Annotate rows
    df["Cluster_number"] = df["ID"].apply(
        lambda ids: min(partition.get(sv, -1) for sv in ids.split(","))
    )

    # (no –1 values remain because every SV in df is in the graph)
    # … then AF stats & clone logic as before …
    df = add_af_stats(df)
    df = assign_clone_ids_per_cluster(df)

    # reorder & return
    df = df[
        [
            "ID","CHROM","CHROM2","POS","END","POS_BKPT","END_BKPT",
            "Read_Count","SV_Count","AF","mean_AF","sd_AF",
            "Sample","CSV_Type","final_combination","any_overlapping_sv",
            "Cluster_number","Clone_ID",
        ]
    ]
    return df, G, partition

def make_all_clusters_html(sv_graph,
                           partition,           # dict {SV_id: SV_Cluster}
                           clone_map,           # dict {SV_id: Clone_ID}
                           df_metadata=None,    # DataFrame with ID→CHROM (optional)
                           outfile="all_clusters_interactive.html"):
    """
    Generates one-file interactive vis-network visualisation.

    Parameters
    ----------
    sv_graph   : NetworkX Graph returned by process_with_modularity
    partition  : dict {SV_id: Cluster_number}
    clone_map  : dict {SV_id: Clone_ID string}
    df_metadata: optional DataFrame that has columns ['ID','CHROM'] so we can
                 pick chromosome per SV more robustly.  If None, chromosome is
                 parsed as the second token in SV id 'Tool.SVTYPE.chr'.
    outfile    : HTML file to write
    """
    # ---------------- palettes ----------------
    chrom_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chrom_palette = dict(zip(
        chrom_list,
        plt.get_cmap("tab20").colors * 2      # 24 distinct colours
    ))

    clone_ids = sorted(set(clone_map.values()))
    clone_palette = {
        cl: mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS)[i % 10]]
        for i, cl in enumerate(clone_ids)
    }

    # helper to fetch chromosome for an SV
    id_to_chrom = {}
    if df_metadata is not None:
        id_to_chrom = dict(zip(df_metadata["ID"], df_metadata["CHROM"]))

    def sv_chrom(svid: str):
        if svid in id_to_chrom:
            return id_to_chrom[svid]
        parts = svid.split(".")
        return parts[1] if len(parts) > 1 else "chr?"

    # ---------------- build network ----------------
    net = Network(height="800px", width="100%", bgcolor="#ffffff")
    net.force_atlas_2based(gravity=-30, spring_length=120, damping=0.9)

    for node in sv_graph.nodes():
        cid   = partition[node]
        clone = clone_map.get(node, "?")
        chrom = sv_chrom(node)

        net.add_node(
            node,
            label=node,
            title=f"Cluster {cid}<br>{clone}<br>{chrom}",
            cluster=cid,
            chrom_col=chrom_palette.get(chrom, "#888888"),
            clone_col=clone_palette.get(clone, "#888888"),
            color=chrom_palette.get(chrom, "#888888")   # default view = chromosome
        )

    for u, v in sv_graph.edges():
        net.add_edge(u, v, color="#bbbbbb")

    base_html = net.generate_html()

    # ---------------- cluster checkboxes ----------------
    box_div = "<div style='position:fixed;left:10px;top:10px;" \
              "background:#fff;border:1px solid #ccc;padding:6px;z-index:2;'>"
    for cid in sorted(set(partition.values())):
        box_div += (f"<label style='display:block;font-size:13px;'>"
                    f"<input type='checkbox' checked "
                    f"onchange='toggleCluster({cid},this.checked)'> "
                    f"Cluster {cid}</label>")
    box_div += "</div>"

    # ---------------- palette radio buttons ----------------
    pal_div = """
    <div style='position:fixed;left:10px;bottom:10px;
                background:#fff;border:1px solid #ccc;
                padding:6px;z-index:2;font-size:13px;'>
      Colour by:<br>
      <label><input type='radio' name='pal' value='chrom'
             onchange="setPalette('chrom')" checked> Chromosome</label><br>
      <label><input type='radio' name='pal' value='clone'
             onchange="setPalette('clone')"> Clone ID</label>
    </div>
    """

    # ---------------- inject JS ----------------
    js_code = """
    <script>
      function toggleCluster(cid, show){
        network.body.data.nodes.forEach(function(n){
          if(n.cluster==cid){
            network.body.data.nodes.update({id:n.id, hidden:!show});
          }
        });
        network.body.data.edges.forEach(function(e){
          const n1 = network.body.data.nodes.get(e.from);
          const n2 = network.body.data.nodes.get(e.to);
          network.body.data.edges.update({id:e.id, hidden:(n1.hidden||n2.hidden)});
        });
      }

      function setPalette(mode){
        network.body.data.nodes.forEach(function(n){
          const col = (mode==='chrom') ? n.chrom_col : n.clone_col;
          network.body.data.nodes.update({id:n.id, color:col});
        });
      }
    </script>
    """

    out_html = base_html.replace("</body>",
                                 box_div + pal_div + js_code + "</body>")
    Path(outfile).write_text(out_html)
    print(f"Interactive graph written → {outfile}")

# build interactive HTML
clone_map = {sv:clone
             for ids, clone in zip(final_df["ID"], final_df["Clone_ID"])
             for sv in ids.split(",")}
