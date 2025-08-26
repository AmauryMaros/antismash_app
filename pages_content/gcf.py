import streamlit as st
import pandas as pd
import matplotlib.colors as mcolors
import json
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

@st.cache_data
def load_data():
    # Load data
    virgo2_metadata = pd.read_csv("data/03_virgo2_metadata.csv")
    # region_summary_original = pd.read_csv("data/04_regions.csv")    
    virgo2_taxakey = pd.read_csv("data/06_VIRGO2_taxaKey.csv")
    # COVERAGE = pd.read_csv("data/coverage.csv")
    BGC_table = pd.read_csv(f"data/01_BGC_table.csv")
    return virgo2_metadata, BGC_table, virgo2_taxakey

virgo2_metadata, BGC_table, virgo2_taxakey = load_data()


def rgb_to_hex(rgb_string):
    '''
    Convert an RGB color string to a hexadecimal color code.
    Parameters:
        rgb_string (str): A string representing an RGB color, e.g., 'rgb(255, 0, 128)'.
    Returns:
        str: The corresponding hexadecimal color code, e.g., '#ff0080'.
    '''
    if rgb_string.startswith('rgb'):
        return mcolors.rgb2hex([int(x)/255 for x in rgb_string.strip('rgb()').split(',')])
    else:
        return rgb_string 

with open(f"data/color_mapping_type.json", "r") as f:
    color_mapping_type = json.load(f) 
color_mapping_type["lanthipeptide-class-iv"] = "magenta"

def page():

    st.title("Gene Cluster Family")

    # ----- Merge full metadata into BGC_table
    merged = BGC_table.merge(virgo2_metadata[["MAG", "FinalTaxonomy"]])

    # ----- Get top 20 families by total count
    N_GCF = st.slider("Number of GCF:", min_value=1, max_value=merged["Family"].nunique(), value=50)
    top_families = merged["Family"].value_counts().head(N_GCF).index

    # ----- Filter merged data to top families
    filtered_merged = merged[merged["Family"].isin(top_families)]

    # ----- Compute order of families by total count (for consistent bar order)
    family_order = filtered_merged.groupby("Family").size().sort_values(ascending=False).index

    # ----- Group for FinalTaxonomy plot
    counts_tax = filtered_merged.groupby(["Family", "FinalTaxonomy"]).size().reset_index(name="count")
    pivot_tax = counts_tax.pivot_table(index="Family",columns="FinalTaxonomy",values="count",fill_value=0).reindex(family_order)

    # ----- Group for Class plot
    counts_class = filtered_merged.groupby(["Family", "Class"]).size().reset_index(name="count")
    pivot_class = counts_class.pivot_table(index="Family",columns="Class",values="count",fill_value=0).reindex(family_order)

    # ----- Define colors for Class
    class_colors = {k: rgb_to_hex(v) if 'rgb' in v else v for k, v in color_mapping_type.items()}
    palette_class = [class_colors.get(cls, "#8c8c8c") for cls in pivot_class.columns]

    # # ----- Define color palettes and labels for taxa
    taxa_color = virgo2_taxakey.set_index("Taxa")["Color"].to_dict()
    palette_tax = [taxa_color.get(taxon, "#8c8c8c") for taxon in pivot_tax.columns]

    # Create subplot layout with 2 rows, 1 column
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True, shared_yaxes=True,
        subplot_titles=("Taxa", "BGC product")
    )

    # ----- First plot: FinalTaxonomy
    for col, color in zip(pivot_tax.columns, palette_tax):
        fig.add_trace(
            go.Bar(
                x=pivot_tax.index,
                y=pivot_tax[col],
                name=col,
                marker_color=color,
                legendgroup="taxa",
                showlegend=False  # we'll control legend visibility
            ),
            row=1, col=1
        )

    # ----- Second plot: Class
    for col, color in zip(pivot_class.columns, palette_class):
        fig.add_trace(
            go.Bar(
                x=pivot_class.index,
                y=pivot_class[col],
                name=col,
                marker_color=color,
                legendgroup="class",
                showlegend=False
            ),
            row=2, col=1
        )

    # Layout adjustments
    fig.update_layout(
        barmode="stack",
        height=800, width=1000,
        legend=dict(
            title="Taxa / BGC product",
            orientation="h",
            yanchor="bottom", y=1.05,
            xanchor="center", x=0.5,
            bordercolor="black", borderwidth=1
        ),
    )

    # Axis labels
    fig.update_yaxes(title_text="Number of BGC", row=1, col=1)
    fig.update_yaxes(title_text="Number of BGC", row=2, col=1)
    # fig.update_xaxes(title_text="", row=2, col=1)

    # Rotate xticks
    fig.update_xaxes(tickangle=-45)

    st.plotly_chart(fig, use_container_width=True)
