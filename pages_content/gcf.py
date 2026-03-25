# pages_content/gcf.py

import streamlit as st
import pandas as pd
import matplotlib.colors as mcolors
import json
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.renderers.default = "browser"
import altair as alt

@st.cache_data
def load_data():
    # Load data
    virgo2_inventory = pd.read_csv("data/MAG_inventory.csv")
    regions = pd.read_csv("data/04_regions.csv")    
    virgo2_taxakey = pd.read_csv("data/VIRGO2_taxaKey_modif.csv")
    # COVERAGE = pd.read_csv("data/coverage.csv")
    BGC_table = pd.read_csv(f"data/01_BGC_table.csv")
    return virgo2_inventory, BGC_table, virgo2_taxakey, regions

virgo2_inventory, BGC_table, virgo2_taxakey, regions = load_data()

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

with open(f"data/accumulation_models.json", "r") as f:
    accumulation_models = json.load(f) 

def page():

    st.header("MAGs selection", divider="grey")

    st.markdown(f"""
    <div class="justified-text">

    BGC detection is analyzed in relation to MAG features. By default, only MAGs with coverage above 10X and completeness greater than 80% are included on this page.

    You can adjust the coverage and completeness thresholds here if needed.
                
    </div>
    """, unsafe_allow_html=True)

    virgo2_metadata = virgo2_inventory[(virgo2_inventory["Coverage"] > 10) & (virgo2_inventory["Complete"] > 80)]

    col1, col2 = st.columns(2)
    with col1:
        min_coverage = st.number_input(
            "Insert a number", value=10, placeholder="Default is 10"
        )
    with col2:
        min_completeness = st.number_input(
            "Insert a number", value=80, placeholder="Default is 80"
        )

    st.write("MAGs with a minimum coverage of", min_coverage, " and a completeness greather than", min_completeness, "% will be considered")

    virgo2_metadata = virgo2_inventory[(virgo2_inventory["Coverage"] > min_coverage) & (virgo2_inventory["Complete"] > min_completeness)]

    st.header("Gene Cluster Families - all", divider="grey")

    # ----- Merge full metadata into BGC_table
    merged = BGC_table.merge(virgo2_metadata[["MAG", "FinalTaxonomy"]])

    singletons = regions[~regions["GBK"].isin(BGC_table["GBK"])]
    singletons = singletons[["MAG", "contig_id", "GBK", "type", "FinalTaxonomy"]].rename(
        columns={"contig_id": "Description", "type": "Class"}
    )
    singletons["Category"] = "singleton"
    singletons["Family"] = [f"singl. {i+1}" for i in range(singletons.shape[0])]
    singletons = singletons[merged.columns]

    merged = pd.concat([merged, singletons], axis=0)

    # ----- Get top 20 families by total count
    n_max_gcf = merged["Family"].nunique()

    # ----- Add FAM/CLASS/TAXA selector

    st.markdown(f"""
    <div class="justified-text">

    By default, only the number of GCFs displayed on the x-axis is controlled. The plot shows the top 50 GCFs ranked by their total number of BGCs.
    """, unsafe_allow_html=True)

    N_GCF = st.slider(f"Number of GCF (over {n_max_gcf}):", min_value=1, max_value=n_max_gcf, value=50)

    st.markdown(f"""
    <div class="justified-text">
                
    You can refine the selection using additional filters:
    # - <b>Number of GCF</b>: Number of GCFs displayed on the x-axis
    - <b>GCF ID</b>: Display specific GCFs (one or multiple) by selecting them directly.
    - <b>BGC Class</b>: Focus on GCFs that contain specific BGC types.
    - <b>Taxa</b>: Explore GCFs associated with a given taxon or a list of taxa.
    
    Note: Only one filter can be applied at a time; filters are not cumulative.
    </div>
    """, unsafe_allow_html=True)


    col1, col2, col3 = st.columns(3)

    with col1:
        FAM_TO_GREP = st.multiselect("Search by GCF ID", options=sorted(merged["Family"].value_counts().index.to_list()), default=None)
    with col2:
        CLASS_TO_GREP = st.multiselect("Search by BGC class", options=sorted(merged["Class"].value_counts().index.to_list()), default=None)
    with col3:
        TAXA_TO_GREP = st.multiselect("Search by taxa", options=sorted(merged["FinalTaxonomy"].value_counts().index.to_list()), default=None)

    top_families = set()


    if FAM_TO_GREP:
        merged_filter_fam = merged[merged["Family"].isin(FAM_TO_GREP)]
        top_families.update(merged_filter_fam["Family"].value_counts().index)

    if CLASS_TO_GREP:
        merged_filter_class = merged[merged["Class"].isin(CLASS_TO_GREP)]
        top_families.update(merged_filter_class["Family"].value_counts().index)

    if TAXA_TO_GREP:
        merged_filter_taxa = merged[merged["FinalTaxonomy"].isin(TAXA_TO_GREP)]
        top_families.update(merged_filter_taxa["Family"].value_counts().index)

    # Fallback if nothing selected
    if not top_families:
        top_families = set(merged["Family"].value_counts().head(N_GCF).index)


    # if FAM_TO_GREP :
    #     top_families = merged["Family"].value_counts().index
    # else :
    #     top_families = merged["Family"].value_counts().head(N_GCF).index

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

    # ----- Define color palettes and labels for taxa
    taxa_color = virgo2_taxakey.set_index("Taxa")["Color"].to_dict()
    palette_tax = [taxa_color.get(taxon, "#8c8c8c") for taxon in pivot_tax.columns]

    # ----- GCF count for hover template
    gcf_count = {}
    for fam in pivot_tax.index:
        tmp = pivot_tax.loc[fam]
        tmp = tmp.sum()
        gcf_count[fam] = tmp

    if FAM_TO_GREP :
        pivot_tax = pivot_tax.loc[FAM_TO_GREP]
        pivot_class = pivot_class.loc[FAM_TO_GREP]

    # ----- Create subplot layout with 2 rows, 1 column
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True, shared_yaxes=True,
        # subplot_titles=("Taxa", "BGC product")
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
                showlegend=False,
                text=[int(gcf_count[fam]) for fam in pivot_tax.index],
                textposition="none",
                hovertemplate = f"Taxa name: {col}<br>" + "Taxa counts: %{y}<br>" + "Total: %{text}<extra></extra>",
                ),
            row=1, col=1
        )
        fig.update_layout(hoverlabel=dict(bgcolor="white",font_size=16, font_color="black"))

    # ----- Second plot: Class
    for col, color in zip(pivot_class.columns, palette_class):
        fig.add_trace(
            go.Bar(
                x=pivot_class.index,
                y=pivot_class[col],
                name=col,
                marker_color=color,
                legendgroup="class",
                showlegend=False,
                text=[int(gcf_count[fam]) for fam in pivot_tax.index],
                textposition="none", 
                hovertemplate=f"Class: {col}<br>" + "Counts: %{y}<br>" + "Total: %{text}<extra></extra>"
            ),
            row=2, col=1
        )
        fig.update_layout(hoverlabel=dict(bgcolor="white",font_size=16, font_color="black"))

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
    fig.update_yaxes(title_text="Number of BGC colored by MAG taxonomy", row=1, col=1)
    fig.update_yaxes(title_text="Number of BGC colored by BGC type", row=2, col=1)
    # fig.update_xaxes(title_text="", row=2, col=1)

    # Rotate xticks
    fig.update_xaxes(tickangle=-45)

    st.plotly_chart(fig, use_container_width=True)

    st.header("Gene Cluster Families - per taxa", divider="grey")

    st.markdown(f"""
    <div class="justified-text">

    This taxa-level analysis shows two complementary views:  
    - The left plot displays accumulation curves of GCFs across species.  
    - The right plot highlights the GCFs grouping BGCs specific to the selected taxa.  

    Double-click an item in the legend to filter the plot interactively.
                
    </div>
    """, unsafe_allow_html=True)


    taxa_selection = st.selectbox("Taxa selection", sorted(accumulation_models.keys()),index=sorted(accumulation_models.keys()).index("Lactobacillus_crispatus"))

    filtered_df = merged[merged["FinalTaxonomy"] == taxa_selection]

    singletons = regions[~regions["GBK"].isin(BGC_table["GBK"])]
    singletons = singletons[["MAG", "contig_id", "GBK", "type", "FinalTaxonomy"]].rename(
        columns={"contig_id": "Description", "type": "Class"}
    )
    singletons["Category"] = "singleton"
    singletons["Family"] = [f"singl. {i+1}" for i in range(singletons.shape[0])]
    singletons = singletons[merged.columns]
    singletons_filtered = singletons[singletons["FinalTaxonomy"] == taxa_selection].copy()

    BGC_PRODUCT_ORDER = ["RiPP-like", "LAP", "NRPS", "lanthipeptide-class-i", "lanthipeptide-class-ii",
                        "lanthipeptide-class-iii", "lanthipeptide-class-iv", "lanthipeptide-class-v",
                        "T3PKS", "resorcinol", "arylpolyene", "RRE-containing"]

    col1, col2 = st.columns(spec=[0.3,0.7])

    with col1:
        st.subheader("Accumulation model")
        curve_df = pd.DataFrame(accumulation_models[taxa_selection])
        curve_df["Taxa"] = taxa_selection
        palette_tax = [taxa_color.get(taxon, "#8c8c8c") for taxon in [taxa_selection]]

        fig1 = px.line(
            curve_df,
            x="Sites",
            y="Richness",
            color="Taxa",
            color_discrete_sequence=palette_tax
        )
        fig1.update_layout(
            xaxis_title="Number of MAGs",
            yaxis_title="Accumulated unique GCF",
            showlegend=False,
            height=400  # set fixed height to align with col2
        )
        st.plotly_chart(fig1, use_container_width=True)

    # with col2:
    #     st.subheader("GCFs composition")
    #     BGC_PRODUCT_ORDER = ["RiPP-like", "LAP", "NRPS", "lanthipeptide-class-i", "lanthipeptide-class-ii",
    #                         "lanthipeptide-class-iii", "lanthipeptide-class-iv", "lanthipeptide-class-v",
    #                         "T3PKS", "resorcinol", "arylpolyene", "RRE-containing"]

    #     filtered_df = merged[merged["FinalTaxonomy"] == taxa_selection]

    #     # singletons = regions[~regions["GBK"].isin(BGC_table["GBK"])]
    #     # singletons = singletons[["MAG", "contig_id", "GBK", "type", "FinalTaxonomy"]].rename(
    #     #     columns={"contig_id": "Description", "type": "Class"}
    #     # )
    #     # singletons["Category"] = "singleton"
    #     # singletons["Family"] = [f"singl. {i+1}" for i in range(singletons.shape[0])]
    #     # singletons = singletons[merged.columns]
    #     singletons_filtered = singletons[singletons["FinalTaxonomy"] == taxa_selection].copy()

    #     combined_df = pd.concat([filtered_df, singletons_filtered], ignore_index=True)
    #     normalized = combined_df.groupby("Family")["Class"].value_counts(normalize=True).rename("proportion").reset_index()
    #     norm_pivot = normalized.pivot(index="Family", columns="Class", values="proportion").fillna(0)

    #     columns_order = [c for c in BGC_PRODUCT_ORDER if c in norm_pivot.columns] + [c for c in norm_pivot.columns if c not in BGC_PRODUCT_ORDER]

    #     gcf_mask = norm_pivot.index.str.startswith("FAM_")
    #     singleton_mask = norm_pivot.index.str.startswith("singl.")
    #     gcf_sorted = norm_pivot[gcf_mask].sort_values(by=columns_order, ascending=False) if gcf_mask.any() else pd.DataFrame(columns=columns_order)
    #     singleton_sorted = norm_pivot[singleton_mask].sort_values(by=columns_order, ascending=False) if singleton_mask.any() else pd.DataFrame(columns=columns_order)

    #     norm_pivot = pd.concat([gcf_sorted, singleton_sorted])
    #     norm_pivot = norm_pivot.reindex(columns=columns_order).fillna(0)
    #     # gcf_order = norm_pivot.index.to_list()
    #     gcf_order = combined_df.groupby("Family")["Class"].value_counts().rename("proportion").reset_index().pivot(index="Family", columns="Class", values="proportion").fillna(0).T.sum().sort_values(ascending=False).index.to_list()
    #     gcf_order = [i for i in gcf_order if "singl" not in i] + [j for j in gcf_order if "singl" in j]
    #     norm_pivot.index = norm_pivot.index.astype(str)

    #     df_to_plot = norm_pivot.reset_index()
    #     if "index" in df_to_plot.columns:
    #         df_to_plot = df_to_plot.rename(columns={"index":"Family"})

    #     df_to_plot = df_to_plot.melt(
    #         id_vars="Family",
    #         value_vars=norm_pivot.columns,
    #         var_name="Class",
    #         value_name="Proportion"
    #     )

    #     df_to_plot["Family"] = pd.Categorical(df_to_plot["Family"], categories=gcf_order, ordered=True)
    #     df_to_plot = df_to_plot.sort_values("Family").reset_index(drop=True)

    #     # Handle "mix" types
    #     df_to_plot["Class"] = df_to_plot["Class"].apply(lambda x: "mix" if ("," in x or "." in x) else x)

    #     # Update class_colors: assign "#8c8c8c" to "mix"
    #     color_mapping_type_filtered = {k: v for k, v in color_mapping_type.items() if k in norm_pivot.columns}
    #     class_colors = {k: rgb_to_hex(v) if 'rgb' in v else v for k, v in color_mapping_type_filtered.items()}

    #     if "mix" in df_to_plot["Class"].unique():
    #         class_colors["mix"] = "#8c8c8c"

    #     # Build Altair color scale
    #     color_scale = alt.Scale(domain=list(class_colors.keys()), range=list(class_colors.values()))

    #     # Create stacked barplot
    #     fig2 = alt.Chart(df_to_plot).mark_bar().encode(
    #         x=alt.X("Family:N", sort=gcf_order, title=""),
    #         y=alt.Y("Proportion:Q", stack="normalize"),
    #         color=alt.Color("Class:N", scale=color_scale, legend=alt.Legend(orient="top", title="Class")),
    #         tooltip=["Family", "Class", "Proportion"]
    #     ).properties(height=400)

    #     st.altair_chart(fig2, use_container_width=True)
    with col2:
        combined_df = pd.concat([filtered_df, singletons_filtered], ignore_index=True)
        normalized = combined_df.groupby("Family")["Class"].value_counts().rename("proportion").reset_index() #normalize=True
        norm_pivot = normalized.pivot(index="Family", columns="Class", values="proportion").fillna(0)

        columns_order = [c for c in BGC_PRODUCT_ORDER if c in norm_pivot.columns] + [c for c in norm_pivot.columns if c not in BGC_PRODUCT_ORDER]

        gcf_mask = norm_pivot.index.str.startswith("FAM_")
        singleton_mask = norm_pivot.index.str.startswith("singl.")
        gcf_sorted = norm_pivot[gcf_mask].sort_values(by=columns_order, ascending=False) if gcf_mask.any() else pd.DataFrame(columns=columns_order)
        singleton_sorted = norm_pivot[singleton_mask].sort_values(by=columns_order, ascending=False) if singleton_mask.any() else pd.DataFrame(columns=columns_order)

        norm_pivot = pd.concat([gcf_sorted, singleton_sorted])
        norm_pivot = norm_pivot.reindex(columns=columns_order).fillna(0)
        gcf_order = norm_pivot.T.sum().sort_values(ascending=False).index.to_list()
        norm_pivot.index = norm_pivot.index.astype(str)

        df_to_plot = norm_pivot.reset_index()
        if "index" in df_to_plot.columns:
            df_to_plot = df_to_plot.rename(columns={"index":"Family"})

        df_to_plot = df_to_plot.melt(
            id_vars="Family",
            value_vars=norm_pivot.columns,
            var_name="Class",
            value_name="Proportion"
        )

        df_to_plot["Family"] = pd.Categorical(df_to_plot["Family"], categories=gcf_order, ordered=True)
        df_to_plot = df_to_plot.sort_values("Family").reset_index(drop=True)

        # Handle "mix" types
        df_to_plot["Class"] = df_to_plot["Class"].apply(lambda x: "mix" if ("," in x or "." in x) else x)

        # Update class_colors: assign "#8c8c8c" to "mix"
        color_mapping_type_filtered = {k: v for k, v in color_mapping_type.items() if k in norm_pivot.columns}
        class_colors = {k: rgb_to_hex(v) if 'rgb' in v else v for k, v in color_mapping_type_filtered.items()}

        if "mix" in df_to_plot["Class"].unique():
            class_colors["mix"] = "#8c8c8c"

        # Build Altair color scale
        # color_scale = alt.Scale(domain=list(class_colors.keys()), range=list(class_colors.values()))

        # # Create stacked barplot
        # fig3 = alt.Chart(df_to_plot).mark_bar().encode(
        #     x=alt.X("Family:N", sort=gcf_order, title=""),
        #     y=alt.Y("Proportion:Q"), #stack="normalize"
        #     color=alt.Color("Class:N", scale=color_scale, legend=alt.Legend(orient="top", title="Class")),
        #     tooltip=["Family", "Class", "Proportion"]
        # ).properties(height=400)
        # st.altair_chart(fig3, use_container_width=True)


        color_scale = dict(zip(class_colors.keys(), class_colors.values()))

        # Create figure
        fig3 = go.Figure()

        # Loop through each Class to create stacked bars
        for cls in df_to_plot["Class"].unique():
            tmp = df_to_plot[df_to_plot["Class"] == cls]

            fig3.add_trace(
                go.Bar(
                    x=tmp["Family"],
                    y=tmp["Proportion"],
                    name=cls,
                    marker_color=color_scale.get(cls),  # safe lookup from your color mapping
                    hovertemplate="Family: %{x}<br>Class: " + cls + "<br>Proportion: %{y}<extra></extra>"
                )
            )

        # Update layout for stacked bars
        fig3.update_layout(
            barmode="stack",
            height=400,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="center",
                x=0.5,
                title="Class"
            ),
            xaxis=dict(title=""),
            yaxis=dict(title="Number of BGCs")
        )

        # Show in Streamlit
        st.plotly_chart(fig3, use_container_width=True)
