# pages_content/bgc.py

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Caching the data loading functions to speed up the Streamlit app
@st.cache_data
def load_data():
    # Load data
    virgo2_metadata_all = pd.read_csv("data/MAG_inventory.csv")
    region_summary_original = pd.read_csv("data/regions.csv")    
    virgo2_taxakey = pd.read_csv("data/06_VIRGO2_taxaKey.csv")
    COVERAGE = pd.read_csv("data/coverage.csv")

    return virgo2_metadata_all, region_summary_original, virgo2_taxakey, COVERAGE

# Load the data using the cached function
virgo2_metadata_all, region_summary_original, virgo2_taxakey, COVERAGE = load_data()

# # Filter MAGs
virgo2_inventory = virgo2_metadata_all[(virgo2_metadata_all["Coverage"] > 10) & (virgo2_metadata_all["Complete"] > 80)]

zero_color = "#708090"
positive_color = "#90EE90"

color_map = {1: positive_color, 0: zero_color}
taxa_color = virgo2_taxakey.set_index("Taxa")["Color"].to_dict()

# Data processing
region_df = region_summary_original.copy()
all_mags = virgo2_inventory['MAG'].unique()
mag_w_antismash_result = region_df['MAG'].unique()
mag_no_antismash_result = [i for i in all_mags if i not in mag_w_antismash_result]
antismash_status = pd.concat([
    pd.DataFrame({"MAG": mag_w_antismash_result, "status": 1}),
    pd.DataFrame({"MAG": mag_no_antismash_result, "status": 0})
], axis=0).sort_values("MAG", ascending=True)

# Merge and process data for display
stack_antismash_status = antismash_status.merge(virgo2_inventory[['MAG', 'FinalTaxonomy']], on='MAG', how='left')
status_counts = stack_antismash_status.groupby(['FinalTaxonomy', 'status']).size().unstack(fill_value=0)

# Sort FinalTaxonomy by total count in descending order
status_counts['Total'] = status_counts.sum(axis=1)
status_counts = status_counts.sort_values(by='Total', ascending=False).drop(columns='Total')
status_counts_long = status_counts.reset_index().melt(id_vars='FinalTaxonomy', var_name='status', value_name='count')

# Streamlit page function
def page():
    st.header("Biosynthetic Gene Clusters - all MAGs", divider='grey')

    # N_MAG = virgo2_metadata_all.shape[0]
    virgo2_inventory = virgo2_metadata_all[(virgo2_metadata_all["Coverage"] > 10) & (virgo2_metadata_all["Complete"] > 80)]

    st.markdown(f"""
    <div class="justified-text">

    BGC detection is analyzed in relation to MAG features. By default, only MAGs with coverage above 10X and completeness greater than 80% are included on this page, representing a number of {virgo2_inventory.shape[0]} MAGs (over {virgo2_metadata_all.shape[0]} MAGs part of the VIRGO2 database).

    You can adjust the coverage and completeness thresholds here if needed.
                
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns(2)
    with col1:
        min_coverage = st.number_input(
            "Coverage", value=10, placeholder="Default is 10"
        )
    with col2:
        min_completeness = st.number_input(
            "Complteness", value=80, placeholder="Default is 80"
        )

    st.write("MAGs with a minimum coverage of", min_coverage, "X and a completeness greather than", min_completeness, "% will be considered")

    virgo2_inventory = virgo2_metadata_all[(virgo2_metadata_all["Coverage"] > min_coverage) & (virgo2_metadata_all["Complete"] > min_completeness)]

    # st.subheader("Proportion of BGC identification - all MAGs", divider='grey')
    col1, col2 = st.columns(spec=[0.3,0.7])
    with col1:
        st.subheader("BGC identification", divider='grey')
        status_counts = antismash_status.loc[lambda df : df["MAG"].isin(virgo2_inventory["MAG"])]['status'].value_counts().reset_index()
        status_counts.columns = ['status', 'count']
        status_counts["status"] = status_counts["status"].replace({0:"No BGC", 1:">1 BGC"})
        status_colors = {"No BGC": zero_color, ">1 BGC": positive_color}
        fig = px.bar(status_counts,x="status",y="count",color="status",color_discrete_map=status_colors)
        fig.update_layout(yaxis_title="Number of MAGs")
        st.plotly_chart(fig)

    with col2:

        st.subheader("MAGs inventory", divider='grey')
        st.dataframe(pd.merge(virgo2_inventory, antismash_status, on='MAG', how='left'))

    st.header("Biosynthetic Gene Clusters - per taxa", divider='grey')

    st.markdown(f"""
    <div class="justified-text">
                
    Use the <b>Taxa selection</b> menu below to explore BGC detection for a specific taxon.

    The visualizations will include:
    <ul>
        <li>MAG detection rate</li>
        <li>Distribution of the number of BGCs per MAG</li>
        <li>Coverage of MAGs belonging to the selected taxon</li>
    </ul>
                
    </div>
    """, unsafe_allow_html=True)

    taxa_selection = st.selectbox("Taxa selection", sorted(virgo2_inventory['FinalTaxonomy'].unique()), index=sorted(virgo2_inventory['FinalTaxonomy'].unique()).index("Lactobacillus_crispatus"))

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("BGC detection rate")


        filtered_data = stack_antismash_status[(stack_antismash_status['FinalTaxonomy'] == taxa_selection) & (stack_antismash_status["MAG"].isin(virgo2_inventory["MAG"]))]
        status_counts = filtered_data.groupby('status').size().reset_index(name='count')
        status_counts["status"] = status_counts["status"].replace({0:"No BGC", 1:">1 BGC"})
        status_colors = {"No BGC": zero_color, ">1 BGC": taxa_color.get(taxa_selection, "#8c8c8c")}
        fig = px.bar(status_counts,x="status",y="count",color="status",text="count",color_discrete_map=status_colors)
        fig.update_layout(xaxis_title="")
        st.plotly_chart(fig, width='stretch')
    
    with col2:
        st.subheader("Number of BGC per MAG")
        # Ensure taxa_selection is treated as a list
        top_taxa = [taxa_selection]
        df = region_df.loc[lambda d: d["FinalTaxonomy"].isin(top_taxa)].groupby(["FinalTaxonomy", "MAG"])["GBK"].count().reset_index()
        df = df.merge(df.groupby("FinalTaxonomy", observed=False)["MAG"].nunique().reset_index().rename(columns={"MAG": "N_MAG"}))
        df = df[df["MAG"].isin(virgo2_inventory["MAG"])].reset_index(drop=True)

        # Color mapping for Plotly
        palette = {k: taxa_color.get(k, "#8c8c8c") for k in top_taxa}

        # Rename special case
        UBA629_label = "C. L. vaginae"
        df["FinalTaxonomy"] = df["FinalTaxonomy"].replace({"UBA629_sp005465875": UBA629_label})
        if "UBA629_sp005465875" in palette:
            palette[UBA629_label] = palette.pop("UBA629_sp005465875")

        # Plotly violin
        fig = px.violin(
            data_frame=df,
            x="FinalTaxonomy",
            y="GBK",
            color="FinalTaxonomy",
            color_discrete_map=palette,
            box=True,  # optional: add boxplot inside
        )

        fig.update_layout(
            xaxis_title="",
            yaxis_title="Number of BGCs per MAG",
            showlegend=False
        )

        st.plotly_chart(fig, width='stretch')

    with col3:
        st.subheader("MAG coverage")
        MIN_COVERAGE = 0
        coverage_taxa_plot = pd.merge(region_df, COVERAGE[COVERAGE["Coverage"] > MIN_COVERAGE], on="MAG", how="inner")
        coverage_taxa_plot = coverage_taxa_plot.loc[lambda df : df["FinalTaxonomy"].isin(top_taxa)]
        # coverage_taxa_plot["FinalTaxonomy"] = pd.Categorical(coverage_taxa_plot["FinalTaxonomy"],categories=species_to_plot,ordered=True)

        # sort by this categorical column to follow the given order
        coverage_taxa_plot = coverage_taxa_plot.sort_values(by="FinalTaxonomy")

        coverage_taxa_plot["FinalTaxonomy"] = coverage_taxa_plot["FinalTaxonomy"].replace({"UBA629_sp005465875":UBA629_label})
        coverage_taxa_plot = coverage_taxa_plot[coverage_taxa_plot["MAG"].isin(virgo2_inventory["MAG"])]

        # fig = px.box(coverage_taxa_plot, x="FinalTaxonomy", y="Coverage", color="FinalTaxonomy", log_y=True, color_discrete_map=palette, points="all")
        # fig.update_layout(
        #     xaxis_title="",
        #     yaxis_title="log(Coverage)",
        #     showlegend=False
        # )
        
        coverage_taxa_plot["log_Coverage"] = np.log10(coverage_taxa_plot["Coverage"] + 1)  # +1 to avoid log(0)
        fig = px.box(
            coverage_taxa_plot, 
            x="FinalTaxonomy", 
            y="log_Coverage", 
            color="FinalTaxonomy", 
            color_discrete_map=palette, 
            points="all"
        )
        fig.update_layout(
            xaxis_title="",
            yaxis_title="log10(Coverage)",
            showlegend=False
        )

        st.plotly_chart(fig, width='stretch')

# Run the page function
if __name__ == "__main__":
    page()


# the_cols = pd.DataFrame({
#     "species": [
#         "G. leopoldii", "G. piotii", "G. sp003585735", "G. sp003585845", "G. spNov1",
#         "G. spNov2", "G. vaginalis A", "G. vaginalis C", "G. vaginalis D", "G. vaginalis E",
#         "G. vaginalis F", "G. vaginalis H", "G. swidsinskii 1", "G. swidsinskii", 
#         "G. vaginalis", "Gardnerella"
#     ],
#     "color": [
#         "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
#         "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#1E3602",
#         "#E6AB02", "#666666"
#     ]
# })

# Functions for displaying data
# def display_antismash_status_pie(antismash_status):
#     status_counts = antismash_status['status'].value_counts().reset_index()
#     status_counts.columns = ['status', 'count']
#     color_map = {1: positive_color, 0: zero_color}
#     status_counts['color'] = status_counts['status'].map(color_map)
    
#     fig = px.pie(
#         status_counts, 
#         values='count', 
#         names='status', 
#         color='status',
#         title='Proportion of BGC identification (MAG)',
#         color_discrete_map=color_map
#     )
#     st.plotly_chart(fig)

# def display_numerical_feature_comparison(mag_inventory, antismash_status):
#     numerical_columns = [col for col in mag_inventory.select_dtypes(include=['float64', 'int64']).columns 
#                          if col not in ['Timepoint', 'FilterContam', 'red_value', 'warnings']]
    
#     n_cols = 4
#     n_rows = (len(numerical_columns) + n_cols - 1) // n_cols
#     fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=numerical_columns)

#     merged_data = mag_inventory.merge(antismash_status, on="MAG", how="left")

#     for i, col in enumerate(numerical_columns):
#         row, col_pos = divmod(i, n_cols)
#         for status in merged_data['status'].unique():
#             filtered_data = merged_data[merged_data['status'] == status]
#             fig.add_trace(
#                 go.Box(y=filtered_data[col], name=f"Status {status}", boxmean=True, marker_color=color_map[status]),
#                 row=row + 1, col=col_pos + 1
#             )
#         fig.update_yaxes(title_text="", row=row + 1, col=col_pos + 1)

#     fig.update_layout(height=500 * n_rows, width=1500, showlegend=False)
#     st.plotly_chart(fig)

# def plot_mean_sequence_length(mean_length_data):

#     fig = go.Figure()

#     for status in mean_length_data['status'].unique():
#         filtered_data = mean_length_data[mean_length_data['status'] == status]
#         fig.add_trace(
#             go.Histogram(
#                 x=filtered_data['length'],
#                 name=f"Status {status}",
#                 marker_color=color_map[status],
#                 nbinsx=200
#             )
#         )

#     fig.update_layout(
#         title="Mean Contig Length per MAG",
#         height=500, width=600,
#         xaxis_title="Contig Length",
#         yaxis_title="Frequency",
#         barmode='overlay',
#         xaxis=dict(range=[-20000, 1000000]),
#         showlegend=True
#     )
#     st.plotly_chart(fig)

# def plot_number_of_sequences(count_length_data):

#     fig = go.Figure()

#     for status in count_length_data['status'].unique():
#         filtered_data = count_length_data[count_length_data['status'] == status]
#         fig.add_trace(
#             go.Histogram(
#                 x=filtered_data['length'],
#                 name=f"Status {status}",
#                 marker_color=color_map[status],
#                 nbinsx=200
#             )
#         )

#     fig.update_layout(
#         title="Number of contigs per MAG",
#         height=500, width=600,
#         xaxis_title="Number of contigs",
#         yaxis_title="Frequency",
#         barmode='overlay',
#         xaxis=dict(range=[-1000, 10000]),
#         showlegend=True
#     )
#     st.plotly_chart(fig)

# def display_taxa_processed(taxa_filter=None):
#     if taxa_filter is not None :
#         to_plot = status_counts_long[status_counts_long['FinalTaxonomy'].str.contains(taxa_filter, case=False, na=False)] 
#     else :
#         to_plot = status_counts_long

#     color_map = {1: "blue", 0: "red"}
#     fig = px.bar(
#         to_plot, 
#         x='count', 
#         y='FinalTaxonomy', 
#         color='status', 
#         title="", 
#         labels={"FinalTaxonomy": "Final Taxonomy", "count": "Count", "status": "Status"}, 
#         height=1000,
#         text_auto=True,
#         color_discrete_map=color_map
#     )

#     st.plotly_chart(fig)