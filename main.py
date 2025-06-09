# main.py
import streamlit as st
import pandas as pd

st.set_page_config(layout='wide')

from pages_content import home, quality, virgo2_compounds, taxa_comparison

# Set up sidebar navigation with "Home" as the default page
st.sidebar.title("Navigation")
# page = st.sidebar.radio("Content", ["Home", "Quality", "Region summary", "VIRGO2 compounds"], index=0, label_visibility='hidden')
page = st.sidebar.radio("Content", ["Home", "Quality", "Taxa comparison"], index=0, label_visibility='hidden') #"VIRGO2 compounds",
# page = st.sidebar.radio("Content", ["Home", "Quality"], index=0, label_visibility='hidden') #"VIRGO2 compounds",



# Display the selected page
if page == "Home":
    home.page()
elif page == "Quality":
    quality.page()
# elif page == "Region summary":
#     region_summary_w_mibig.page()
# elif page == "VIRGO2 compounds":
#     virgo2_compounds.page()
elif page == "Taxa comparison":
    taxa_comparison.page()  
# elif page == "Compare species":
#     species_compare.page()