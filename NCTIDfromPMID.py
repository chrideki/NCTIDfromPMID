#!/usr/bin/env python
# coding: utf-8

from Bio import Entrez
import pandas as pd
import time
import csv
import requests
import xml.etree.ElementTree as ET
import re
import streamlit as st

def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

def parse_nct_from_etree(article_elem):
    pmid = article_elem.findtext(".//PMID")
    nct_ids = []

    # 1. OtherID / SI field
    other_ids = article_elem.findall(".//OtherID")
    for oid in other_ids:
        if oid.text and "NCT" in oid.text:
            nct_ids.append(oid.text)

    # 2. DataBankList
    for databank_list in article_elem.findall(".//DataBankList"):
        for databank in databank_list.findall(".//DataBank"):
            db_name = databank.findtext("DataBankName")
            if db_name == "ClinicalTrials.gov":
                for acc in databank.findall(".//AccessionNumber"):
                    if acc.text:
                        nct_ids.append(acc.text)

    # 3. Abstract search
    if not nct_ids:
        abstract_texts = article_elem.findall(".//AbstractText")
        abstract = " ".join([at.text for at in abstract_texts if at.text])
        found = re.findall(r"NCT\d{8}", abstract)
        if found:
            nct_ids.extend(found)

    return pmid, list(set(nct_ids)) if nct_ids else None


def fetch_nct_ids_combined(pmids, batch_size=200, email="christian.goldoni@gmail.com"):
    Entrez.email = email
    results = {}

    for batch in chunk_list(pmids, batch_size):
        # Fetch raw XML
        handle = Entrez.efetch(db="pubmed", id=",".join(batch), retmode="xml")
        xml_data = handle.read()
        handle.close()

        # Try Entrez.read() first
        try:
            records = Entrez.read(xml_data)
            articles = records["PubmedArticle"]
            found_pmids = set()

            for article in articles:
                pmid = str(article["MedlineCitation"]["PMID"])
                found_pmids.add(pmid)
                nct_ids = []

                # 1. OtherID
                if "OtherID" in article["MedlineCitation"]:
                    other_ids = article["MedlineCitation"]["OtherID"]
                    if not isinstance(other_ids, list):
                        other_ids = [other_ids]
                    nct_ids = [oid for oid in other_ids if "NCT" in oid]

                # 2. DataBankList
                if not nct_ids and "Article" in article["MedlineCitation"]:
                    article_section = article["MedlineCitation"]["Article"]
                    if "DataBankList" in article_section:
                        data_bank_list = article_section["DataBankList"]
                        if isinstance(data_bank_list, list):
                            data_bank_list = data_bank_list[0]

                        if "DataBank" in data_bank_list:
                            databanks = data_bank_list["DataBank"]
                            if not isinstance(databanks, list):
                                databanks = [databanks]
                            for databank in databanks:
                                if databank["DataBankName"] == "ClinicalTrials.gov":
                                    acc_nums = databank["AccessionNumberList"]["AccessionNumber"]
                                    if not isinstance(acc_nums, list):
                                        acc_nums = [acc_nums]
                                    nct_ids.extend(acc_nums)

                # 3. Abstract
                if not nct_ids:
                    if "Abstract" in article["MedlineCitation"]["Article"]:
                        abstract = article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
                        if isinstance(abstract, list):
                            abstract = " ".join(str(a) for a in abstract)
                        else:
                            abstract = str(abstract)
                        nct_ids = re.findall(r"NCT\d{8}", abstract)

                results[pmid] = list(set(nct_ids)) if nct_ids else None

        except Exception:
            # If Entrez.read() fails, fallback to ElementTree parse
            articles = ET.fromstring(xml_data).findall(".//PubmedArticle")
            found_pmids = set()
            for art in articles:
                pmid, nct_ids = parse_nct_from_etree(art)
                found_pmids.add(pmid)
                results[pmid] = nct_ids

        # For any PMIDs missing in results, parse manually with ElementTree
        missing_pmids = [pmid for pmid in batch if pmid not in results]
        
        if missing_pmids:
            root = ET.fromstring(xml_data)
            for art in root.findall(".//PubmedArticle"):
                pmid = art.findtext(".//PMID")
                if pmid in missing_pmids:
                    pmid, nct_ids = parse_nct_from_etree(art)
                    results[pmid] = nct_ids

        # Ensure all PMIDs accounted for (add None if missing)
        for pmid in batch:
            if pmid not in results:
                results[pmid] = None

    return results


st.title("PMID to NCT ID Lookup")

uploaded_file = st.file_uploader("Upload Excel with PMIDs")

if uploaded_file:
    df = pd.read_excel(uploaded_file)
    pmid_list = df['PMID'].astype(str).tolist()

    results = fetch_nct_ids_combined(pmid_list, email="christian.goldoni@gmail.com")
    out_df = pd.DataFrame([(pmid, ", ".join(ids) if ids else None) for pmid, ids in results.items()],
                          columns=["PMID", "NCT_IDs"])
    st.dataframe(out_df)

    st.write("PMIDs in the original Excel:", len(pmid_list))
    st.write("Unique PMIDs in the original Excel:", len(results))
    st.write("NCTIDs found:", out_df['NCT_IDs'].notna().sum())
    
    st.download_button("Download Results", out_df.to_csv(index=False), file_name="nct_results.csv")
    






