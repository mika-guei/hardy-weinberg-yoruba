import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Hardy-Weinberg - SNP Yoruba", layout="wide")

st.title("🧬 Análise de Frequências Genotípicas - Hardy-Weinberg")
st.write("Visualização interativa dos dados do cromossomo 21 (Yoruba, 1000 Genomes Project).")

# --- Upload do arquivo ---
uploaded_file = st.file_uploader("Envie seu arquivo .txt com os SNPs", type=["txt", "csv"])

if uploaded_file is not None:
    # --- Leitura do arquivo ---
    dados = pd.read_csv(uploaded_file, delim_whitespace=True)
    
    st.subheader("📄 Primeiras linhas dos dados:")
    st.dataframe(dados.head())

    n_individuos = 108
    n_alelos = n_individuos * 2

    # --- Cálculos ---
    dados["p_ref"] = (2 * dados["ref.ref"] + dados["ref.alt"]) / n_alelos
    dados["q_alt"] = 1 - dados["p_ref"]

    dados["freq_refref_obs"] = dados["ref.ref"] / n_individuos
    dados["freq_refalt_obs"] = dados["ref.alt"] / n_individuos
    dados["freq_altalt_obs"] = dados["alt.alt"] / n_individuos

    dados["freq_refref_exp"] = (1 - dados["q_alt"])**2
    dados["freq_refalt_exp"] = 2 * dados["q_alt"] * (1 - dados["q_alt"])
    dados["freq_altalt_exp"] = dados["q_alt"]**2

    # --- Filtro opcional de SNPs ---
    num_snps = st.slider("Quantos SNPs exibir no gráfico?", 100, len(dados), 1000, step=100)

    # --- Gráfico ---
    q = np.linspace(0, 1, 100)
    fig, ax = plt.subplots(figsize=(9, 7))

    ax.plot(q, (1-q)**2, label="ref.ref esperado", color="blue")
    ax.plot(q, 2*q*(1-q), label="ref.alt esperado", color="green")
    ax.plot(q, q**2, label="alt.alt esperado", color="red")

    # Pontos observados
    subset = dados.sample(num_snps)
    ax.scatter(subset["q_alt"], subset["freq_refref_obs"], color="blue", alpha=0.5, s=12)
    ax.scatter(subset["q_alt"], subset["freq_refalt_obs"], color="green", alpha=0.5, s=12)
    ax.scatter(subset["q_alt"], subset["freq_altalt_obs"], color="red", alpha=0.5, s=12)

    ax.set_xlabel("Frequência alélica (q - alelo alternativo)")
    ax.set_ylabel("Frequência genotípica")
    ax.set_title("Frequências genotípicas observadas vs esperadas (Hardy-Weinberg)\nPopulação Yoruba - Cromossomo 21")
    ax.legend()
    ax.grid(alpha=0.3)

    st.pyplot(fig)

    # --- Estatísticas simples ---
    st.subheader("📊 Estatísticas de exemplo:")
    st.write("Distribuição das frequências alélicas (q):")
    st.bar_chart(dados["q_alt"])

else:
    st.info("👆 Envie o arquivo `.txt` para começar a análise.")
