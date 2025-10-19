import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Hardy-Weinberg - Comparação braços cromossomo 21", layout="wide")

st.title("🧬 Comparação Hardy-Weinberg - Braço curto vs Braço longo (Yoruba)")
st.write("Analise comparativa das frequências genotípicas observadas vs esperadas em duas regiões do cromossomo 21.")

# --- Upload dos dois arquivos ---
col_up1, col_up2 = st.columns(2)
with col_up1:
    file_short = st.file_uploader("📂 Envie o arquivo do **braço curto**", type=["txt", "csv"], key="short")
with col_up2:
    file_long = st.file_uploader("📂 Envie o arquivo do **braço longo**", type=["txt", "csv"], key="long")

def processar_dados(arquivo):
    dados = pd.read_csv(arquivo, delim_whitespace=True)
    n_individuos = 108
    n_alelos = n_individuos * 2

    dados["ref"] = (2 * dados["ref.ref"] + dados["ref.alt"]) / n_alelos
    dados["alt"] = 1 - dados["ref"]

    dados["freq_refref_obs"] = dados["ref.ref"] / n_individuos
    dados["freq_refalt_obs"] = dados["ref.alt"] / n_individuos
    dados["freq_altalt_obs"] = dados["alt.alt"] / n_individuos

    dados["freq_refref_exp"] = (1 - dados["alt"])**2
    dados["freq_refalt_exp"] = 2 * dados["alt"] * (1 - dados["alt"])
    dados["freq_altalt_exp"] = dados["alt"]**2

    return dados

# --- Processar e plotar ---
if file_short is not None and file_long is not None:
    dados_curto = processar_dados(file_short)
    dados_longo = processar_dados(file_long)

    num_snps = st.slider("Quantos SNPs exibir em cada gráfico?", 100, len(dados_curto), 1000, step=100)

    q = np.linspace(0, 1, 100)

    # --- Gráfico 1: Braço curto ---
    fig1, ax1 = plt.subplots(figsize=(8,6))
    ax1.plot(q, (1-q)**2, label="ref.ref esperado", color="blue")
    ax1.plot(q, 2*q*(1-q), label="ref.alt esperado", color="green")
    ax1.plot(q, q**2, label="alt.alt esperado", color="red")

    subset_c = dados_curto.sample(num_snps)
    ax1.scatter(subset_c["alt"], subset_c["freq_refref_obs"], color="blue", alpha=0.5, s=12)
    ax1.scatter(subset_c["alt"], subset_c["freq_refalt_obs"], color="green", alpha=0.5, s=12)
    ax1.scatter(subset_c["alt"], subset_c["freq_altalt_obs"], color="red", alpha=0.5, s=12)
    ax1.set_title("Braço curto do cromossomo 21")
    ax1.set_xlabel("Frequência alélica (alt - alelo alternativo)")
    ax1.set_ylabel("Frequência genotípica")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # --- Gráfico 2: Braço longo ---
    fig2, ax2 = plt.subplots(figsize=(8,6))
    ax2.plot(q, (1-q)**2, label="ref.ref esperado", color="blue")
    ax2.plot(q, 2*q*(1-q), label="ref.alt esperado", color="green")
    ax2.plot(q, q**2, label="alt.alt esperado", color="red")

    subset_l = dados_longo.sample(num_snps)
    ax2.scatter(subset_l["alt"], subset_l["freq_refref_obs"], color="blue", alpha=0.5, s=12)
    ax2.scatter(subset_l["alt"], subset_l["freq_refalt_obs"], color="green", alpha=0.5, s=12)
    ax2.scatter(subset_l["alt"], subset_l["freq_altalt_obs"], color="red", alpha=0.5, s=12)
    ax2.set_title("Braço longo do cromossomo 21")
    ax2.set_xlabel("Frequência alélica (alt - alelo alternativo)")
    ax2.legend()
    ax2.grid(alpha=0.3)

    # --- Mostrar lado a lado ---
    col1, col2 = st.columns(2)
    with col1:
        st.pyplot(fig1)
    with col2:
        st.pyplot(fig2)

else:
    st.info("👆 Envie os dois arquivos (.txt) para gerar os gráficos comparativos.")
