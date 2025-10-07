from Bio import Entrez
import datetime
import requests
import os

Entrez.email = os.getenv("NCBI_EMAIL")
BOT_TOKEN = os.getenv("TELEGRAM_BOT_TOKEN")
CHAT_ID = os.getenv("TELEGRAM_CHAT_ID")

# ----- CONFIGURAÇÃO DA BUSCA -----
TERMO = (
    '("Intensive Care Units, Pediatric"[MeSH Terms] OR "Intensive Care Units, Neonatal"[MeSH Terms]) '
    'AND ("Critical Illness"[MeSH Terms] OR "Sepsis"[MeSH Terms] OR "Multiple Organ Failure"[MeSH Terms])'
)
HOJE = datetime.date.today()
ONTEM = HOJE - datetime.timedelta(days=1)

def buscar_pmids():
    handle = Entrez.esearch(
        db="pubmed",
        term=TERMO,
        datetype="pdat",
        mindate=ONTEM.strftime("%Y/%m/%d"),
        maxdate=HOJE.strftime("%Y/%m/%d"),
        retmax=5,  # quantos artigos buscar por dia
        sort="pub+date"
    )
    record = Entrez.read(handle)
    return record.get("IdList", [])

def detalhes_artigo(pmid):
    """Busca título e abstract do artigo"""
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    artigo = records["PubmedArticle"][0]
    titulo = artigo["MedlineCitation"]["Article"]["ArticleTitle"]
    resumo = ""
    if "Abstract" in artigo["MedlineCitation"]["Article"]:
        parts = artigo["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
        resumo = " ".join(parts)
    return titulo, resumo

def traduzir_para_pt(texto):
    # usa API simples de tradução gratuita (MyMemory)
    try:
        resp = requests.get(
            "https://api.mymemory.translated.net/get",
            params={"q": texto, "langpair": "en|pt"}
        )
        data = resp.json()
        return data["responseData"]["translatedText"]
    except:
        return texto  # se falhar, devolve original

def enviar_telegram(msg):
    url = f"https://api.telegram.org/bot{BOT_TOKEN}/sendMessage"
    r = requests.post(url, json={"chat_id": CHAT_ID, "text": msg, "disable_web_page_preview": True}, timeout=30)
    r.raise_for_status()

if __name__ == "__main__":
    pmids = buscar_pmids()
    if not pmids:
        enviar_telegram("Nenhum artigo novo nas últimas 24h.")
    else:
        mensagens = []
        for pmid in pmids:
            titulo, resumo = detalhes_artigo(pmid)
            resumo_pt = traduzir_para_pt(resumo) if resumo else "Sem resumo disponível."
            mensagens.append(f"**{titulo}**\n{resumo_pt}\n🔗 https://pubmed.ncbi.nlm.nih.gov/{pmid}/\n")
        texto = "🧠 *Novos artigos encontrados:*\n\n" + "\n\n".join(mensagens)
        # O Telegram usa Markdown, então vamos marcar em Markdown simples:
        enviar_telegram(texto)
    print("✅ Mensagem enviada pelo Telegram!")
