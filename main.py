from Bio import Entrez
import datetime
import requests
import os
import gspread
from oauth2client.service_account import ServiceAccountCredentials

Entrez.email = os.getenv("NCBI_EMAIL")
BOT_TOKEN = os.getenv("TELEGRAM_BOT_TOKEN")

# ---- CONFIGURAÇÃO GOOGLE SHEETS ----
SHEET_NAME = "Assinantes_PubMed"

def get_assinantes():
    creds_json = os.getenv("GOOGLE_CREDENTIALS")
    if not creds_json:
        raise RuntimeError("Credenciais do Google não encontradas.")
    import json
    creds_dict = json.loads(creds_json)
    scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_dict(creds_dict, scope)
    gc = gspread.authorize(creds)
    sheet = gc.open(SHEET_NAME).sheet1
    rows = sheet.get_all_records()
    return [str(r["telegram_id"]) for r in rows if r.get("telegram_id")]

# ---- CONFIGURAÇÃO PUBMED ----
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
        retmax=5,
        sort="pub+date"
    )
    record = Entrez.read(handle)
    return record.get("IdList", [])

def detalhes_artigo(pmid):
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="xml")
    records = Entrez.read(handle)
    artigo = records["PubmedArticle"][0]["MedlineCitation"]["Article"]

    titulo = artigo.get("ArticleTitle", "Sem título")
    resumo = ""
    if "Abstract" in artigo and "AbstractText" in artigo["Abstract"]:
        partes = artigo["Abstract"]["AbstractText"]
        resumo = " ".join(str(p) for p in partes)
    return titulo, resumo

def traduzir_para_pt(texto):
    if not texto:
        return ""
    try:
        resp = requests.get(
            "https://api.mymemory.translated.net/get",
            params={"q": texto, "langpair": "en|pt"},
            timeout=30
        )
        return resp.json().get("responseData", {}).get("translatedText", texto)
    except:
        return texto

def enviar_telegram(chat_id, msg):
    url = f"https://api.telegram.org/bot{BOT_TOKEN}/sendMessage"
    requests.post(
        url,
        json={"chat_id": chat_id, "text": msg, "parse_mode": "Markdown", "disable_web_page_preview": True},
        timeout=30
    ).raise_for_status()

if __name__ == "__main__":
    assinantes = get_assinantes()
    pmids = buscar_pmids()
    if not pmids:
        for user in assinantes:
            enviar_telegram(user, "Nenhum artigo novo nas últimas 24h.")
    else:
        mensagens = []
        for pmid in pmids:
            titulo, resumo = detalhes_artigo(pmid)
            if resumo:
                resumo_pt = traduzir_para_pt(resumo)
            else:
                resumo_pt = "_Resumo não disponível._"
            mensagens.append(f"*{titulo}*\n{resumo_pt}\n🔗 https://pubmed.ncbi.nlm.nih.gov/{pmid}/\n")
        texto = "🧠 *Novos artigos encontrados:*\n\n" + "\n\n".join(mensagens)
        for user in assinantes:
            enviar_telegram(user, texto)
    print("✅ Mensagens enviadas!")
