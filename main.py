from Bio import Entrez
import datetime
import requests
import os

Entrez.email = os.getenv("NCBI_EMAIL")
TERMO = '("pediatric intensive care" OR "neonatal intensive care") AND ("critically ill" OR sepsis)'
HOJE = datetime.date.today()
ONTEM = HOJE - datetime.timedelta(days=1)

BOT_TOKEN = os.getenv("TELEGRAM_BOT_TOKEN")
CHAT_ID = os.getenv("TELEGRAM_CHAT_ID")

def buscar_artigos():
    handle = Entrez.esearch(
        db="pubmed",
        term=TERMO,
        datetype="pdat",
        mindate=ONTEM.strftime("%Y/%m/%d"),
        maxdate=HOJE.strftime("%Y/%m/%d"),
        retmax=10,
        sort="pub+date"
    )
    record = Entrez.read(handle)
    return record.get("IdList", [])

def enviar_telegram(msg):
    url = f"https://api.telegram.org/bot{BOT_TOKEN}/sendMessage"
    r = requests.post(url, json={"chat_id": CHAT_ID, "text": msg, "disable_web_page_preview": True}, timeout=30)
    r.raise_for_status()

if __name__ == "__main__":
    artigos = buscar_artigos()
    if artigos:
        texto = f"🧠 Novos artigos ({len(artigos)}):\n" + "\n".join([f"https://pubmed.ncbi.nlm.nih.gov/{a}/" for a in artigos])
    else:
        texto = "Nenhum artigo novo nas últimas 24h."
    enviar_telegram(texto)
    print("✅ Mensagem enviada pelo Telegram!")
