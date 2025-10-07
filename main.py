from Bio import Entrez
import datetime
import smtplib
from email.mime.text import MIMEText
import os

# Configurações do PubMed
Entrez.email = os.getenv("NCBI_EMAIL")  # variável de ambiente no GitHub
TERMO = '("pediatric intensive care" OR "neonatal intensive care") AND ("critically ill" OR sepsis)'
HOJE = datetime.date.today()
ONTEM = HOJE - datetime.timedelta(days=1)

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
    return record["IdList"]

def enviar_email(mensagem):
    remetente = os.getenv("EMAIL_FROM")
    senha = os.getenv("EMAIL_PASSWORD")
    destinatario = os.getenv("EMAIL_TO")

    msg = MIMEText(mensagem)
    msg["Subject"] = "🔬 Novos artigos PubMed - UTI Pediátrica"
    msg["From"] = remetente
    msg["To"] = destinatario

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(remetente, senha)
        server.send_message(msg)

if __name__ == "__main__":
    artigos = buscar_artigos()
    if artigos:
        texto = f"🧠 Novos artigos encontrados ({len(artigos)}):\n" + "\n".join([f"https://pubmed.ncbi.nlm.nih.gov/{a}/" for a in artigos])
    else:
        texto = "Nenhum artigo novo nas últimas 24h."
    enviar_email(texto)
    print("✅ E-mail enviado com sucesso!")
