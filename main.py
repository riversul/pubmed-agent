from Bio import Entrez
import datetime

# Configurações
Entrez.email = "vitorgabriel03@hotmail.com"  # coloque seu email (obrigatório pela API)
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
    ids = record["IdList"]
    return ids

if __name__ == "__main__":
    artigos = buscar_artigos()
    if artigos:
        print(f"🧠 Novos artigos encontrados ({len(artigos)}): {', '.join(artigos)}")
    else:
        print("Nenhum artigo novo nas últimas 24h.")
