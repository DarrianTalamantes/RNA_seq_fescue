# Purpose: Find all the meanings behind my GO terms


import requests

def fetch_go_term_description(go_term):
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_term}"
    response = requests.get(url, headers={"Accept": "application/json"})
    if response.status_code == 200:
        data = response.json()
        return data['results'][0]['definition']['text']
    else:
        return "Description not found"

go_terms = ["GO:0008150", "GO:0005576", "GO:0005515"]  # Replace with your GO terms

for go_term in go_terms:
    description = fetch_go_term_description(go_term)
    print(f"{go_term}: {description}")