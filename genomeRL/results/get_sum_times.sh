#!/bin/bash

# Verifica che sia stato passato un file come argomento
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 file_di_input"
    exit 1
fi

# Variabile per tenere traccia della somma totale
somma_totale=0

# Legge il file riga per riga
while IFS=, read -r id tempo
do
    # Aggiunge il valore di tempo alla somma totale
    somma_totale=$(echo "$somma_totale + $tempo" | bc)
done < "$1"

# Stampa la somma totale
echo "Somma totale dei tempi: $somma_totale"

