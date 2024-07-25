




#!/bin/bash

# Verifica che sia stato passato un file come argomento
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 file_di_input"
    exit 1
fi

# Variabili per tenere traccia della somma totale e del conteggio
somma_totale=0
conteggio=0

# Legge il file riga per riga
while IFS=, read -r id tempo
do
    # Aggiunge il valore di tempo alla somma totale
    somma_totale=$(echo "$somma_totale + $tempo" | bc)
    # Incrementa il conteggio
    conteggio=$((conteggio + 1))
done < "$1"

# Calcola la media
if [ $conteggio -ne 0 ]; then

    media=$(echo "scale=4; $somma_totale / $conteggio" | bc)
    echo "AVG: $media"
else
    echo "Nessun dato da processare"
fi























