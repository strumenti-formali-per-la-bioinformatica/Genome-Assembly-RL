# Reinforcement Learning Genome Assembly

## Problema 
L'assemblaggio di un genoma è uno dei compiti più importanti e computazionalmente complessi in bioinformatica, che consiste nel ricostruire un genoma a partire da porzioni più piccole chiamate reads, che corrispondono rispettivamente a sequenze di nucleotidi (A, C, T, G) per DNA ed RNA ed a sequenze di amminoacidi nel caso delle proteine. 
Ricostruire un genoma è estremamente utile: i sequenziatori  non ci consentono di "leggere" un genoma interamente per via della loro lunghezza, anche nel caso di microrganismi. Le reads utilizzate vengono prodotte tramite lo shotgun sequencing: il genoma viene diviso randomicamente in segmenti più piccoli sottoforma di sequenze di testo. Compito dell'assemblatore sarà ricongiungere i frammenti nel corretto ordine analizzando la loro sovrapposizione o overlap.
L'assemblaggio di un genoma è perciò un problema non ancora risolto, i cui risultati possono variare dipendentemente dall'assemblatore che si utilizza e che richiede ancora il supporto di esperti del settore per la configurazione. 
Sono nati così approcci alternativi al problema, in particolare per questo progetto viene considerato il Reinforcement Learning (RL). 

## Obiettivo
Questo lavoro si concentra sull'utilizzo di Reinforcement Learning (RL) per l'assemblamento di un genoma de novo. Nel contesto del RL, un agente effettua osservazioni in un ambiente e prende decisioni ricevendo ricompense o reward positive o negative. L'obiettivo ultimo è quello di costruire un agente che, avendo a disposizione una serie di reads genomiche, riesca autonomamente a ricomporre il genoma completo. A partire da un modello che combina RL e un algoritmo genetico, si propone un approccio RL end-to-end, dove viene utilizzato non solo per la composizione di un pool di cromosomi per un approccio evolutionary-based, ma anche nelle fasi di mutazione e crossing-over. 

# Utilizzo

## Esperimenti
Il progetto richiede l'installazione di Python3 sul sistema d'utilizzo.

Per eseguire il genome assembly con RL su un singolo genoma, utilizzare il seguente comando:
```sh
python3 main.py episodes genome seed
```
Dove: 
- episodes indica il numero di episodi da adottare;
- genome specifica il genoma da scegliere nel dataset;
- seed specifica un eventuale seed di randomness, per utilizzare lo stesso seed in caso di esperimenti di confronto.

Per riprodurre il progetto su tutti i genomi del dataset, senza salvare risultati:
```sh
./reproduce.sh
```

Per riprodurre il progetto su tutti i genomi salvando i risultati in file per una successiva analisi utilizzare il comando:
```sh
./sC.sh 
```

I file contenenti i risultati avranno il seguente nome:
```sh
results_episodes_genome_seed.txt
```
Dove il cambio "episodes" varia a seconda degli episodi richiesti e i campi "genome" e "seed" variano a seconda dei genomi e semi casuali di riferimento.

L'output dei risultati include le stampe del programma, che possono essere analizzate per valutare le prestazioni e l'accuratezza dell'assemblaggio genomico.

## Risultati
Una volta ottenuti i risultati e trasferiti i file nella cartella "/genomeRL/results/" è possibile analizzarli attraverso i seguenti comandi:
Per ottenere il tempo di esecuzione totale:
```sh
./get_time.sh > tempi.txt
```

I risultati avranno il seguente formato:
```sh
genome,tempo (s)
```

Per ottenere il tempo totale dell'esecuzione, su tutti i genomi, utilizzare il comando seguente, passando in input il file dei tempi generato in precedenza:
```sh
./get_sum_times.sh tempi.txt
```

Per ottenere la Reward Measure, che misura la distanza in termini reward del risultato rispetto al genoma originale:
```sh
./get_pm.sh > pm.txt
```
I risultati avranno il seguente formato:
```sh
genome,reward_measure
```

Per ottenere tutti i semi, in caso di riproduzione dell'esperimento più volte:
```sh
./get_seeds.sh
```
I risultati sono strutturati nel formato seguente:
```sh
genome,seed
```





